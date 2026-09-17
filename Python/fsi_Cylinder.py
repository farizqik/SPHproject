#!/usr/bin/env python3
"""2D vortex shedding / one-degree-of-freedom vortex-induced vibration.

Install: python -m pip install numpy matplotlib pillow numba
Numba is optional but recommended for speed.
Spyder: enter %matplotlib qt in the IPython console, then press F5.
Live frames default to every 0.05 tU/D; --frame-dt 0.025 gives more frames.
Bilinear colour interpolation affects display only, not the numerical solution.
GIF frames default to every 0.10 tU/D and play at 25 fps. Use --gif-dt 0.05
for finer temporal sampling. Longer GIFs require more encoding time and memory.
Numba disk caching is disabled for compatibility with Spyder runfile.
Run:     python fsi_cylinder.py
Fixed:   python fsi_cylinder.py --fixed
Movie:   python fsi_cylinder.py --gif --duration 160
Headless: python fsi_cylinder.py --no-live
Quick check: python fsi_cylinder.py --duration 2 --no-live

D2Q9 TRT lattice Boltzmann + smooth volume penalization + Guo forcing.
A rigid cylinder translates in y, with x fixed. Mass and forces are per unit
span. It is coupled to
    m*y'' + c*y' + k*y = lift.
All internal units are lattice units (dx = dt = rho_ref = 1).
Reported time is t*U/D, displacement y/D, vorticity omega*D/U.
For a physical D_phys, U_phys: dx_phys=D_phys/D and
dt_phys=U*dx_phys/U_phys; physical nu=U_phys*D_phys/Re.
Re = U*D/nu; Ur = U/(fn*D); m* = m/(rho*pi*D**2/4).
Domain: 24D by 10D, cylinder at x=5D. Uniform velocity inlet,
zero-gradient outlet, periodic transverse boundaries (an array of cylinders
spaced 10D apart, approximating external flow). No physical top/bottom walls.

This is an educational, weakly coupled solver, not a benchmark-validated FSI
package. The diffuse surface changes effective diameter; refine D, enlarge
width/downstream extent, vary eta, and compare fixed-body Cd and Strouhal
before quantitative use. High-Re 3D turbulence and flexible deformation are
not modeled. Explicit FSI coupling is unsuitable for very low mass ratios.
The force includes an approximate fictitious interior-fluid momentum
correction; its derivative and moving-mask discretization need refinement.
The movie never prescribes a shedding frequency or sinusoidal lift.

Reference for the implemented fluid forcing:
Guo, Zheng & Shi (2002), Phys. Rev. E 65, 046308.
https://doi.org/10.1103/PhysRevE.65.046308
"""
from pathlib import Path
import argparse
import csv
import json
import time
import numpy as np

# Lattice directions: rest, E, N, W, S, NE, NW, SW, SE.
C = np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1],
              [1,1],[-1,1],[-1,-1],[1,-1]], dtype=int)
W = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36])
CX = C[:,0,None,None]
CY = C[:,1,None,None]
OPP = np.array([0,3,4,1,2,7,8,5,6])


def equilibrium(rho, ux, uy):
    cu = CX*ux + CY*uy
    return W[:,None,None]*rho*(1 + 3*cu + 4.5*cu**2
                                      - 1.5*(ux**2+uy**2))


def moments(f):
    rho = f.sum(axis=0)
    return rho, (f*CX).sum(axis=0), (f*CY).sum(axis=0)


def cylinder_mask(X, Y, xc, yc, diameter):
    # One lattice-cell transition centered on the nominal radius.
    distance = np.sqrt((X-xc)**2 + (Y-yc)**2)
    return 0.5*(1 - np.tanh((distance-diameter/2)/0.75))


def fluid_step(f, mask, cylinder_v, inlet_u, tau, eta):
    rho, jx, jy = moments(f)
    # Implicit solution of rho*u = j + F/2, F=beta*(u_s-u).
    beta = rho*mask/eta
    ux = jx/(rho + 0.5*beta)
    uy = (jy + 0.5*beta*cylinder_v)/(rho + 0.5*beta)
    fx = -beta*ux
    fy = beta*(cylinder_v-uy)
    raw_force = -np.array([fx.sum(), fy.sum()])
    # Momentum of the fluid occupying the smooth solid region.
    interior_momentum = np.array([(mask*rho*ux).sum(),
                                  (mask*rho*uy).sum()])
    cu = CX*ux + CY*uy
    source = W[:,None,None]*(3*((CX-ux)*fx+(CY-uy)*fy)
                                      + 9*cu*(CX*fx+CY*fy))
    # Two-relaxation-time collision: viscous even and stabilizing odd modes.
    # Magic parameter Lambda=(tau_plus-.5)*(tau_minus-.5)=3/16.
    tau_minus = 0.5 + (3/16)/(tau-0.5)
    delta = f-equilibrium(rho,ux,uy)
    even = 0.5*(delta+delta[OPP])
    odd = 0.5*(delta-delta[OPP])
    source_even = 0.5*(source+source[OPP])
    source_odd = 0.5*(source-source[OPP])
    post = (f-even/tau-odd/tau_minus
            +(1-0.5/tau)*source_even+(1-0.5/tau_minus)*source_odd)
    streamed = np.empty_like(f)
    for q, (ex,ey) in enumerate(C):
        # y is periodic; x wrap is overwritten by inlet/outlet reconstruction.
        streamed[q] = np.roll(post[q], (ey,ex), axis=(0,1))

    # Zou/He velocity inlet: reconstruct ONLY the populations entering at x=0.
    a = streamed[:,:,0]
    r = (a[0]+a[2]+a[4]+2*(a[3]+a[6]+a[7]))/(1-inlet_u)
    a[1] = a[3] + (2/3)*r*inlet_u
    a[5] = a[7] + (1/6)*r*inlet_u + 0.5*(a[4]-a[2])
    a[8] = a[6] + (1/6)*r*inlet_u + 0.5*(a[2]-a[4])
    # Simple open outlet. Only unknown left-going populations are extrapolated.
    for q in (3,6,7):
        streamed[q,:,-1] = streamed[q,:,-2]
    return streamed, raw_force, interior_momentum, rho, ux, uy


# Optional compiled version of the same update. No parallel reductions or
# fast-math approximations: useful for comparing against the NumPy fallback.
try:
    from numba import njit
except ImportError:
    njit = None


def compiled_fluid_step(f, mask, cylinder_v, inlet_u, tau, eta):
    _, ny, nx = f.shape
    streamed = np.empty_like(f)
    density = np.empty((ny,nx))
    ux_field = np.empty((ny,nx))
    uy_field = np.empty((ny,nx))
    eq = np.empty(9)
    source = np.empty(9)
    delta = np.empty(9)
    rawx = rawy = px = py = 0.0
    tau_minus = 0.5+(3/16)/(tau-0.5)
    for iy in range(ny):
        for ix in range(nx):
            rho = jx = jy = 0.0
            for q in range(9):
                rho += f[q,iy,ix]
                jx += C[q,0]*f[q,iy,ix]
                jy += C[q,1]*f[q,iy,ix]
            beta = rho*mask[iy,ix]/eta
            ux = jx/(rho+0.5*beta)
            uy = (jy+0.5*beta*cylinder_v)/(rho+0.5*beta)
            fx = -beta*ux
            fy = beta*(cylinder_v-uy)
            rawx -= fx
            rawy -= fy
            px += mask[iy,ix]*rho*ux
            py += mask[iy,ix]*rho*uy
            density[iy,ix] = rho
            ux_field[iy,ix] = ux
            uy_field[iy,ix] = uy
            for q in range(9):
                ex,ey = C[q,0],C[q,1]
                cu = ex*ux+ey*uy
                eq[q] = W[q]*rho*(1+3*cu+4.5*cu*cu-1.5*(ux*ux+uy*uy))
                delta[q] = f[q,iy,ix]-eq[q]
                source[q] = W[q]*(3*((ex-ux)*fx+(ey-uy)*fy)+9*cu*(ex*fx+ey*fy))
            for q in range(9):
                opposite = OPP[q]
                even = 0.5*(delta[q]+delta[opposite])
                odd = 0.5*(delta[q]-delta[opposite])
                se = 0.5*(source[q]+source[opposite])
                so = 0.5*(source[q]-source[opposite])
                post = f[q,iy,ix]-even/tau-odd/tau_minus+(1-0.5/tau)*se+(1-0.5/tau_minus)*so
                desty = (iy+C[q,1]) % ny
                destx = (ix+C[q,0]) % nx
                streamed[q,desty,destx] = post
    for iy in range(ny):
        r = (streamed[0,iy,0]+streamed[2,iy,0]+streamed[4,iy,0]
             +2*(streamed[3,iy,0]+streamed[6,iy,0]+streamed[7,iy,0]))/(1-inlet_u)
        streamed[1,iy,0] = streamed[3,iy,0]+(2/3)*r*inlet_u
        streamed[5,iy,0] = streamed[7,iy,0]+(1/6)*r*inlet_u+0.5*(streamed[4,iy,0]-streamed[2,iy,0])
        streamed[8,iy,0] = streamed[6,iy,0]+(1/6)*r*inlet_u+0.5*(streamed[2,iy,0]-streamed[4,iy,0])
        for q in (3,6,7):
            streamed[q,iy,nx-1] = streamed[q,iy,nx-2]
    return (streamed,np.array([rawx,rawy]),np.array([px,py]),
            density,ux_field,uy_field)


if njit is not None:
    compiled_fluid_step = njit(cache=False)(compiled_fluid_step)


def oscillator_step(y, v, force, mass, damping, stiffness):
    # RK4, dt=1, fluid force held constant over this structural substep.
    def rhs(a,b):
        return b, (force-damping*b-stiffness*a)/mass
    a1,b1 = rhs(y,v)
    a2,b2 = rhs(y+a1/2,v+b1/2)
    a3,b3 = rhs(y+a2/2,v+b2/2)
    a4,b4 = rhs(y+a3,v+b3)
    return y+(a1+2*a2+2*a3+a4)/6, v+(b1+2*b2+2*b3+b4)/6


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--numpy', action='store_true', help='Disable optional Numba acceleration')
    p.add_argument('--fixed', action='store_true', help='Lock cylinder position')
    p.add_argument('--diameter', type=int, default=16, help='Grid cells per D (refine to 24/32)')
    p.add_argument('--re', type=float, default=100)
    p.add_argument('--ur', type=float, default=5, help='Reduced velocity U/(fn D)')
    p.add_argument('--mass-ratio', type=float, default=10)
    p.add_argument('--zeta', type=float, default=0.01, help='Structural damping ratio')
    p.add_argument('--duration', type=float, default=120, help='Final nondimensional time tU/D')
    p.add_argument('--width', type=float, default=10, help='Domain height / D')
    p.add_argument('--length', type=float, default=24, help='Domain length / D')
    p.add_argument('--eta', type=float, default=0.5, help='Penalization relaxation time')
    p.add_argument('--no-live', action='store_true')
    p.add_argument('--gif', action='store_true', help='Save GIF (needs Pillow)')
    p.add_argument('--frame-dt', type=float, default=0.05,
                   help='Live frame interval in tU/D; smaller gives smoother motion')
    p.add_argument('--gif-dt', type=float, default=0.10,
                   help='GIF frame interval in tU/D')
    p.add_argument('--fps', type=int, default=25, help='GIF playback frames per second')
    p.add_argument('--output', default='fsi_results')
    a = p.parse_args()
    if (a.diameter < 12 or a.re <= 0 or a.ur <= 0 or a.mass_ratio < 2
            or a.zeta < 0 or a.duration <= 0 or a.width < 6
            or a.length < 16 or a.eta < 0.1
            or a.frame_dt <= 0 or a.gif_dt <= 0 or a.fps < 1):
        p.error('Require frame-dt/gif-dt>0, fps>=1; D>=12, Re/Ur/duration>0, m*>=2, zeta>=0, width>=6D, length>=16D, eta>=0.1')
    return a


def main():
    args = parse_args()
    args.gif = True       # Enable the existing frame-capture code
    args.gif_dt = 0.025   # Capture more closely spaced simulation frames
    args.fps = 60
    import matplotlib
    if args.no_live:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    D, U = args.diameter, 0.05
    nu = U*D/args.re
    tau = 0.5 + 3*nu
    if tau < 0.52:
        raise ValueError('Viscous tau < 0.52: increase diameter or reduce Re.')
    nx, ny = round(args.length*D), round(args.width*D)
    X,Y = np.meshgrid(np.arange(nx), np.arange(ny))
    xc, yc0 = 5*D, ny/2
    mass = args.mass_ratio*np.pi*D**2/4
    omega_n = 2*np.pi*U/(args.ur*D)
    stiffness = mass*omega_n**2
    damping = 2*args.zeta*mass*omega_n
    steps = int(np.ceil(args.duration*D/U))
    ramp_steps = 5*D/U
    # Small localized transverse disturbance breaks exact wake symmetry.
    ux = np.zeros((ny,nx))
    uy = 0.002*U*np.exp(-((X-(xc+2*D))/(2*D))**2
                           -((Y-yc0)/(0.7*D))**2)
    f = equilibrium(np.ones_like(ux),ux,uy)
    y, v = 0.01*D if not args.fixed else 0.0, 0.0
    previous_p = None
    out = Path(args.output)
    out.mkdir(parents=True,exist_ok=True)
    config = vars(args).copy()
    config.update(U=U,nu=nu,tau=tau,mass=mass,k=stiffness,c=damping,
                  nx=nx,ny=ny,dt=1,force_correction='raw + d(interior momentum)/dt')
    (out/'parameters.json').write_text(json.dumps(config,indent=2))
    print(f'Re={args.re:g}, D={D}, grid={nx}x{ny}, tau={tau:.4f}, steps={steps}')
    stepper = compiled_fluid_step if njit is not None and not args.numpy else fluid_step
    print('Backend:', 'Numba (first call compiles)' if stepper is compiled_fluid_step else 'NumPy (slower)')
    # Collect diagnostics at ~0.05 convective-time intervals.
    sample_stride = max(1,round(0.05*D/U))
    frame_stride = max(1,round(args.frame_dt*D/U))
    gif_stride = max(1,round(args.gif_dt*D/U))
    records = []
    # Keep compressed frames on disk while solving, rather than RGB images in RAM.
    import tempfile
    frame_store = tempfile.TemporaryDirectory(prefix='fsi_frames_')
    frame_paths = []
    fig, (ax, trace) = plt.subplots(2,1,figsize=(11,6),
                  gridspec_kw={'height_ratios':[3,1]}, constrained_layout=True)
    image = ax.imshow(np.zeros_like(ux),origin='lower',
                    extent=(-0.5/D,(nx-0.5)/D,-0.5/D,(ny-0.5)/D),
                    cmap='RdBu_r',vmin=-2,vmax=2,aspect='equal',
                    interpolation='bilinear',resample=True)
    circle = Circle((xc/D,(yc0+y)/D),0.5,color='black')
    ax.add_patch(circle)
    ax.set(xlabel='x/D',ylabel='y/D')
    fig.colorbar(image,ax=ax,label='Vorticity omega D/U')
    curve, = trace.plot([],[],color='black')
    trace.set(xlabel='t U/D',ylabel='Displacement / D',xlim=(0,args.duration),ylim=(-1,1))
    ax.set_title(f'{"Fixed cylinder" if args.fixed else "Coupled cylinder FSI"} | Re={args.re:g}')
    status = ax.text(0.02,0.95,'Starting...',transform=ax.transAxes,
                     va='top',bbox=dict(facecolor='white',alpha=0.85,edgecolor='none'))
    # Blitting redraws only moving artists; the axes/colorbar stay in the background.
    artists = [image,circle,curve,status]
    use_blit = not args.no_live and fig.canvas.supports_blit
    background = [None]
    def cache_background(event):
        if use_blit:
            background[0] = fig.canvas.copy_from_bbox(fig.bbox)
    if use_blit:
        for artist in artists:
            artist.set_animated(True)
        fig.canvas.mpl_connect('draw_event',cache_background)
    if not args.no_live:
        plt.ion()
        plt.show()
    fig.canvas.draw()
    # Freeze layout after the initial draw to avoid repeated layout calculations.
    fig.set_layout_engine(None)
    started = time.perf_counter()
    last_progress = started
    last_events = started
    with (out/'history.csv').open('w',newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['tU_D','y_D','v_U','Cd','Cl','inlet_U_fraction','rho_min','rho_max'])
        for n in range(steps):
            # Cylinder state and forces refer to time n; advance structure last.
            mask = cylinder_mask(X,Y,xc,yc0+y,D)
            inlet = U*min(1.0,(n+1)/ramp_steps)
            f, raw, interior_p, rho, ux, uy = stepper(f,mask,v,inlet,tau,args.eta)
            force = raw.copy()
            if previous_p is not None:
                force += interior_p-previous_p  # lattice dt=1
            previous_p = interior_p
            if n % sample_stride == 0 or n == steps-1:
                if (not np.all(np.isfinite(f)) or rho.min()<0.8 or rho.max()>1.2
                        or np.max(np.hypot(ux,uy))>0.25):
                    raise RuntimeError(f'Flow unstable at step {n}; refine grid, increase mass/damping, or reduce Re.')
                if abs(y)>ny/2-2*D:
                    raise RuntimeError('Cylinder too close to periodic domain edge; enlarge --width.')
                row = [n*U/D,y/D,v/U,force[0]/(0.5*U**2*D),
                       force[1]/(0.5*U**2*D),inlet/U,rho.min(),rho.max()]
                records.append(row)
                writer.writerow(row)
            live_frame = not args.no_live and n % frame_stride == 0
            gif_frame = args.gif and n % gif_stride == 0
            if live_frame or gif_frame or n == steps-1:
                # Periodic central derivative in y; standard central in x.
                vort = (np.gradient(uy,axis=1)
                        -(np.roll(ux,-1,axis=0)-np.roll(ux,1,axis=0))/2)*D/U
                image.set_data(np.ma.masked_where(mask>0.5,vort))
                circle.center = (xc/D,(yc0+y)/D)
                status.set_text(f'tU/D = {n*U/D:.2f}    y/D = {y/D:+.3f}')
                data = np.asarray(records)
                curve.set_data(data[:,0],data[:,1])
                ymax = max(1,np.max(np.abs(data[:,1]))*1.2)
                if ymax > trace.get_ylim()[1]:
                    trace.set_ylim(-ymax,ymax)
                    background[0] = None
                if not args.no_live:
                    if use_blit:
                        if background[0] is None:
                            fig.canvas.draw()
                        fig.canvas.restore_region(background[0])
                        for artist in artists:
                            artist.axes.draw_artist(artist)
                        fig.canvas.blit(fig.bbox)
                    else:
                        fig.canvas.draw()
                    fig.canvas.flush_events()
                if args.gif and (gif_frame or n == steps-1):
                    # A complete draw includes every artist in the exported frame.
                    for artist in artists:
                        artist.set_animated(False)
                    fig.canvas.draw()
                    from PIL import Image
                    frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3])
                    path = Path(frame_store.name)/f'{len(frame_paths):06d}.png'
                    frame.save(path)
                    frame_paths.append(path)
                    if use_blit:
                        for artist in artists:
                            artist.set_animated(True)
                        background[0] = None
            now = time.perf_counter()
            if not args.no_live and now-last_events > 0.1:
                fig.canvas.flush_events()
                last_events = now
            if n == 0 or n == steps-1 or now-last_progress > 2:
                print(f' tU/D={n*U/D:6.2f}; y/D={y/D:+.4f}; elapsed {now-started:.1f}s',flush=True)
                last_progress = now
            if not args.fixed:
                y,v = oscillator_step(y,v,force[1],mass,damping,stiffness)
    # Save fields and cylinder geometry at the last evaluated time, n*dt.
    last_y = records[-1][1]*D
    np.savez_compressed(out/'final_field.npz',ux=ux,uy=uy,rho=rho,
                        mask=mask,D=D,U=U,xc=xc,yc=yc0+last_y,time=records[-1][0])
    for artist in artists:
        artist.set_animated(False)
    fig.savefig(out/'wake_and_motion.png',dpi=160)
    data = np.asarray(records)
    fig2, axes = plt.subplots(3,1,figsize=(9,7),sharex=True,constrained_layout=True)
    for axis,column,label in zip(axes,[1,3,4],['y/D','Cd','Cl']):
        axis.plot(data[:,0],data[:,column],lw=1)
        axis.set_ylabel(label)
        axis.grid(alpha=0.3)
    axes[-1].set_xlabel('t U/D')
    fig2.savefig(out/'response.png',dpi=160)
    # Rough dominant lift frequency from final half, excluding startup ramp.
    tail = data[data[:,0]>max(10,0.5*args.duration)]
    if len(tail)>64:
        # Final irregular sample, if any, is excluded from FFT.
        uniform = tail[:-1]
        signal = uniform[:,4]-uniform[:,4].mean()
        spectrum = np.abs(np.fft.rfft(signal*np.hanning(len(signal))))
        freqs = np.fft.rfftfreq(len(signal),sample_stride*U/D)
        dominant = freqs[1+np.argmax(spectrum[1:])]
        print(f'Dominant lift fD/U ~ {dominant:.3f}; interpret only after checking periodic steady response.')
        print(f'Last-half mean Cd={tail[:,3].mean():.3f}, max |y/D|={np.abs(tail[:,1]).max():.3f}')
        if frame_paths:
            import imageio.v2 as imageio
    
            video_path = out / 'fsi_wake_60fps.mp4'
    
            print('Saving 60 fps video...', flush=True)
    
            with imageio.get_writer(
                str(video_path),
                format='FFMPEG',
                mode='I',
                fps=args.fps,
                codec='libx264',
                quality=8,
                macro_block_size=2
            ) as writer:
                for path in frame_paths:
                    writer.append_data(imageio.imread(path))
    
            print(f'Video saved: {video_path.resolve()}')
    frame_store.cleanup()
    print(f'Saved results to {out.resolve()}')
    if not args.no_live:
        plt.ioff()
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    main()
