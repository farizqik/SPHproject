"""Shared data reader / native-polar streamfunction reconstruction for EWCSPH.
Requires numpy, pandas, scipy, matplotlib. Python 3.10+.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from functools import lru_cache
import csv, re
import numpy as np
import pandas as pd
from scipy.linalg import solve_banded

COLUMNS = ['ID','x','y','empty1','rho','drhodt','pressure','empty2',
           'u','v','velocity','empty3','dudt','dvdt','vorticity','Type']
TIME_RE = re.compile(r'EWCSPH_hdp_([-+\d.eE]+)_t_([-+\d.eE]+)\.csv$', re.I)
FIELD_NAMES = {
    'velocity': r'$|\mathbf{u}|/U_\infty$',
    'u': r'$u/U_\infty$', 'v': r'$v/U_\infty$',
    'rho': r'$\rho/\rho_0$', 'pressure': r'$C_p$',
    'vorticity': r'$\omega D/U_\infty$',
    'drhodt': r'$\partial\rho/\partial t$ (kg m$^{-3}$ s$^{-1}$)',
    'dudt': r'$\partial u/\partial t$ (m s$^{-2}$)',
    'dvdt': r'$\partial v/\partial t$ (m s$^{-2}$)'
}

def load_metadata(folder):
    p=Path(folder)/'SimulationConstants.csv'
    if not p.is_file(): raise FileNotFoundError(f'Missing {p}')
    m={}
    with p.open(encoding='utf-8-sig',newline='') as f:
        for line in csv.reader(f):
            if len(line)>=2 and line[0].strip(): m[line[0].strip()]=line[1].strip()
    return m

def numeric(meta,key,default=None):
    try: return float(meta[key])
    except (KeyError, ValueError,TypeError):
        if default is not None:return float(default)
        raise ValueError(f'Missing/invalid {key} in SimulationConstants.csv')

def read_particles(path):
    # Six header/metadata lines precede the particle data in C++ CSV output.
    frame=pd.read_csv(path,skiprows=6,header=None,names=COLUMNS,
                      usecols=list(range(16)),low_memory=False)
    for k in COLUMNS:
        if k not in ('empty1','empty2','empty3','Type'):
            frame[k]=pd.to_numeric(frame[k],errors='coerce')
    frame['Type']=frame['Type'].astype(str).str.strip().str.lower()
    return frame

@dataclass
class FrameInfo:
    time: float
    file: Path

class Case:
    def __init__(self,folder,max_frames=1201,hdp=None):
        self.folder=Path(folder).expanduser().resolve()
        self.meta=load_metadata(self.folder)
        self.D=numeric(self.meta,'Dcylinder')
        self.U=numeric(self.meta,'inletVelocity')
        self.rho0=numeric(self.meta,'rho0')
        self.nt=int(numeric(self.meta,'NTheta'))
        self.Re=numeric(self.meta,'Reynolds',0)
        self.phi=np.deg2rad(numeric(self.meta,'flowAngleDeg',0))
        # Do not infer free-slip from viscosity or Reynolds number.
        bc=str(self.meta.get('boundaryCondition','')).strip().lower().replace('_','-').replace(' ','-')
        scenario=str(self.meta.get('testScenario','')).strip().lower()
        self.free_slip = (bc in {'free-slip','freeslip','slip'} or
                          (not bc and scenario in {'potential','potential-flow','inviscid'}))
        self.has_potential_reference = self.free_slip
        if self.D<=0 or self.U<=0 or self.nt<4: raise ValueError('Invalid D, U, or NTheta')
        items=[]
        for path in self.folder.glob('EWCSPH_hdp_*_t_*.csv'):
            mat=TIME_RE.fullmatch(path.name)
            if mat:
                try:items.append((float(mat[1]),float(mat[2]),path))
                except ValueError:pass
        if not items:raise ValueError(f'No EWCSPH timestep CSV files in {self.folder}')
        self.h_options=sorted({x[0] for x in items})
        self.hdp=hdp if hdp is not None else self.h_options[0]
        items=sorted((x for x in items if abs(x[0]-self.hdp)<1e-10*max(1,abs(self.hdp))),key=lambda x:x[1])
        seen=set();self.frames=[]
        for _,t,p in items:
            if t not in seen:self.frames.append(FrameInfo(t,p));seen.add(t)
        self.frames=self.frames[:max_frames]
        self.times=np.array([f.time for f in self.frames])
        self._first=self.read(0)
        types=self._first['Type'].to_numpy()
        self.fluid=np.flatnonzero(types=='fluid')
        self.wall=np.flatnonzero(types=='boundary')
        if not len(self.fluid) or len(self.fluid)%self.nt:raise ValueError('Fluid rows do not form complete polar rings')
        self.nr=len(self.fluid)//self.nt
        if self.nr<2:raise ValueError('At least two fluid rings required')
        x=self._first['x'].to_numpy();y=self._first['y'].to_numpy()
        self.center=np.array([np.mean(x[self.wall]),np.mean(y[self.wall])]) if len(self.wall) else np.array([np.mean(x[self.fluid[:self.nt]]),np.mean(y[self.fluid[:self.nt]])])
        self.x=(x-self.center[0])/self.D;self.y=(y-self.center[1])/self.D
        theta=np.mod(np.arctan2(self.y[self.fluid[:self.nt]],self.x[self.fluid[:self.nt]]),2*np.pi)
        self.order=np.argsort(theta)
        self.theta=theta[self.order]
        rr=np.hypot(self.x[self.fluid],self.y[self.fluid]).reshape((self.nr,self.nt)).T
        self.radii=rr.mean(axis=0)
        if np.any(np.diff(self.radii)<=0):raise ValueError('Fluid ring radii must increase')
        self.wall_r=float(np.mean(np.hypot(self.x[self.wall],self.y[self.wall]))) if len(self.wall) else self.radii[0]
        self.R=np.r_[self.wall_r,self.radii] if self.wall_r<self.radii[0] else self.radii.copy()
        self.has_wall=self.wall_r<self.radii[0]
        self.X=(self.x[self.fluid].reshape((self.nr,self.nt)).T)[self.order]
        self.Y=(self.y[self.fluid].reshape((self.nr,self.nt)).T)[self.order]
        self.Xplot=np.vstack([self.X,self.X[0]])
        self.Yplot=np.vstack([self.Y,self.Y[0]])
        self.Xwall=self.x[self.wall];self.Ywall=self.y[self.wall]
        if self.has_wall:self._build_poisson()

    @lru_cache(maxsize=24)
    def read(self,i):
        return read_particles(self.frames[int(i)].file)

    @lru_cache(maxsize=1)
    def potential_reference(self):
        """Analytical inviscid, irrotational flow past an impermeable cylinder.

        Values are returned in the *same nondimensional convention* as values().
        The velocity reference is steady and is independent of frame number.
        Only available when C++ metadata explicitly identifies a free-slip test.
        """
        if not self.has_potential_reference:
            raise ValueError('Analytical cylinder reference is reserved for free-slip simulations')
        theta=np.arctan2(self.y,self.x)
        rr=np.hypot(self.x,self.y)
        # Cylinder radius in units of D (wall_r normally 0.5).
        a=self.wall_r
        r_safe=np.maximum(rr,a)
        q=(a/r_safe)**2
        beta=theta-self.phi
        ur=(1-q)*np.cos(beta)
        ut=-(1+q)*np.sin(beta)
        ux=ur*np.cos(theta)-ut*np.sin(theta)
        uy=ur*np.sin(theta)+ut*np.cos(theta)
        speed2=ux*ux+uy*uy
        cp=1-speed2
        pressure=0.5*self.rho0*self.U**2*cp
        # Reconstruct density only when barotropic EOS constants are recorded.
        B=numeric(self.meta,'B',np.nan)
        gamma=numeric(self.meta,'gammaEOS',np.nan)
        if np.isfinite(B) and B>0 and np.isfinite(gamma) and gamma>0:
            base=np.maximum(1.0+pressure/B,1.e-12)
            rho_ratio=base**(1.0/gamma)
        else:
            rho_ratio=np.ones_like(cp)
        return {'u':ux,'v':uy,'velocity':np.sqrt(speed2),
                'pressure':cp,'vorticity':np.zeros_like(cp),
                'rho':rho_ratio,'drhodt':np.zeros_like(cp),
                'dudt':np.zeros_like(cp),'dvdt':np.zeros_like(cp)}

    def reference_values(self,field):
        return self.potential_reference()[field]

    def nearest(self,t): return int(np.argmin(np.abs(self.times-t)))

    def values(self,i,field):
        d=self.read(i)
        if field=='pressure':v=d['pressure'].to_numpy()/(0.5*self.rho0*self.U**2)
        elif field=='rho':v=d['rho'].to_numpy()/self.rho0
        elif field=='vorticity':v=d['vorticity'].to_numpy()*self.D/self.U
        elif field in ('u','v','velocity'):v=d[field].to_numpy()/self.U
        else:v=d[field].to_numpy()
        return v

    def polar(self,array):
        return array[self.fluid].reshape((self.nr,self.nt)).T[self.order]

    def _build_poisson(self):
        # Fixed geometry, periodic theta; solve each azimuthal Fourier mode by
        # a radial tridiagonal system. Same discretisation as MATLAB viewer.
        r=self.R;nt=self.nt;drtheta=2*np.pi/nt
        if np.max(np.abs(np.diff(np.r_[self.theta,self.theta[0]+2*np.pi])-drtheta))>0.05*drtheta:
            raise ValueError('Streamfunction solver requires nearly uniform angular spacing')
        ri=r[1:-1];rm=r[:-2];rp=r[2:]
        drcell=(rp-rm)/2
        self.aa=(ri+rm)/2/(ri*drcell*(ri-rm))
        self.cc=(ri+rp)/2/(ri*drcell*(rp-ri))
        self.freq=np.fft.fftfreq(nt,d=1/nt).astype(int)
        eig=4*np.sin(np.pi*np.arange(nt)/nt)**2/drtheta**2
        self.diag=-(self.aa[None,:]+self.cc[None,:])-eig[:,None]/(ri[None,:]**2)
        self.lower=self.aa.copy();self.upper=self.cc.copy()

    def _curl(self,u,v):
        theta=self.theta[:,None];rr=self.radii
        ur=u*np.cos(theta)+v*np.sin(theta)
        ut=-u*np.sin(theta)+v*np.cos(theta)
        a=ut*rr[None,:]
        da=np.gradient(a,rr,axis=1,edge_order=2 if len(rr)>2 else 1)
        durdt=(np.roll(ur,-1,axis=0)-np.roll(ur,1,axis=0))/(2*2*np.pi/self.nt)
        return (da-durdt)/rr[None,:]

    @lru_cache(maxsize=12)
    def streamfunction(self,i,vorticity_source='velocity'):
        if not self.has_wall: return None
        d=self.read(i)
        u=self.polar(d['u'].to_numpy())/self.U
        v=self.polar(d['v'].to_numpy())/self.U
        if vorticity_source=='solver':omega=self.polar(d['vorticity'].to_numpy())*self.D/self.U
        else:omega=self._curl(u,v)
        if not np.all(np.isfinite(omega)):return None
        ur=u[:,-1]*np.cos(self.theta)+v[:,-1]*np.sin(self.theta)
        ur=ur-ur.mean()
        flux_hat=np.fft.fft(self.R[-1]*ur)
        psi_hat=np.zeros(self.nt,dtype=complex)
        usable=self.freq!=0
        if self.nt%2==0:usable[self.nt//2]=False
        psi_hat[usable]=flux_hat[usable]/(1j*self.freq[usable])
        psi_outer=np.fft.ifft(psi_hat).real
        rhs=np.fft.fft(-omega[:,:-1],axis=0)
        rhs[:,-1]-=self.cc[-1]*np.fft.fft(psi_outer)
        n=rhs.shape[1]
        solution=np.zeros_like(rhs)
        for m in range(self.nt):
            ab=np.zeros((3,n),dtype=float)
            ab[1]=self.diag[m]
            if n>1:
                ab[0,1:]=self.upper[:-1]
                ab[2,:-1]=self.lower[1:]
            solution[m]=solve_banded((1,1),ab,rhs[m],check_finite=False)
        return np.fft.ifft(np.column_stack([np.zeros(self.nt),solution,np.fft.fft(psi_outer)]),axis=0).real

    def contour_geometry(self,psi):
        a=np.r_[self.theta,self.theta[0]+2*np.pi]
        ps=np.vstack([psi,psi[0]])
        rs=np.broadcast_to(self.R[None,:],ps.shape)
        return rs*np.cos(a[:,None]),rs*np.sin(a[:,None]),ps

    def initial_stats(self,field='vorticity'):
        samples=np.unique(np.linspace(0,len(self.frames)-1,min(5,len(self.frames)),dtype=int))
        chunks=[]
        for i in samples:
            z=self.values(int(i),field)[self.fluid];chunks.append(z[np.isfinite(z)][::max(1,len(z)//30000)])
        dat=np.concatenate(chunks);lo,hi=np.percentile(dat,[5,95])
        if field=='vorticity':hi=max(abs(lo),abs(hi));lo=-hi
        elif field=='velocity':lo=0
        if hi<=lo:hi=lo+1
        return (float(lo),float(hi))

    @lru_cache(maxsize=1)
    def summary_history(self):
        vals=np.full((len(self.frames),7),np.nan)
        for k,f in enumerate(self.frames):
            try:
                with f.file.open(encoding='utf-8-sig') as fh:
                    next(fh);next(fh);line=next(fh)
                nums=[float(x) for x in line.strip().split(',')[:7]]
                vals[k,:len(nums)]=nums
            except (OSError,ValueError,StopIteration):pass
        return vals

    @lru_cache(maxsize=1)
    def global_forces(self):
        files=list(self.folder.glob('PotentialFlowGlobalComparison_hdp_*.csv'))
        if not files:return None
        hit=[]
        for p in files:
            try:
                h=float(p.stem.split('hdp_')[-1])
                if abs(h-self.hdp)<1e-10*max(1,self.hdp):hit.append(p)
            except ValueError:pass
        try:return pd.read_csv(hit[0]) if hit else None
        except Exception:return None

    @lru_cache(maxsize=1)
    def local_comparison(self):
        paths=list(self.folder.glob('PotentialFlowLocalComparison_hdp_*.csv'))
        for p in paths:
            try:
                h=float(p.stem.split('hdp_')[-1])
                if abs(h-self.hdp)<1e-10*max(1,self.hdp):return pd.read_csv(p)
            except ValueError:pass
        return None


def radial_profile(case,index,field,angle_deg):
    """Interpolate along theta between real SPH rings; avoids nearest-angle jumps."""
    theta=np.mod(np.deg2rad(float(angle_deg)),2*np.pi)
    a=case.theta
    k=int(np.searchsorted(a,theta))
    left=(k-1)%case.nt;right=k%case.nt
    al=float(a[left]);ar=float(a[right])
    if right==0:ar+=2*np.pi
    if theta<al:theta+=2*np.pi
    w=(theta-al)/(ar-al)
    zz=case.polar(case.values(index,field))
    return (1-w)*zz[left,:]+w*zz[right,:]

def sample_line(case,index,field,A,B,n=300):
    """Fast polar interpolation on the original rings without rebuilding Delaunay."""
    from scipy.interpolate import RegularGridInterpolator
    A=np.asarray(A,dtype=float);B=np.asarray(B,dtype=float)
    distance=float(np.linalg.norm(B-A));ss=np.linspace(0,1,n)
    XY=A[None,:]+ss[:,None]*(B-A)[None,:]
    r=np.hypot(XY[:,0],XY[:,1]);ang=np.mod(np.arctan2(XY[:,1],XY[:,0]),2*np.pi)
    z=case.polar(case.values(index,field))
    theta=np.r_[case.theta[-1]-2*np.pi,case.theta,case.theta[0]+2*np.pi]
    z=np.vstack([z[-1],z,z[0]])
    interp=RegularGridInterpolator((theta,case.radii),z,bounds_error=False,fill_value=np.nan)
    result=interp(np.column_stack([ang,r]));result[r<=case.wall_r]=np.nan
    return ss*distance,result
