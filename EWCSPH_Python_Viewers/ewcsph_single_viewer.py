"""Interactive single-case EWCSPH Python viewer — run: python ewcsph_single_viewer.py"""
from __future__ import annotations
import argparse
import sys
import subprocess
from pathlib import Path
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import imageio.v2 as imageio
from ewcsph_core import Case,FIELD_NAMES,radial_profile,sample_line
from ewcsph_gui import add_figure,draw_field,folder_information,journal_axes,pretty_column,launch_viewer

class SingleViewer:
    def __init__(self,root,folder):
        self.root=root;self.case=Case(folder);self.i=0;self.playing=False;self.after_id=None
        self.root.title('EWCSPH — Single Simulation Viewer');self.root.geometry('1350x850')
        self.field=tk.StringVar(value='velocity');self.psi=tk.BooleanVar(value=True)
        self.spacing=tk.StringVar(value='0.06');self.source=tk.StringVar(value='velocity')
        self.near=tk.BooleanVar(value=True);self.fps=tk.StringVar(value='12')
        self.show_wall=tk.BooleanVar(value=True)
        self.profile_angle=tk.StringVar(value='90');self.profile_field=tk.StringVar(value='velocity')
        self.line_points=[];self.line_active=False
        self.limits={};self._drawing=False;self._suspend_slider=False;self._last_aux={};self._play_count=0
        self.force_group=tk.StringVar(value="Drag coefficient")
        self.local_quantity=tk.StringVar(value="Cp")
        self.local_mode=tk.StringVar(value="Angle at current time")
        self.local_angle=tk.StringVar(value="90")
        self.local_layer=tk.StringVar(value="Boundary (0)")
        self.show_contours_play=tk.BooleanVar(value=False);self.play_stride=tk.StringVar(value='2')
        self.window_width=tk.StringVar(value='1350');self.window_height=tk.StringVar(value='850')
        self.tangential_mode=tk.StringVar(value='Multiple');self.tangential_angles={a:tk.BooleanVar(value=(a in (0,45,90,135,180))) for a in range(0,360,15)}
        self.exporting=False
        top=ttk.Frame(root,padding=6);top.pack(fill='x')
        ttk.Label(top,text='Field').pack(side='left')
        ttk.Combobox(top,textvariable=self.field,values=list(FIELD_NAMES),width=13,state='readonly').pack(side='left',padx=3)
        ttk.Checkbutton(top,text='Streamfunction contours',variable=self.psi,command=self.draw).pack(side='left',padx=8)
        ttk.Label(top,text='Δψ/(U∞D)').pack(side='left')
        x=ttk.Entry(top,textvariable=self.spacing,width=7);x.pack(side='left',padx=3);x.bind('<Return>',lambda e:self.draw())
        ttk.Checkbutton(top,text='Near-zero levels',variable=self.near,command=self.draw).pack(side='left',padx=5)
        ttk.Label(top,text='Vorticity:').pack(side='left',padx=(10,2))
        ttk.Combobox(top,textvariable=self.source,values=['velocity','solver'],state='readonly',width=9).pack(side='left')
        ttk.Checkbutton(top,text='Cylinder',variable=self.show_wall,command=self.draw).pack(side='left',padx=5)
        ttk.Checkbutton(top,text='Contours during playback',variable=self.show_contours_play).pack(side='left',padx=5)
        ttk.Button(top,text='Open New Viewer',command=self.open_new_viewer).pack(side='right',padx=4)
        ttk.Button(top,text='Change Folder',command=self.change_folder).pack(side='right',padx=4)
        ttk.Button(top,text='Folder info',command=lambda:folder_information(root,[self.case.folder])).pack(side='right')
        ttk.Button(top,text='Window size',command=self.window_settings).pack(side='right',padx=5)
        ttk.Button(top,text='Save MP4',command=self.save_mp4).pack(side='right',padx=5)
        nav=ttk.Frame(root,padding=5);nav.pack(fill='x')
        ttk.Button(nav,text='|<',command=lambda:self.set_frame(0)).pack(side='left')
        ttk.Button(nav,text='<',command=lambda:self.set_frame(self.i-1)).pack(side='left')
        self.play_button=ttk.Button(nav,text='Play',command=self.toggle_play);self.play_button.pack(side='left',padx=4)
        ttk.Button(nav,text='>',command=lambda:self.set_frame(self.i+1)).pack(side='left')
        ttk.Button(nav,text='>|',command=lambda:self.set_frame(len(self.case.frames)-1)).pack(side='left')
        self.slider=tk.Scale(nav,from_=0,to=len(self.case.frames)-1,orient='horizontal',showvalue=False,command=lambda v:self.set_frame(int(float(v)),from_slider=True))
        self.slider.pack(side='left',fill='x',expand=True,padx=15)
        ttk.Label(nav,text='FPS').pack(side='left');ttk.Entry(nav,textvariable=self.fps,width=5).pack(side='left',padx=3)
        ttk.Label(nav,text='Frame step').pack(side='left',padx=(8,2))
        ttk.Combobox(nav,textvariable=self.play_stride,values=['1','2','4','8'],width=3,state='readonly').pack(side='left')
        self.status=ttk.Label(nav,text='');self.status.pack(side='left',padx=10)
        self.tabs=ttk.Notebook(root);self.tabs.pack(fill='both',expand=True)
        self.tabmap={}
        tab_names=['Flow Field','L2 vs Time','Hydrodynamic Forces','Local Comparison',
                   'Tangential vs Radius','Radial Profiles','Particle Distribution','Line Profile']
        if self.case.has_potential_reference:
            tab_names.insert(1,'Potential-Flow Validation')
        for name in tab_names:
            f=ttk.Frame(self.tabs);self.tabs.add(f,text=name);self.tabmap[name]=f
        self.fig=Figure(figsize=(10,6),dpi=100)
        self.ax=self.fig.add_subplot(111)
        self.ax_ref=None
        self._reference_field=None
        if self.case.has_potential_reference:
            # The second panel is only created for explicitly free-slip data.
            self.ax.set_position([.07,.14,.36,.76])
            self.ax_ref=self.fig.add_axes([.55,.14,.36,.76],sharex=self.ax,sharey=self.ax)
        else:
            self.fig.subplots_adjust(left=.11,right=.91,bottom=.11,top=.93)
        _,self.canvas=add_figure(self.tabmap['Flow Field'],self.fig)
        self.cid=self.canvas.mpl_connect('button_press_event',self.click_flow)
        self.aux={}
        for name in list(self.tabmap)[1:]:
            fig=Figure(figsize=(10,5),dpi=90)
            ax=fig.add_subplot(111)
            _,c=add_figure(self.tabmap[name],fig)
            self.aux[name]=(fig,ax,c)
        self.tabs.bind('<<NotebookTabChanged>>',lambda e:self.draw())
        self.field.trace_add('write',lambda *_:self.draw())
        self.profile_field.trace_add('write',lambda *_:self.draw())
        self.source.trace_add('write',lambda *_:self.draw())
        self._create_analysis_controls()
        self._create_matlab_controls()
        if self.case.has_potential_reference:
            self.root.title('EWCSPH — Single Simulation Viewer | Free-slip / potential-flow reference')
        self.root.protocol('WM_DELETE_WINDOW',self.close)
        self.draw()

    def _create_analysis_controls(self):
        f=self.tabmap['Radial Profiles'];bar=ttk.Frame(f);bar.pack(fill='x',before=self.aux['Radial Profiles'][2].get_tk_widget())
        ttk.Label(bar,text='θ (°)').pack(side='left',padx=5)
        ttk.Entry(bar,textvariable=self.profile_angle,width=10).pack(side='left')
        ttk.Label(bar,text='Variable').pack(side='left',padx=5)
        ttk.Combobox(bar,textvariable=self.profile_field,values=list(FIELD_NAMES),state='readonly',width=16).pack(side='left')
        ttk.Button(bar,text='Update',command=self.draw).pack(side='left',padx=5)
        f=self.tabmap['Line Profile'];bar=ttk.Frame(f);bar.pack(fill='x',before=self.aux['Line Profile'][2].get_tk_widget())
        self.line_label=ttk.Label(bar,text='Click “Pick A and B”, then click two points on Flow Field')
        self.line_label.pack(side='left',padx=6)
        ttk.Button(bar,text='Pick A and B',command=self.start_line).pack(side='left',padx=5)
        ttk.Button(bar,text='Clear',command=self.clear_line).pack(side='left')

    def _create_matlab_controls(self):
        # MATLAB-like field selectors for force and local comparison tabs.
        f=self.tabmap['Hydrodynamic Forces']
        bar=ttk.Frame(f);bar.pack(fill='x',before=self.aux['Hydrodynamic Forces'][2].get_tk_widget())
        ttk.Label(bar,text='Global quantity:').pack(side='left',padx=6)
        ttk.Combobox(bar,textvariable=self.force_group,state='readonly',width=29,
            values=['Force Fx','Force Fy','Force magnitude','Drag coefficient',
                    'Lift coefficient','Resultant force coefficient']).pack(side='left',padx=4)
        ttk.Button(bar,text='Update',command=self.draw).pack(side='left')
        f=self.tabmap['Local Comparison']
        bar=ttk.Frame(f);bar.pack(fill='x',before=self.aux['Local Comparison'][2].get_tk_widget())
        ttk.Label(bar,text='Plot against:').pack(side='left',padx=3)
        ttk.Combobox(bar,textvariable=self.local_mode,width=23,state='readonly',
                     values=['Angle at current time','Time at selected angle']).pack(side='left',padx=3)
        ttk.Label(bar,text='Parameter:').pack(side='left',padx=3)
        vals=['Cp','Pressure','u velocity','v velocity','Radial velocity','Tangential velocity',
              'Velocity magnitude','Wall shear stress','Local drag force','Local lift force']
        ttk.Combobox(bar,textvariable=self.local_quantity,width=19,state='readonly',values=vals).pack(side='left',padx=3)
        ttk.Label(bar,text='Layer:').pack(side='left',padx=3)
        ttk.Combobox(bar,textvariable=self.local_layer,width=17,state='readonly',
                     values=['Boundary (0)','Fluid ring 1 (1)','Fluid ring 2 (2)']).pack(side='left',padx=3)
        ttk.Label(bar,text='θ (°):').pack(side='left',padx=3)
        ttk.Entry(bar,textvariable=self.local_angle,width=6).pack(side='left')
        ttk.Button(bar,text='Update',command=self.draw).pack(side='left',padx=3)
        f=self.tabmap['Tangential vs Radius']
        bar=ttk.Frame(f);bar.pack(fill='x',before=self.aux['Tangential vs Radius'][2].get_tk_widget())
        ttk.Label(bar,text='θ (°):').pack(side='left',padx=6)
        ttk.Combobox(bar,textvariable=self.tangential_mode,values=['Single','Multiple'],state='readonly',width=10).pack(side='left')
        ttk.Button(bar,text='Select angles (checklist)',command=self.select_tangential_angles).pack(side='left',padx=7)
        ttk.Button(bar,text='Update',command=self.draw).pack(side='left')

    def _draw_potential_reference(self,field,limits):
        """Update a persistent analytical field, using the SPH field's colour scale."""
        import matplotlib.tri as mtri
        ax=self.ax_ref
        if ax is None:return
        case=self.case
        vals=case.polar(case.reference_values(field)).ravel()
        tri=mtri.Triangulation(case.X.ravel(),case.Y.ravel(),case.triangles)
        if not hasattr(ax,'_reference_mesh'):
            from matplotlib.patches import Circle
            mesh=ax.tripcolor(tri,vals,shading='gouraud',cmap='RdBu_r' if field=='vorticity' else 'turbo',rasterized=True)
            ax._reference_mesh=mesh
            ax.add_patch(Circle((0,0),case.wall_r,facecolor='white',edgecolor='black',linewidth=.9,zorder=4))
            ax.set_aspect('equal',adjustable='box')
            ax.set(xlabel=r'$x/D$',ylabel=r'$y/D$')
            journal_axes(ax)
        else:
            mesh=ax._reference_mesh
            mesh.set_array(vals)
            mesh.set_cmap('RdBu_r' if field=='vorticity' else 'turbo')
        mesh.set_clim(*limits)
        ax.set_title('Analytical potential flow',pad=8)
        ax.set_xlim(self.ax.get_xlim());ax.set_ylim(self.ax.get_ylim())

    def _potential_validation(self,fig):
        """Current-frame SPH vs analytical profiles; no reference for no-slip."""
        case=self.case
        fig.clear()
        axs=fig.subplots(2,2)
        d=case.read(self.i)
        angles=np.rad2deg(case.theta)
        a=case.wall_r
        u=case.polar(d['u'].to_numpy())/case.U
        v=case.polar(d['v'].to_numpy())/case.U
        ur=u*np.cos(case.theta[:,None])+v*np.sin(case.theta[:,None])
        ut=-u*np.sin(case.theta[:,None])+v*np.cos(case.theta[:,None])
        refu=case.polar(case.reference_values('u'))
        refv=case.polar(case.reference_values('v'))
        refut=-refu*np.sin(case.theta[:,None])+refv*np.cos(case.theta[:,None])
        # Compare first FLUID ring Cp with potential Cp at identical locations.
        cp=case.polar(case.values(self.i,'pressure'))[:,0]
        cp_ref=case.polar(case.reference_values('pressure'))[:,0]
        axs[0,0].plot(angles,cp,'.',markersize=3,label='SPH (first fluid ring)')
        axs[0,0].plot(angles,cp_ref,'-',linewidth=1.25,label='Potential flow (same radius)')
        axs[0,0].set(xlabel=r'$\theta$ (°)',ylabel=r'$C_p$',title='Pressure coefficient')
        axs[0,0].legend(fontsize=7)
        # Radial tangential velocity at angular positions.
        for angle,color in [(90,'tab:blue'),(270,'tab:orange')]:
            k=int(np.argmin(np.abs((angles-angle+180)%360-180)))
            axs[0,1].plot(case.radii-a,ut[k],color=color,lw=1.2,label=fr'SPH ${angle}^\circ$')
            axs[0,1].plot(case.radii-a,refut[k],color=color,ls='--',lw=1.1,label=fr'Exact ${angle}^\circ$')
        axs[0,1].set(xlabel=r'$(r-r_0)/D$',ylabel=r'$u_\theta/U_\infty$',title='Tangential velocity')
        axs[0,1].legend(fontsize=7,ncol=2)
        dv=np.sqrt((u-refu)**2+(v-refv)**2)
        axs[1,0].plot(case.radii-a,np.sqrt(np.nanmean(dv**2,axis=0)),color='tab:red',lw=1.3)
        axs[1,0].set(xlabel=r'$(r-r_0)/D$',ylabel=r'$\mathrm{RMS}(|\Delta\mathbf{u}|/U_\infty)$',title='Velocity error by radius')
        errs=case.summary_history()
        for column,label in [(0,r'$L_2(|\mathbf{u}|)$'),(1,r'$L_2(u)$'),(2,r'$L_2(v)$')]:
            y=errs[:,column]
            mask=np.isfinite(y)&(y>0)
            if np.any(mask):axs[1,1].semilogy(case.times[mask],y[mask],label=label,lw=1.1)
        axs[1,1].axvline(case.times[self.i],color='black',ls=':',lw=.8)
        axs[1,1].set(xlabel=r'$t$ (s)',ylabel=r'$L_2$',title='C++ potential-flow diagnostics')
        if axs[1,1].lines:axs[1,1].legend(fontsize=7)
        for ax in axs.flat:journal_axes(ax,grid=True)
        fig.tight_layout(pad=2)

    def _local_chart(self,ax):
        import pandas as pd
        case=self.case;df=case.local_comparison()
        if df is None or df.empty:
            ax.text(.5,.5,'PotentialFlowLocalComparison CSV not available',ha='center',transform=ax.transAxes)
            return
        cols=list(df.columns);low={c:c.lower().replace(' ','').replace('_','') for c in cols}
        def find_col(*words):
            for term in words:
                for col,norm in low.items():
                    if norm==term.lower().replace('_','').replace(' ',''):return col
            for term in words:
                for col,norm in low.items():
                    if term.lower().replace('_','').replace(' ','') in norm:return col
            return None
        q=self.local_quantity.get()
        names={
            'Cp':('CpSPH','Cp','pressureCoefficient'),
            'Pressure':('pressureSPH','pSPH','pressure'),
            'u velocity':('uSPH','uxSPH','u'),
            'v velocity':('vSPH','vySPH','v'),
            'Radial velocity':('urSPH','uRadialSPH'),
            'Tangential velocity':('utSPH','uTangentialSPH'),
            'Velocity magnitude':('velocitySPH','speedSPH','vMagSPH'),
            'Wall shear stress':('tauWallSPH','wallShearStress'),
            'Local drag force':('dFxSPH','FxSPH'),
            'Local lift force':('dFySPH','FySPH')}
        vcol=find_col(*(names.get(q,(q,))))
        if vcol is None:
            numeric=list(df.select_dtypes('number').columns)
            vcol=numeric[-1] if numeric else None
        if vcol is None:
            ax.text(.5,.5,'No numeric values in local-comparison CSV',ha='center',transform=ax.transAxes);return
        timecol=find_col('t','time');angcol=find_col('angleDeg','thetaDeg','angle','theta');layercol=find_col('layer','ring','particleLayer')
        view=df
        if layercol:
            layer={'Boundary (0)':0,'Fluid ring 1 (1)':1,'Fluid ring 2 (2)':2}.get(self.local_layer.get(),0)
            cand=view[pd.to_numeric(view[layercol],errors='coerce')==layer]
            if not cand.empty:view=cand
        mode=self.local_mode.get()
        if mode=='Angle at current time':
            if timecol:
                tarr=pd.to_numeric(view[timecol],errors='coerce')
                if tarr.notna().any():view=view[np.isclose(tarr,tarr.iloc[(tarr-self.case.times[self.i]).abs().argmin()])]
            if angcol:
                view=view.sort_values(angcol)
                ax.plot(view[angcol],view[vcol],'.-',markersize=3);ax.set_xlabel(pretty_column(angcol))
            else:
                ax.plot(np.arange(len(view)),view[vcol]);ax.set_xlabel('Local record')
        else:
            if angcol:
                aa=pd.to_numeric(view[angcol],errors='coerce')
                target=float(self.local_angle.get())
                dif=np.abs((aa-target+180)%360-180)
                good=np.isfinite(dif)
                if good.any():view=view[np.isclose(dif,np.min(dif[good]),atol=1e-4)]
            if timecol:
                view=view.sort_values(timecol);ax.plot(view[timecol],view[vcol],'.-');ax.set_xlabel(r'$t$ (s)')
            else:ax.plot(np.arange(len(view)),view[vcol]);ax.set_xlabel('Local record')
        ax.set_ylabel(pretty_column(vcol) if pretty_column(vcol)!=vcol else vcol + self._local_unit(vcol));ax.set_title(f'{q} — {mode} | {self.local_layer.get()}')

    def _force_chart(self,ax):
        df=self.case.global_forces()
        if df is None or df.empty:
            ax.text(.5,.5,'PotentialFlowGlobalComparison CSV not available',ha='center',transform=ax.transAxes);return
        group=self.force_group.get()
        pairs={
          'Force Fx':['FxPressureSPH','FxViscousSPH','FxSPH','FxExact'],
          'Force Fy':['FyPressureSPH','FyViscousSPH','FySPH','FyExact'],
          'Force magnitude':['FPressureSPH','FViscousSPH','FSPH','FExact'],
          'Drag coefficient':['CDPressureSPH','CDViscousSPH','CDSPH','CDExact'],
          'Lift coefficient':['CLPressureSPH','CLViscousSPH','CLSPH','CLExact'],
          'Resultant force coefficient':['CFPressureSPH','CFViscousSPH','CFSPH','CFExact']}
        tt=df['t'] if 't' in df else np.arange(len(df))
        for key in pairs[group]:
            if key.endswith('Exact') and not self.case.has_potential_reference:
                continue
            if key in df:ax.plot(tt,df[key],label=key)
        ax.axvline(self.case.times[self.i],c='black',ls='--',lw=.8)
        if ax.lines:ax.legend(fontsize=8)
        ax.set(xlabel=r'$t$ (s)',ylabel=group+(' (N/m)' if group.startswith('Force') else ''),title=group)
    @staticmethod
    def _local_unit(column):
        c=column.lower()
        if 'coeff' in c or c.startswith(('cp','cd','cl','cf','cu','cv')): return ''
        if 'pressure' in c or c.startswith('psph'): return ' [Pa]'
        if 'stress' in c or 'tau' in c: return ' [Pa]'
        if 'force' in c or c.startswith(('fx','fy','dfx','dfy')): return ' [N/m]'
        if 'vorticity' in c: return ' [1/s]'
        if 'velocity' in c or 'speed' in c or c.startswith(('usph','vsph','ut','ur')): return ' [m/s]'
        return ''

    def _stop_for_reload(self):
        self.playing=False
        if self.after_id is not None:
            self.root.after_cancel(self.after_id);self.after_id=None

    def change_folder(self):
        """Replace the current dataset, leaving it intact on cancel/validation error."""
        if self.exporting:return
        folder=filedialog.askdirectory(parent=self.root,title='Change EWCSPH result folder',initialdir=str(Path(self.case.folder).parent),mustexist=True)
        if not folder:return
        try:
            new_case=Case(folder)   # validate FIRST, before touching the old viewer
        except Exception as exc:
            messagebox.showerror('Cannot load selected result folder',str(exc),parent=self.root)
            return
        self._stop_for_reload()
        old_size=self.root.geometry()
        for child in list(self.root.winfo_children()):child.destroy()
        self.root.configure(menu='')
        self.__init__(self.root,folder)
        self.root.geometry(old_size)
        print('[EWCSPH] Changed folder to:',new_case.folder,flush=True)

    def open_new_viewer(self):
        """Run a second independent process; first viewer remains interactive."""
        folder=filedialog.askdirectory(parent=self.root,title='Open new EWCSPH viewer',initialdir=str(Path(self.case.folder).parent),mustexist=True)
        if not folder:return
        try:Case(folder)  # validate input before launching
        except Exception as exc:
            messagebox.showerror('Invalid result folder',str(exc),parent=self.root);return
        try:launch_viewer(__file__,[folder])
        except Exception as exc:messagebox.showerror('New viewer failed',str(exc),parent=self.root)

    def window_settings(self):
        win=tk.Toplevel(self.root);win.title('Viewer window size');win.resizable(False,False)
        win.transient(self.root)
        for j,(name,var) in enumerate((('Width (pixels)',self.window_width),('Height (pixels)',self.window_height))):
            ttk.Label(win,text=name).grid(row=j,column=0,padx=12,pady=8,sticky='w')
            ttk.Entry(win,textvariable=var,width=10).grid(row=j,column=1,padx=12,pady=8)
        def apply():
            try:
                w=int(self.window_width.get());h=int(self.window_height.get())
                if w<700 or h<450:raise ValueError('Minimum window is 700 × 450 pixels')
                self.root.geometry(f'{w}x{h}');win.destroy()
            except ValueError as e:messagebox.showerror('Invalid size',str(e),parent=win)
        ttk.Button(win,text='Apply',command=apply).grid(row=2,column=0,columnspan=2,pady=10)

    def select_tangential_angles(self):
        win=tk.Toplevel(self.root);win.title('Tangential profile angles');win.geometry('475x450')
        win.transient(self.root)
        ttk.Label(win,text='Select angles (degrees); Single mode uses the first checked angle.').pack(anchor='w',padx=9,pady=8)
        grid=ttk.Frame(win);grid.pack(fill='both',expand=True,padx=10)
        for i,(a,var) in enumerate(self.tangential_angles.items()):
            ttk.Checkbutton(grid,text=f'{a}°',variable=var).grid(row=i//4,column=i%4,sticky='w',padx=10,pady=4)
        def finish(): win.destroy();self.draw()
        ttk.Button(win,text='Apply',command=finish).pack(pady=9)

    def save_mp4(self):
        filename=filedialog.asksaveasfilename(parent=self.root,title='Export EWCSPH animation',defaultextension='.mp4',filetypes=[('MP4 video','*.mp4')])
        if not filename:return
        if self.exporting:return
        old_i=self.i;old_tab=self.tabs.select();old_play=self.playing
        self.playing=False;self.exporting=True
        if self.after_id is not None:self.root.after_cancel(self.after_id);self.after_id=None
        self.tabs.select(self.tabmap['Flow Field'])
        dialog=tk.Toplevel(self.root);dialog.title('Exporting animation');dialog.transient(self.root)
        label=ttk.Label(dialog,text='Starting MP4 export…');label.pack(padx=24,pady=18)
        dialog.update()
        try:
            fps=max(1,min(60,int(float(self.fps.get()))))
            # imageio-ffmpeg supplies the encoder; avoid a separate system FFmpeg.
            with imageio.get_writer(filename,fps=fps,codec='libx264',quality=8,macro_block_size=2) as writer:
                for idx in range(len(self.case.frames)):
                    self.i=idx;self.draw();self.canvas.draw()
                    rgba=np.asarray(self.canvas.buffer_rgba())
                    writer.append_data(np.ascontiguousarray(rgba[:,:,:3]))
                    if idx%10==0:
                        label.config(text=f'Exporting frame {idx+1}/{len(self.case.frames)}')
                        dialog.update()
            messagebox.showinfo('MP4 saved',filename,parent=self.root)
        except Exception as e:messagebox.showerror('Export failed',str(e),parent=self.root)
        finally:
            self.exporting=False;self.i=old_i;self.tabs.select(old_tab)
            if dialog.winfo_exists():dialog.destroy()
            self.playing=False;self.play_button.config(text='Play');self.set_frame(old_i)

    def start_line(self):
        self.line_points=[];self.line_active=True
        self.line_label.config(text='Click point A, then point B in Flow Field')
        self.tabs.select(self.tabmap['Flow Field'])
    def clear_line(self):self.line_points=[];self.line_active=False;self.draw()
    def click_flow(self,event):
        if not self.line_active or event.inaxes!=self.ax or event.xdata is None:return
        self.line_points.append((event.xdata,event.ydata))
        if len(self.line_points)==2:
            self.line_active=False;self.line_label.config(text='Line selected (updates with frame)')
            self.tabs.select(self.tabmap['Line Profile'])
        self.draw()
    def set_frame(self,i,from_slider=False):
        j=int(np.clip(i,0,len(self.case.frames)-1))
        if j==self.i and from_slider:return
        self.i=j
        if not from_slider:
            self._suspend_slider=True;self.slider.set(self.i);self._suspend_slider=False
        if not self._suspend_slider:self.draw()
    def toggle_play(self):
        self.playing=not self.playing
        self.play_button.config(text='Pause' if self.playing else 'Play')
        if self.playing:self.play_step()
        elif self.after_id is not None:
            self.root.after_cancel(self.after_id);self.after_id=None;self.draw()
    def play_step(self):
        if not self.playing:return
        if self.i>=len(self.case.frames)-1:
            self.playing=False;self.play_button.config(text='Play');self.draw();return
        self._play_count+=1
        self.set_frame(self.i+max(1,int(self.play_stride.get())))
        try:delay=max(30,int(1000/float(self.fps.get())))
        except ValueError:delay=90
        self.after_id=self.root.after(delay,self.play_step)
    def close(self):
        self.playing=False
        if self.after_id is not None:self.root.after_cancel(self.after_id)
        self.root.destroy()

    def draw(self):
        if not hasattr(self,'canvas') or self._drawing:return
        self._drawing=True
        field=self.field.get();case=self.case
        self.status.configure(text=f'Frame {self.i+1}/{len(case.times)} | h/Δp={case.hdp:g}')
        current=self.tabs.tab(self.tabs.select(),'text')
        try:
            if current=='Flow Field':
                # Keep triangulation, colourbar and artists alive across frames.
                try:step=float(self.spacing.get())
                except ValueError:step=.06
                if field not in self.limits:self.limits[field]=case.initial_stats(field)
                m=draw_field(self.ax,case,self.i,field,self.psi.get(),step,self.near.get(),self.source.get(),self.limits[field],self.show_wall.get(),update_contours=(not self.playing or self.show_contours_play.get()))
                if not hasattr(self,'_colorbar'):
                    self._colorbar=self.fig.colorbar(m,ax=self.ax,label=FIELD_NAMES[field],fraction=.03,pad=.03)
                else:self._colorbar.update_normal(m);self._colorbar.set_label(FIELD_NAMES[field])
                self.ax.set_title((f'SPH  |  ' if self.ax_ref is not None else '') +
                                  fr'$Re={case.Re:g}$   |   $t={case.times[self.i]:.5g}\,\mathrm{{s}}$',pad=8)
                if self.ax_ref is not None:
                    self._draw_potential_reference(field,self.limits[field])
                if not hasattr(self,'line_artist'):
                    (self.line_artist,)=self.ax.plot([],[],'o--',color='lime',lw=1.4)
                if self.line_points:
                    p=np.asarray(self.line_points);self.line_artist.set_data(p[:,0],p[:,1])
                else:self.line_artist.set_data([],[])
                self.canvas.draw_idle()
            else:
                fig,ax,canvas=self.aux[current]
                if current not in ('L2 vs Time','Potential-Flow Validation'):ax.clear()
                if current=='Potential-Flow Validation':
                    self._potential_validation(fig)
                elif current=='L2 vs Time':
                    if not hasattr(fig,'_l2_cursor_lines'):
                        L=case.summary_history();fig.clear();markers=[]
                        for k,(label,idx) in enumerate([('Density',3),('u velocity',1),('v velocity',2),('Velocity vector',0)]):
                            a=fig.add_subplot(2,2,k+1)
                            y=L[:,idx];mask=np.isfinite(y)&(y>0)
                            if mask.any():a.semilogy(case.times[mask],y[mask],lw=1.2)
                            mark=a.axvline(case.times[self.i],color='black',ls='--',lw=.8)
                            markers.append(mark)
                            a.set(xlabel=r'$t$ (s)',ylabel=r'$L_2$ norm / reference deviation',title=label);journal_axes(a)
                            a.grid(alpha=.2)
                        fig.tight_layout();fig._l2_cursor_lines=markers
                    else:
                        for mark in fig._l2_cursor_lines:
                            mark.set_xdata([case.times[self.i]]*2)
                elif current=='Hydrodynamic Forces':
                    self._force_chart(ax)
                elif current=='Local Comparison':
                    self._local_chart(ax)
                elif current=='Tangential vs Radius':
                    d=case.read(self.i);u=case.polar(d['u'].to_numpy())/case.U;v=case.polar(d['v'].to_numpy())/case.U
                    angles=[a for a,vv in self.tangential_angles.items() if vv.get()]
                    if self.tangential_mode.get()=='Single': angles=angles[:1]
                    for ang in angles:
                        k=np.argmin(np.abs(np.angle(np.exp(1j*(case.theta-np.deg2rad(ang))))))
                        ut=-u[k]*np.sin(case.theta[k])+v[k]*np.cos(case.theta[k])
                        line,=ax.plot(np.abs(ut),case.radii-case.wall_r,label=f'SPH {ang}°')
                        if case.has_potential_reference:
                            refu=case.polar(case.reference_values('u'))
                            refv=case.polar(case.reference_values('v'))
                            refut=-refu[k]*np.sin(case.theta[k])+refv[k]*np.cos(case.theta[k])
                            ax.plot(np.abs(refut),case.radii-case.wall_r,'--',color=line.get_color(),
                                    linewidth=1.1,label=f'Potential {ang}°')
                    ax.set(xlabel=r'$|u_\theta|/U_\infty$',ylabel=r'$(r-r_0)/D$',title='Tangential velocity profile')
                    if angles:ax.legend()
                    else:ax.text(.5,.5,'Choose at least one angle',ha='center',transform=ax.transAxes)
                elif current=='Radial Profiles':
                    try:ang=float(self.profile_angle.get())
                    except ValueError:ang=90
                    z=radial_profile(case,self.i,self.profile_field.get(),ang)
                    ax.plot(case.radii-case.wall_r,z,'.-',ms=2,label='SPH')
                    if case.has_potential_reference:
                        # Same angular interpolation as the SPH radial profile.
                        zz=case.polar(case.reference_values(self.profile_field.get()))
                        theta=np.mod(np.deg2rad(ang),2*np.pi)
                        k=int(np.searchsorted(case.theta,theta))
                        left=(k-1)%case.nt;right=k%case.nt
                        al=float(case.theta[left]);ar=float(case.theta[right])
                        if right==0:ar+=2*np.pi
                        if theta<al:theta+=2*np.pi
                        weight=(theta-al)/(ar-al)
                        exact=(1-weight)*zz[left]+weight*zz[right]
                        ax.plot(case.radii-case.wall_r,exact,'--',lw=1.3,label='Potential flow')
                        ax.legend()
                    ax.set(xlabel=r'$(r-r_0)/D$',ylabel=FIELD_NAMES[self.profile_field.get()],title=fr'Radial profile, $\theta={ang:g}^\circ$')
                elif current=='Particle Distribution':
                    ax.scatter(case.x[case.fluid],case.y[case.fluid],s=1,alpha=.65)
                    ax.scatter(case.Xwall,case.Ywall,s=3,c='black');ax.set_aspect('equal')
                    ax.set(xlabel=r'$x/D$',ylabel=r'$y/D$',title='SPH particle distribution')
                elif current=='Line Profile':
                    if len(self.line_points)<2:
                        ax.text(.5,.5,'Select two points on Flow Field to sample a line.',ha='center',transform=ax.transAxes)
                    else:
                        A,B=np.asarray(self.line_points)
                        dist,val=sample_line(case,self.i,field,A,B,n=300)
                        ax.plot(dist,val)
                        ax.set(xlabel=r'$s/D$',ylabel=FIELD_NAMES[field],title='SPH field along A–B')
                if current!='Potential-Flow Validation':
                    journal_axes(ax,grid=True);fig.tight_layout(pad=1.5)
                canvas.draw_idle()
        except Exception as exc:
            messagebox.showerror('Plot error',f'{type(exc).__name__}: {exc}',parent=self.root)
        finally:
            self._drawing=False

def pd_is_numeric(series):
    import pandas as pd
    return pd.api.types.is_numeric_dtype(series)

def main():
    """Robust startup: make folder chooser visible and report errors."""
    import sys
    import traceback
    from tkinter import messagebox
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", nargs="?")
    args, _ = parser.parse_known_args()
    print("[EWCSPH] Starting viewer...", flush=True)
    root = None
    try:
        root = tk.Tk()
        root.title("EWCSPH - Choose simulation results")
        root.geometry("540x100")
        root.lift()
        root.attributes("-topmost", True)
        root.update()
        root.after(500, lambda: root.attributes("-topmost", False) if root.winfo_exists() else None)
        if args.folder:
            folder = args.folder
        else:
            print("[EWCSPH] Opening folder picker. Check Alt+Tab if necessary.", flush=True)
            folder = filedialog.askdirectory(parent=root, title="Select EWCSPH result folder", mustexist=True)
        if not folder:
            print("[EWCSPH] No folder selected; viewer closed.", flush=True)
            root.destroy()
            return
        print(f"[EWCSPH] Loading: {folder}", flush=True)
        root.deiconify()
        SingleViewer(root, folder)
        root.lift()
        root.mainloop()
    except Exception:
        details = traceback.format_exc()
        print("[EWCSPH] STARTUP ERROR:\n" + details, file=sys.stderr, flush=True)
        try:
            if root is not None and root.winfo_exists():
                root.deiconify()
                messagebox.showerror("EWCSPH startup error", details[-4500:], parent=root)
                root.destroy()
        except Exception:
            pass
        raise

if __name__ == "__main__":
    main()
