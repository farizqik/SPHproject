"""Reusable Tk/Matplotlib widgets for EWCSPH Python viewers."""
from __future__ import annotations
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from ewcsph_core import FIELD_NAMES
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter

rcParams.update({
    "font.family":"DejaVu Sans", "font.size":10, "axes.labelsize":11,
    "axes.titlesize":11, "axes.linewidth":0.8, "xtick.direction":"in",
    "ytick.direction":"in", "xtick.minor.visible":True, "ytick.minor.visible":True,
    "legend.frameon":False, "savefig.dpi":220
})

def journal_axes(ax, grid=False):
    """Publication-style axes with mathematics typeset in labels."""
    ax.tick_params(which="major",direction="in",top=True,right=True,length=5,width=.8)
    ax.tick_params(which="minor",direction="in",top=True,right=True,length=2.5,width=.6)
    for axis in (ax.xaxis,ax.yaxis):
        formatter=ScalarFormatter(useOffset=False,useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-3,3))
        axis.set_major_formatter(formatter)
    if grid:ax.grid(True,alpha=.16,lw=.55)
    else:ax.grid(False)
    for spine in ax.spines.values():spine.set_linewidth(.8)

def pretty_column(name):
    """Human-readable mathematical notation for raw C++ CSV columns."""
    name=str(name)
    low=name.lower().replace("_","").replace(" ","")
    if low in ("thetadeg","angledeg","angle","theta","beta","betadeg"):return r"$\theta$ (°)"
    if low in ("t","time","times"):return r"$t$ (s)"
    if low.startswith("cpsph") or low=="cp":return r"$C_p$"
    if low.startswith("cd"):return r"$C_D$"
    if low.startswith("cl"):return r"$C_L$"
    if low.startswith("cf"):return r"$C_F$"
    if low in ("usph","uxsph","u"):return r"$u$ (m/s)"
    if low in ("vsph","vysph","v"):return r"$v$ (m/s)"
    if low.startswith("ursph") or low.startswith("uradial"):return r"$u_r$ (m/s)"
    if low.startswith("utsph") or low.startswith("utangential"):return r"$u_\theta$ (m/s)"
    if low.startswith("pressure") or low.startswith("psph"):return r"$p$ (Pa)"
    if low.startswith("timestep"):return "Time step"
    if low.startswith("vorticity"):return r"$\omega$ (s$^{-1}$)"
    return name

def launch_viewer(script,folders):
    """Launch a truly independent GUI even when the first was run in Spyder."""
    import subprocess,sys
    from pathlib import Path
    subprocess.Popen([sys.executable,str(Path(script).resolve()),*[str(x) for x in folders]],
                     cwd=str(Path(script).resolve().parent),
                     creationflags=getattr(subprocess,"CREATE_NEW_PROCESS_GROUP",0))

FIELDS=list(FIELD_NAMES)

def folder_information(root,folders):
    dlg=tk.Toplevel(root);dlg.title('EWCSPH — Full simulation folder paths');dlg.geometry('900x290')
    ttk.Label(dlg,text='Full simulation folder paths (select and copy):').pack(anchor='w',padx=12,pady=(12,5))
    txt=tk.Text(dlg,wrap='none',height=9)
    txt.pack(fill='both',expand=True,padx=12)
    txt.insert('1.0','\n'.join(f'{i+1}. {Path(f).resolve()}' for i,f in enumerate(folders)))
    txt.configure(state='normal')
    def copy():
        dlg.clipboard_clear();dlg.clipboard_append('\n'.join(str(Path(f).resolve()) for f in folders))
    controls=ttk.Frame(dlg);controls.pack(fill='x',padx=12,pady=8)
    ttk.Button(controls,text='Copy paths',command=copy).pack(side='left')
    ttk.Button(controls,text='Hide',command=dlg.withdraw).pack(side='right')
    return dlg

def levels_for(psi,spacing=0.06,near=True):
    if spacing<=0: raise ValueError('Streamfunction spacing must be positive')
    finite=psi[np.isfinite(psi)]
    if finite.size==0:return np.array([])
    mn,mx=finite.min(),finite.max()
    n=min(500,int(np.ceil(max(abs(mn),abs(mx))/spacing))+1)
    lev=np.arange(-n,n+1)*spacing
    if near:
        extra=np.array([.0005,.001,.002,.004,.008,.015,.03,.05])
        lev=np.unique(np.r_[lev,-extra,extra])
    return lev[(lev>=mn)&(lev<=mx)]

def draw_field(ax,case,index,field='vorticity',show_stream=True,spacing=.06,
               near_levels=True,source='velocity',limits=None,show_wall=True,
               update_contours=True):
    """Persistent mesh: update CData only, reuse triangulation and contours."""
    from matplotlib.patches import Circle
    if not hasattr(case,'triangles'):
        jj=np.arange(case.nr-1);kk=np.arange(case.nt)
        here=(kk[:,None]*case.nr+jj[None,:]).ravel()
        nxt=(((kk+1)%case.nt)[:,None]*case.nr+jj[None,:]).ravel()
        case.triangles=np.vstack([np.column_stack([here,here+1,nxt]),
                                   np.column_stack([here+1,nxt+1,nxt])])
    state=getattr(ax,'_ewcsph_state',None)
    cmap='RdBu_r' if field=='vorticity' else 'turbo'
    z=case.polar(case.values(index,field)).ravel()
    if state is None or state['case'] is not case:
        ax.clear()
        mesh=ax.tripcolor(case.X.ravel(),case.Y.ravel(),z,triangles=case.triangles,
                          shading='gouraud',cmap=cmap,rasterized=True)
        wall=Circle((0,0),case.wall_r,facecolor='#292d34',edgecolor='black',lw=.4,zorder=8)
        ax.add_patch(wall)
        state={'case':case,'mesh':mesh,'wall':wall,'contour':None,'contour_key':None}
        ax._ewcsph_state=state
        ax.set_aspect('equal',adjustable='box');ax.set_xlabel(r'$x/D$');ax.set_ylabel(r'$y/D$');journal_axes(ax)
    else:
        mesh=state['mesh'];mesh.set_array(z);mesh.set_cmap(cmap)
    if limits is not None:mesh.set_clim(*limits)
    state['wall'].set_visible(bool(show_wall and len(case.wall)))
    key=(int(index),float(spacing),bool(near_levels),str(source)) if show_stream else None
    if state['contour'] is not None and not show_stream:
        state['contour'].set_visible(False)
    elif show_stream and case.has_wall:
        if state['contour_key'] != key and update_contours:
            if state['contour'] is not None:state['contour'].remove()
            state['contour']=None
            psi=case.streamfunction(index,source)
            if psi is not None:
                XX,YY,PP=case.contour_geometry(psi)
                lev=levels_for(PP,spacing,near_levels)
                if len(lev)>1:
                    state['contour']=ax.contour(XX,YY,PP,levels=lev,colors='black',linewidths=.55)
            state['contour_key']=key
        if state['contour'] is not None:state['contour'].set_visible(True)
    return mesh

def add_figure(parent,fig=None,toolbar=True):
    fig=fig if fig is not None else Figure(figsize=(10,5),dpi=95)
    canvas=FigureCanvasTkAgg(fig,master=parent)
    canvas.get_tk_widget().pack(fill='both',expand=True)
    if toolbar:
        bar=NavigationToolbar2Tk(canvas,parent,pack_toolbar=False)
        bar.update();bar.pack(side='bottom',fill='x')
    return fig,canvas

class MultiFolderDialog(tk.Toplevel):
    """Direct multi-directory browser with Ctrl/Shift selection, plus an Add queue.
    Unlike OS file pickers, lets user select multiple DIRECTORIES, even on
    platforms where the native directory picker is single-select only.
    """
    def __init__(self,parent):
        super().__init__(parent);self.title('Select EWCSPH result folders (Ctrl / Shift)')
        self.geometry('870x590');self.transient(parent);self.result=[]
        self.current=Path.home();self.queued=[]
        top=ttk.Frame(self);top.pack(fill='x',padx=10,pady=10)
        ttk.Label(top,text='Folder:').pack(side='left')
        self.path_var=tk.StringVar(value=str(self.current))
        entry=ttk.Entry(top,textvariable=self.path_var);entry.pack(side='left',expand=True,fill='x',padx=8)
        entry.bind('<Return>',lambda e:self.go())
        ttk.Button(top,text='Go',command=self.go).pack(side='left')
        ttk.Button(top,text='Up',command=lambda:self.enter(self.current.parent)).pack(side='left',padx=4)
        ttk.Button(top,text='Browse…',command=self.pick).pack(side='left')
        ttk.Label(self,text='Select result directories with Ctrl/Shift. Double-click to enter a directory.').pack(anchor='w',padx=12)
        frm=ttk.Frame(self);frm.pack(fill='both',expand=True,padx=10,pady=5)
        self.box=tk.Listbox(frm,selectmode=tk.EXTENDED,exportselection=False,font=('Segoe UI',10))
        scroll=ttk.Scrollbar(frm,orient='vertical',command=self.box.yview)
        self.box.configure(yscrollcommand=scroll.set)
        self.box.pack(side='left',fill='both',expand=True);scroll.pack(side='right',fill='y')
        self.box.bind('<Double-Button-1>',self.open_selected)
        ttk.Label(self,text='Selected / added folders:').pack(anchor='w',padx=12)
        self.chosen=tk.Listbox(self,height=4,font=('Segoe UI',9));self.chosen.pack(fill='x',padx=10,pady=4)
        bottom=ttk.Frame(self);bottom.pack(fill='x',padx=10,pady=9)
        ttk.Button(bottom,text='Add selected',command=self.add).pack(side='left')
        ttk.Button(bottom,text='Remove queued',command=self.remove).pack(side='left',padx=6)
        ttk.Button(bottom,text='Use selected / added folders',command=self.confirm).pack(side='right')
        ttk.Button(bottom,text='Cancel',command=self.destroy).pack(side='right',padx=8)
        self.entries=[];self.enter(self.current)
        self.grab_set();self.wait_window(self)

    def pick(self):
        folder=filedialog.askdirectory(parent=self,title='Navigate to folder',initialdir=str(self.current))
        if folder:self.enter(Path(folder))
    def go(self):self.enter(Path(self.path_var.get()).expanduser())
    def enter(self,p):
        p=Path(p)
        if not p.is_dir():messagebox.showerror('Folder not found',str(p),parent=self);return
        try:items=sorted((f for f in p.iterdir() if f.is_dir()),key=lambda x:x.name.casefold())
        except PermissionError:items=[]
        self.current=p;self.path_var.set(str(p));self.entries=items;self.box.delete(0,'end')
        for item in items:
            mark='[RESULT] ' if (item/'SimulationConstants.csv').is_file() else '         '
            self.box.insert('end',mark+item.name)
    def open_selected(self,event=None):
        ind=self.box.curselection()
        if ind:self.enter(self.entries[ind[0]])
    def add(self):
        selection=[self.entries[i] for i in self.box.curselection()]
        if not selection and (self.current/'SimulationConstants.csv').is_file():selection=[self.current]
        invalid=[p for p in selection if not (p/'SimulationConstants.csv').is_file()]
        for p in selection:
            if p not in invalid and p not in self.queued:self.queued.append(p)
        self.chosen.delete(0,'end')
        for p in self.queued:self.chosen.insert('end',str(p))
        if invalid:messagebox.showwarning('Not a result folder',f'{len(invalid)} folder(s) lack SimulationConstants.csv.',parent=self)
    def remove(self):
        idx=set(self.chosen.curselection());self.queued=[p for i,p in enumerate(self.queued) if i not in idx]
        self.chosen.delete(0,'end')
        for p in self.queued:self.chosen.insert('end',str(p))
    def confirm(self):
        self.add()
        if not self.queued:messagebox.showinfo('Select folders','Select one or more result folders.',parent=self);return
        self.result=list(self.queued);self.destroy()
