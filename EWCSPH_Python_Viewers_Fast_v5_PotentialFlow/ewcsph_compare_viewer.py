"""Synchronized multi-case EWCSPH comparison — run: python ewcsph_compare_viewer.py"""
from __future__ import annotations
import argparse
import sys
import subprocess
from pathlib import Path
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import imageio.v2 as imageio
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
from ewcsph_core import Case,FIELD_NAMES
from ewcsph_gui import MultiFolderDialog,add_figure,draw_field,folder_information,journal_axes,launch_viewer

class CompareViewer:
    def __init__(self,root,folders):
        self.root=root;root.title('EWCSPH — Multi-Case Comparison');root.geometry('1450x920')
        self.cases=[Case(folder) for folder in folders]
        self.time_start=max(c.times[0] for c in self.cases)
        self.time_end=min(c.times[-1] for c in self.cases)
        if self.time_end<self.time_start:raise ValueError('Cases have no common saved time interval')
        self.times=np.arange(np.ceil(self.time_start*10-1e-9),np.floor(self.time_end*10+1e-9)+1)*.1
        if not len(self.times):self.times=np.array([self.time_start])
        self.times=self.times[:1201]
        self.i=0;self.playing=False;self.after_id=None;self.exporting=False
        self.field=tk.StringVar(value='vorticity');self.psi=tk.BooleanVar(value=True)
        self.near=tk.BooleanVar(value=True);self.spacing=tk.StringVar(value='0.06')
        self.source=tk.StringVar(value='velocity');self.fps=tk.StringVar(value='10')
        self.layout=tk.StringVar(value='Auto');self.show_wall=tk.BooleanVar(value=True)
        self.xmin=tk.StringVar(value='-10');self.xmax=tk.StringVar(value='10')
        self.ymin=tk.StringVar(value='-1.3');self.ymax=tk.StringVar(value='1.3')
        self.cmin=tk.StringVar(value='');self.cmax=tk.StringVar(value='')
        self.window_width=tk.StringVar(value='1450');self.window_height=tk.StringVar(value='920')
        self.limits={};self._drawing=False;self._suspend_slider=False
        self.show_contours_play=tk.BooleanVar(value=False);self.play_stride=tk.StringVar(value='4')
        self._axes=[];self._shape=None;self._bar=None;self._cax=None;self._settings_after=None
        menus=tk.Menu(root);root.configure(menu=menus)
        filemenu=tk.Menu(menus,tearoff=False);menus.add_cascade(label='File',menu=filemenu)
        filemenu.add_command(label='Change result folders…',command=self.change_folders)
        filemenu.add_command(label='Open new comparison…',command=self.open_new_viewer)
        filemenu.add_separator()
        filemenu.add_command(label='Export PNG',command=self.save_png)
        filemenu.add_command(label='Export MP4',command=self.save_mp4)
        menus.add_command(label='Show full folder names',command=lambda:folder_information(root,[c.folder for c in self.cases]))
        menus.add_command(label='Window size',command=self.window_settings)
        # Three distinct control rows; no overflowing horizontal pack that covers plots.
        top=ttk.Frame(root,padding=(8,5));top.pack(side='top',fill='x')
        ttk.Label(top,text='Field').grid(row=0,column=0,sticky='w',padx=3)
        ttk.Combobox(top,textvariable=self.field,values=list(FIELD_NAMES),state='readonly',width=16).grid(row=0,column=1,padx=3)
        ttk.Checkbutton(top,text='Streamfunction contours',variable=self.psi,command=self.draw).grid(row=0,column=2,padx=8)
        ttk.Label(top,text='Δψ/(U∞D)').grid(row=0,column=3)
        self.spacing_entry=ttk.Entry(top,textvariable=self.spacing,width=7);self.spacing_entry.grid(row=0,column=4,padx=4)
        self.spacing_entry.bind('<Return>',lambda e:self.draw())
        ttk.Checkbutton(top,text='Near-zero levels',variable=self.near,command=self.draw).grid(row=0,column=5,padx=6)
        ttk.Label(top,text='ω source').grid(row=0,column=6,padx=3)
        ttk.Combobox(top,textvariable=self.source,values=['velocity','solver'],state='readonly',width=10).grid(row=0,column=7,padx=3)
        ttk.Checkbutton(top,text='Cylinder',variable=self.show_wall,command=self.draw).grid(row=0,column=8,padx=7)
        ttk.Label(top,text='Layout').grid(row=0,column=9,padx=4)
        ttk.Combobox(top,textvariable=self.layout,values=['Auto','Column','Row','Grid'],state='readonly',width=9).grid(row=0,column=10,padx=3)
        ttk.Button(top,text='Folder key',command=lambda:folder_information(root,[c.folder for c in self.cases])).grid(row=0,column=11,padx=5)
        nav=ttk.Frame(root,padding=(8,5));nav.pack(side='top',fill='x')
        for name,callback in [('|<',lambda:self.set_frame(0)),('<',lambda:self.set_frame(self.i-1))]:
            ttk.Button(nav,text=name,command=callback,width=4).pack(side='left')
        self.playbtn=ttk.Button(nav,text='Play',command=self.toggle,width=7);self.playbtn.pack(side='left',padx=5)
        for name,callback in [('>',lambda:self.set_frame(self.i+1)),('>|',lambda:self.set_frame(len(self.times)-1))]:
            ttk.Button(nav,text=name,command=callback,width=4).pack(side='left')
        self.slider=tk.Scale(nav,from_=0,to=len(self.times)-1,orient='horizontal',showvalue=False,
                             command=lambda x:self.set_frame(int(float(x)),from_slider=True))
        self.slider.pack(side='left',fill='x',expand=True,padx=9)
        ttk.Label(nav,text='FPS').pack(side='left');ttk.Entry(nav,textvariable=self.fps,width=5).pack(side='left',padx=4)
        ttk.Label(nav,text='Frame step').pack(side='left',padx=(8,2))
        ttk.Combobox(nav,textvariable=self.play_stride,values=['1','2','4','8','16'],width=4,state='readonly').pack(side='left')
        ttk.Checkbutton(nav,text='Contours while playing',variable=self.show_contours_play).pack(side='left',padx=7)
        ttk.Button(nav,text='Save MP4',command=self.save_mp4).pack(side='left',padx=5)
        limbar=ttk.Frame(root,padding=(8,5));limbar.pack(side='top',fill='x')
        for tag,var in [('x min',self.xmin),('x max',self.xmax),('y min',self.ymin),('y max',self.ymax),('Colour min',self.cmin),('Colour max',self.cmax)]:
            ttk.Label(limbar,text=tag).pack(side='left',padx=(7,2))
            entry=ttk.Entry(limbar,textvariable=var,width=8);entry.pack(side='left')
            entry.bind('<Return>',lambda e:self.draw())
        ttk.Button(limbar,text='Apply settings',command=self.draw).pack(side='left',padx=8)
        ttk.Button(limbar,text='Window size',command=self.window_settings).pack(side='left')
        ttk.Button(limbar,text='Change Folders',command=self.change_folders).pack(side='left',padx=6)
        ttk.Button(limbar,text='Open New Viewer',command=self.open_new_viewer).pack(side='left',padx=6)
        ttk.Button(limbar,text='Save PNG',command=self.save_png).pack(side='right')
        self.status=ttk.Label(root,text='Loading case geometry…',padding=(10,3))
        self.status.pack(side='bottom',fill='x')
        # Scrollable plotting region: multiple stacked panels never superimpose.
        body=ttk.Frame(root);body.pack(side='top',fill='both',expand=True)
        self.scroll=tk.Canvas(body,highlightthickness=0)
        scrollbar=ttk.Scrollbar(body,orient='vertical',command=self.scroll.yview)
        self.scroll.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right',fill='y');self.scroll.pack(side='left',fill='both',expand=True)
        self.plot_holder=ttk.Frame(self.scroll)
        self.plot_win=self.scroll.create_window((0,0),window=self.plot_holder,anchor='nw')
        self.plot_holder.bind('<Configure>',lambda e:self.scroll.configure(scrollregion=self.scroll.bbox('all')))
        self.scroll.bind('<Configure>',self._resize_scroll)
        self.fig=Figure(figsize=(12,7),dpi=100)
        self.canvas=FigureCanvasTkAgg(self.fig,master=self.plot_holder)
        self.canvas.get_tk_widget().pack(side='top',fill='both',expand=True)
        toolbar=NavigationToolbar2Tk(self.canvas,self.plot_holder,pack_toolbar=False)
        toolbar.update();toolbar.pack(side='bottom',fill='x')
        self.field.trace_add('write',lambda *_:self.draw())
        self.layout.trace_add('write',lambda *_:self.draw())
        self.source.trace_add('write',lambda *_:self.draw())
        self.root.protocol('WM_DELETE_WINDOW',self.close)
        self.draw()

    def _resize_scroll(self,event):
        self.scroll.itemconfigure(self.plot_win,width=max(100,event.width))
        self._update_figure_size()

    def _update_figure_size(self):
        if not self._shape:return
        rows,cols=self._shape
        w=max(760,self.scroll.winfo_width()-20)
        h=max(self.scroll.winfo_height()-20,rows*280+60)
        self.fig.set_size_inches(w/self.fig.dpi,h/self.fig.dpi,forward=False)
        self.canvas.get_tk_widget().configure(width=w,height=h)
        self.scroll.configure(scrollregion=self.scroll.bbox('all'))

    def _choose_folders(self):
        dialog=MultiFolderDialog(self.root)
        return list(dialog.result or [])

    @staticmethod
    def _check_folders(folders):
        cases=[Case(str(p)) for p in folders]
        if not cases:raise ValueError('Select at least one EWCSPH result folder')
        if max(c.times[0] for c in cases)>min(c.times[-1] for c in cases):
            raise ValueError('Selected result folders have no common simulation time interval')
        return cases

    def _stop_for_reload(self):
        self.playing=False
        if self.after_id is not None:
            self.root.after_cancel(self.after_id);self.after_id=None

    def change_folders(self):
        """Reload datasets in-place while preserving the window dimensions."""
        if self.exporting:return
        folders=self._choose_folders()
        if not folders:return
        try:self._check_folders(folders)
        except Exception as exc:
            messagebox.showerror('Invalid comparison folders',str(exc),parent=self.root);return
        old_size=self.root.geometry()
        self._stop_for_reload()
        for child in list(self.root.winfo_children()):child.destroy()
        self.root.configure(menu='')
        self.__init__(self.root,folders)
        self.root.geometry(old_size)
        print('[EWCSPH] Comparison folders changed:',folders,flush=True)

    def open_new_viewer(self):
        """Open a wholly independent comparison process."""
        folders=self._choose_folders()
        if not folders:return
        try:self._check_folders(folders)
        except Exception as exc:
            messagebox.showerror('Invalid comparison folders',str(exc),parent=self.root);return
        try:launch_viewer(__file__,folders)
        except Exception as exc:messagebox.showerror('Could not open new viewer',str(exc),parent=self.root)

    def window_settings(self):
        dlg=tk.Toplevel(self.root);dlg.title('Comparison window size');dlg.transient(self.root)
        for i,(name,var) in enumerate([('Width (pixels)',self.window_width),('Height (pixels)',self.window_height)]):
            ttk.Label(dlg,text=name).grid(row=i,column=0,padx=10,pady=7)
            ttk.Entry(dlg,textvariable=var,width=9).grid(row=i,column=1,padx=10,pady=7)
        def apply():
            try:
                w,h=int(self.window_width.get()),int(self.window_height.get())
                if w<850 or h<500:raise ValueError('Minimum 850 × 500 pixels')
                self.root.geometry(f'{w}x{h}');dlg.destroy()
            except ValueError as exc:messagebox.showerror('Invalid window size',str(exc),parent=dlg)
        ttk.Button(dlg,text='Apply',command=apply).grid(row=2,column=0,columnspan=2,pady=10)

    def set_frame(self,i,from_slider=False):
        j=int(np.clip(i,0,len(self.times)-1))
        if from_slider and j==self.i:return
        self.i=j
        if not from_slider:
            self._suspend_slider=True;self.slider.set(self.i);self._suspend_slider=False
        if not self._suspend_slider:self.draw()
    def toggle(self):
        self.playing=not self.playing
        self.playbtn.config(text='Pause' if self.playing else 'Play')
        if self.playing:self.step()
        elif self.after_id is not None:
            self.root.after_cancel(self.after_id);self.after_id=None;self.draw()
    def step(self):
        if not self.playing:return
        if self.i>=len(self.times)-1:
            self.playing=False;self.playbtn.config(text='Play');self.draw();return
        self.set_frame(self.i+max(1,int(self.play_stride.get())))
        try:delay=max(30,int(1000/float(self.fps.get())))
        except ValueError:delay=100
        self.after_id=self.root.after(delay,self.step)
    def close(self):
        self.playing=False
        if self.after_id is not None:self.root.after_cancel(self.after_id)
        self.root.destroy()
    def draw(self):
        if not hasattr(self,'canvas') or self._drawing:return
        self._drawing=True
        t=self.times[self.i];field=self.field.get();n=len(self.cases)
        self.status.config(text=f'Time {t:.4g} s | {self.i+1}/{len(self.times)} | {n} cases')
        try:
            if field not in self.limits:
                all_lims=[c.initial_stats(field) for c in self.cases]
                self.limits[field]=(min(x[0] for x in all_lims),max(x[1] for x in all_lims))
                if field=='vorticity':
                    maximum=max(abs(x) for x in self.limits[field]);self.limits[field]=(-maximum,maximum)
            lim=self.limits[field]
            cmin=self.cmin.get().strip();cmax=self.cmax.get().strip()
            if cmin and cmax:
                lim=(float(cmin),float(cmax))
                if not np.isfinite(lim).all() or lim[0]>=lim[1]:raise ValueError('Colour limits must be increasing finite numbers')
            elif cmin or cmax:raise ValueError('Supply both colour limits, or leave both empty for automatic scaling')
            step=float(self.spacing.get())
            if not np.isfinite(step) or step<=0:raise ValueError('Contour spacing must be greater than zero')
            xv=(float(self.xmin.get()),float(self.xmax.get()))
            yv=(float(self.ymin.get()),float(self.ymax.get()))
            if xv[0]>=xv[1] or yv[0]>=yv[1]:raise ValueError('Invalid view limits')
            layout=self.layout.get();asp=(xv[1]-xv[0])/(yv[1]-yv[0])
            if layout=='Row' or (layout=='Auto' and asp<.77):nr,nc=1,n
            elif layout=='Grid' or (layout=='Auto' and .77<=asp<=1.3):nc=int(np.ceil(np.sqrt(n)));nr=int(np.ceil(n/nc))
            else:nr,nc=n,1
            shape=(nr,nc)
            if shape!=self._shape:
                self.fig.clear()
                gs=self.fig.add_gridspec(nr,nc,left=.095,right=.89,top=.945,bottom=.075,hspace=.45,wspace=.24)
                axes=np.array([[self.fig.add_subplot(gs[r,c]) for c in range(nc)] for r in range(nr)],dtype=object)
                self._axes=axes;self._shape=shape
                self._cax=self.fig.add_axes([.92,.22,.018,.57]);self._bar=None
                self._update_figure_size()
            else:axes=self._axes
            for j,c in enumerate(self.cases):
                ax=axes[j//nc,j%nc];ii=c.nearest(t)
                draw_field(ax,c,ii,field,self.psi.get(),step,self.near.get(),self.source.get(),lim,self.show_wall.get(),
                           update_contours=(not self.playing or self.show_contours_play.get()))
                ax.set_xlim(xv);ax.set_ylim(yv)
                ax.set_title(fr'Case {j+1}: $Re={c.Re:g}$   |   $t={c.times[ii]:.4g}\,\mathrm{{s}}$   |   $h/\Delta p={c.hdp:g}$',fontsize=9,pad=8)
                ax.set_aspect('equal',adjustable='box');journal_axes(ax)
            for j in range(n,nr*nc):axes[j//nc,j%nc].set_visible(False)
            cmap='RdBu_r' if field=='vorticity' else 'turbo'
            if self._bar is None:
                self._sm=ScalarMappable(norm=Normalize(*lim),cmap=cmap)
                self._bar=self.fig.colorbar(self._sm,cax=self._cax,label=FIELD_NAMES[field])
            else:
                self._sm.set_norm(Normalize(*lim));self._sm.set_cmap(cmap)
                self._bar.update_normal(self._sm);self._bar.set_label(FIELD_NAMES[field])
            self.canvas.draw_idle()
        except Exception as e:
            if not self.playing and not self.exporting:messagebox.showerror('Plot error',f'{type(e).__name__}: {e}',parent=self.root)
            else:print(f'[EWCSPH] Plot error: {type(e).__name__}: {e}')
        finally:self._drawing=False
    def save_png(self):
        filename=filedialog.asksaveasfilename(parent=self.root,title='Export comparison figure',defaultextension='.png',filetypes=[('PNG','*.png')])
        if filename:self.fig.savefig(filename,dpi=180)
    def save_mp4(self):
        filename=filedialog.asksaveasfilename(parent=self.root,title='Save comparison animation',defaultextension='.mp4',filetypes=[('MP4','*.mp4')])
        if not filename or self.exporting:return
        self.exporting=True;previous=self.i;self.playing=False
        if self.after_id is not None:self.root.after_cancel(self.after_id);self.after_id=None
        dialog=tk.Toplevel(self.root);dialog.title('Exporting comparison');dialog.transient(self.root)
        progress=ttk.Label(dialog,text='Preparing MP4…');progress.pack(padx=28,pady=20)
        dialog.update()
        try:
            fps=max(1,min(60,int(float(self.fps.get()))))
            with imageio.get_writer(filename,format='FFMPEG',fps=fps,codec='libx264',quality=8,macro_block_size=2) as writer:
                for idx in range(len(self.times)):
                    self.i=idx;self.draw();self.canvas.draw()
                    rgb=np.ascontiguousarray(np.asarray(self.canvas.buffer_rgba())[:,:,:3])
                    writer.append_data(rgb)
                    if idx%10==0:
                        progress.config(text=f'Writing frame {idx+1}/{len(self.times)}')
                        dialog.update()
            messagebox.showinfo('MP4 saved',filename,parent=self.root)
        except Exception as e:messagebox.showerror('MP4 export failed',str(e),parent=self.root)
        finally:
            self.exporting=False;self.i=previous;self.playing=False;self.playbtn.config(text='Play')
            if dialog.winfo_exists():dialog.destroy()
            self.set_frame(previous);self.draw()

def main():
    """Robust startup: make folder chooser visible and report errors."""
    import sys
    import traceback
    from tkinter import messagebox
    parser = argparse.ArgumentParser()
    parser.add_argument("folders", nargs="*")
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
        if args.folders:
            folders = args.folders
        else:
            print("[EWCSPH] Opening multi-folder picker. Check Alt+Tab if necessary.", flush=True)
            dialog = MultiFolderDialog(root)
            folders = dialog.result
        if not folders:
            print("[EWCSPH] No folders selected; viewer closed.", flush=True)
            root.destroy()
            return
        print("[EWCSPH] Loading folders: " + ", ".join(map(str, folders)), flush=True)
        root.deiconify()
        CompareViewer(root, folders)
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
