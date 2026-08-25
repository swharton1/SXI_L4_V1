#This code will analyse my fits to the GAC images that I rebinned. 
#This is a cleaned up version for sharing. It has had some old functions taken out relating to the GAC challenge.  
#It will just plot the summary plot. 

import numpy as np 
import matplotlib.pyplot as plt 
import os
import glob
import datetime as dt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Circle
import matplotlib.dates as dates 
from matplotlib.ticker import MultipleLocator 
from matplotlib.patches import ConnectionPatch
import string 


from . import read_fits_result 
from . import calc_xgse 
from . import quality_flag as qf 
from . import make_image_axes
from . import read_cmap
from . import get_earth 
from . import subplot_label 

class analyse_gac():
    '''This will read in the GAC fitting results and plot them.''' 
    
    def __init__(self,  path='/data/smile/shared/sims/GAC_V2_Sam/GAC_V2.2_1x1_300s/', tag='original'):
        '''This reads in the data files.
        
        Parameters
        ----------
        path - absolute path of where the fitted fits files are kept. 
        tag - which set of fits to use. 
                Original - Analyse fits as described in the CMEM2 paper (def). 
                https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JA034829
                CTSMAP - Analyse fits that include the background model. 
                CXFOV - Analyse fits with my alternative background subtraction method. 
        
        ''' 
        
        self.path = path 
        self.tag = tag
        filestring=f'SMILE*_fitted_{tag}_cmem2g_mpfit_BATSRUS_hybrid1.fits'
        
        #List the filenames. 
        self.filenames = sorted(glob.glob1(self.path, filestring))
        
        #Loop through each filename and read in the results. 
        self.results = [read_fits_result.results(filename=f, path=path, get_ecube=False, bayane=True) for f in self.filenames] 
        
        #Extract key values. Start with Rcmem. 
        self.rcmem = np.array([r.rcmem for r in self.results]) 
        self.rcmem_err = np.array([r.rcmem_err for r in self.results]) 
        
       
        #Aim Point.
        self.target_loc = np.array([r.target_loc for r in self.results])  
        self.aim = np.array([r.target_loc[0] for r in self.results]) 
        
        #Satellite location in GSE. 
        self.smile_loc = np.array([r.smile_loc for r in self.results]) 
        
        #Get the boundaries of the FOV. 
        self.x_fov = np.array([calc_xgse.get_image(smile_loc=s) for s in self.smile_loc]) 
        self.fov_max = np.array([x.x_gse[-1] for x in self.x_fov])
        self.fov_min = np.array([x.x_gse[0] for x in self.x_fov])
        
        #Time. 
        self.time_string = [r.time_string for r in self.results] 
        self.ut = np.array([dt.datetime.strptime(t, '%Y-%m-%dT%H:%M:%S.%f') for t in self.time_string]) 
        
        #Get optimisation time. 
        self.opt_times = np.array([r.opt_time for r in self.results]) 
        
        #Total count rate of each image. 
        self.total_count_rate = np.array([r.total_count_rate for r in self.results]) 

        #Get total background rate. 
        self.bkg_count_rate = np.array([r.bkg_count_rate for r in self.results]) 
        
        #Extract the resolutions as well for the varied one. 
        self.xres = np.array([r.xdeg_sep for r in self.results]) 
        self.yres = np.array([r.ydeg_sep for r in self.results]) 
        self.expos = np.array([r.exposure for r in self.results]) 
        
        #Extract quality flags. 
        self.qf = np.array([r.qf for r in self.results]) 
        
        #Extract specific parameters. 
        if tag == 'ctsmap2':
            #Get scaling parameter. 
            self.s = np.array([r.params_best[-1] for r in self.results]) 
        
        #Make quality flag with SNR instead of count rate. 
        self.cxfov = [r.original for r in self.results]
        self.bkgmap = [r.bkgmap for r in self.results] 
        
        self.qf2 = [qf.calc_quality_flag_2(self.target_loc[i], self.smile_loc[i], self.cxfov[i], self.bkgmap[i]) for i in range(len(self.qf))]
        
        #Extract Cusp ID. 
        self.cusp = [r.cusp for r in self.results] 
        
        #Convert to number for plotting. 
        self.cusp_no = [1 if r else 0 for r in self.cusp] 
        

    def plot_results(self, cmap='lundi', vmin=0, vmax=150, i0=50, isep=80, save=False, add_qf=False, add_cusp=False, figname=None):
        '''This will create a large multipanel plot showing the results throughout the orbit. 
        
        Parameters
        ----------
        cmap - def = 'lundi' for the sample images.
        vmin - minimum counts in images. def = 0. 
        vmax - maximum counts in images. def = 20.
        i0 - which image in the sequence to start with. def = 6. 
        isep - how often to plot an image. def = every sixth. 
        save - boolean to save the plot to a standard name. 
        add_qf - boolean to add the quality flag plot. 
        add_cusp - boolean to add Bayane's cusp ID result. 
        figname - option to manually specify the figure name. 
        
        ''' 
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()
        
        #Sort scaling out for colour maps. 
        #This makes it counts per degree squared instead of per pixel.  
        scaling = self.xres*self.yres
        

            
        #Make a figure. 
        plt.close("all")
        fig = plt.figure(figsize=(6,8))
        if add_qf:
            gs = GridSpec(6,6,figure=fig) 
        else:
            gs = GridSpec(5,6,figure=fig) 
        
        
        fig.subplots_adjust(hspace=0.3, top=0.90, wspace=0.4)     
        
        letters = string.ascii_lowercase 
        
        #Make some subplots at the top for images. 
        ax1 = fig.add_subplot(gs[0, 0]) 
        ax2 = fig.add_subplot(gs[0, 1]) 
        ax3 = fig.add_subplot(gs[0, 2]) 
        ax4 = fig.add_subplot(gs[0, 3]) 
        ax5 = fig.add_subplot(gs[0, 4]) 
        ax6 = fig.add_subplot(gs[0, 5]) 
        image_ax = [ax1, ax2, ax3, ax4, ax5, ax6]
        
        #Make a wide subplot for key results. 
        ax7 = fig.add_subplot(gs[1:3,:])
        ax7b = fig.add_subplot(gs[3,:])
        
        #Make a wide subplot for the quality flag.  
        if add_qf: 
            ax8 = fig.add_subplot(gs[4,:])
        
            #Make three subplots to showcase the orbit. 
            ax9 = fig.add_subplot(gs[5,0:2])
            ax10 = fig.add_subplot(gs[5,2:4])
            ax11 = fig.add_subplot(gs[5,4:6])     
        
        else:
            #Make three subplots to showcase the orbit. 
            ax9 = fig.add_subplot(gs[4,0:2])
            ax10 = fig.add_subplot(gs[4,2:4])
            ax11 = fig.add_subplot(gs[4,4:6])  
            
            
        image_no = np.arange(6)*isep + i0

        
        for i, ival in enumerate(image_no):
            array = self.results[ival].original
            xdeg_min = self.results[ival].xdeg_min
            ydeg_min = self.results[ival].ydeg_min 
            n_pixels = self.results[ival].n_pixels
            m_pixels = self.results[ival].m_pixels
            label = f'({letters[i]})'#+self.ut[ival].strftime('%Y-%m-%d\n%H:%M')
            
            #Make axis one. 
            make_image_axes.make_image_axes(image_ax[i], array/scaling[ival], xdeg_min, ydeg_min, n_pixels, m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title=None, add_cbar=False)
            image_ax[i].set_xticks([])
            image_ax[i].set_yticks([])
            image_ax[i].set_ylabel('')
            image_ax[i].set_xlabel('')
            image_ax[i].set_title(label, fontsize=10) 
        

        
        
        #Add CMEM results and FOV info. 
        ax7.errorbar(self.ut, self.rcmem, yerr=self.rcmem_err, c='k', capsize=1, lw=0.5)
        ax7.plot(self.ut, self.aim, c='grey')
        ax7.plot(self.ut, self.fov_max, c='grey', lw=0.1)
        ax7.plot(self.ut, self.fov_min, c='grey', lw=0.1)
        ax7.set_xlim(self.ut[0], self.ut[-1]) 
        ax7.set_ylabel(r'X$_{GSE}$ [R$_E$]')
        ax7.set_ylim(4,16)
        ax7.yaxis.set_major_locator(MultipleLocator(1))
        
        #Sort out polygons to show the edges of the FOV. 
        xlims = ax7.get_xlim()
        time = np.linspace(xlims[0], xlims[1], self.ut.size) 
        xvals = list(time)
        xvals.append(xlims[1])
        xvals.append(xlims[0]) 
        
        yvals = list(self.fov_min)
        yvals.append(0)
        yvals.append(0) 
        
        yvals2 = list(self.fov_max)
        yvals2.append(30)
        yvals2.append(30)
        
        #Bottom Polygon
        verts = [[xvals[i],yvals[i]] for i in range(len(xvals))]
        polygon = Polygon(verts, closed=True, edgecolor='grey', facecolor='grey', alpha=0.5) 
        ax7.add_patch(polygon) 
        
        #Top Polygon 
        verts2 = [[xvals[i],yvals2[i]] for i in range(len(xvals))]
        polygon2 = Polygon(verts2, closed=True, edgecolor='grey', facecolor='grey', alpha=0.5) 
        ax7.add_patch(polygon2) 
        
        #Add vertical lines to show locations of images.
        ylims = ax7.get_ylim()
        for i, ival in enumerate(image_no):
            ax7.plot([self.ut[ival],self.ut[ival]], ylims, c='purple', lw=1, alpha=1, zorder=3)
            ax7.set_ylim(ylims) 
        
        if add_cusp: 
            
            ax7.plot(self.ut, self.cusp_no+ax7.get_ylim()[1], clip_on=False, c='b') 
            fig.text(0.04, 0.78, 'Cusp', fontsize=8) 
        
        #Add connecting lines between the images and their point in time on the main axis. 
        for i, ival in enumerate(image_no):
            xyA = (0, -13.5)
            xyB = (self.ut[ival], ax7.get_ylim()[1])
        
            con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA="data", coordsB="data",
                      axesA=image_ax[i], axesB=ax7, color="purple")
            image_ax[i].add_artist(con)
        
        
        
        
            
        subplot_label.add_label(ax7, text=f'({letters[6]})', backcolour='w', textcolour='k', width=0.05, zorder=2, fontsize=10, alpha=0.8) 
        
        
        #Add Count Rate plot. total_count_rate is just for the foreground. 
        ax7b.plot(self.ut, self.total_count_rate, c='b', label='CXFOV') 
        ax7b.plot(self.ut, self.bkg_count_rate, c='navy', label='BKG')
        ax7b.set_xlim(self.ut[0], self.ut[-1]) 
        ax7b.set_ylabel('Count Rate '+r'[s$^{-1}$]')
        ax7b.yaxis.set_major_locator(MultipleLocator(100))
        ax7b.set_ylim(0,)
        ax7b.legend(fontsize=8)
        subplot_label.add_label(ax7b, text=f'({letters[7]})', backcolour='w', textcolour='k', width=0.05, height=0.2, zorder=2, fontsize=10, alpha=0.8) 
        

        
        #Add vertical lines to show locations of images.
        ylims = ax7b.get_ylim()
        for i, ival in enumerate(image_no):
            ax7b.plot([self.ut[ival],self.ut[ival]], ylims, c='purple', lw=1, alpha=1)
            ax7b.set_ylim(ylims)     
            
        #Sort time axes. 
        t_form = dates.DateFormatter('%H:%M')
        ax7.xaxis.set_major_formatter(t_form)
        ax7b.xaxis.set_major_formatter(t_form)
        ax7.grid()
        ax7b.grid()
        
        #Now add the quality flags. 
        if add_qf:
        
            #Replace with new QF decomposition plot. 
            ax8 = qf.interpret_qf_list(self.ut, self.qf2, ax1=ax8) 
            ax8.set_ylabel('QF2', color='g') 
            subplot_label.add_label(ax8, text=f'({letters[8]})', backcolour='w', textcolour='k', width=0.05, zorder=2, fontsize=10, alpha=0.8) 
            
            #Redo y axis labels here as they don't fit. Put descriptions at bottom. 
            yticklabels = [r'$\lambda^{a}$', r'$\lambda^{b}$', r'SNR$^{a}$', r'SNR$^{b}$']
            ax8.set_yticklabels(yticklabels)            
            
            fig.text(0.1, 0.03, r'$\lambda^{a}$: 25$^{\circ}$ < $\lambda$ $\leq$ 50$^{\circ}$', ha='left', fontsize=8)
            fig.text(0.3, 0.03, r'$\lambda^{b}$: $\lambda$ > 50$^{\circ}$', ha='left', fontsize=8)
            fig.text(0.45, 0.03, r'SNR$^{a}$: 0.33 < SNR $\leq$ 0.67', ha='left', fontsize=8)
            fig.text(0.75, 0.03, r'SNR$^{b}$: SNR < 0.33', ha='left', fontsize=8)
            #Add vertical lines to show locations of images.
            ylims = ax8.get_ylim()
            for i, ival in enumerate(image_no):
                ax8.plot([self.ut[ival],self.ut[ival]], ylims, c='purple', lw=1, alpha=1)
                ax8.set_ylim(ylims) 
            
        #Now sort out the orbit axes. 
        xpos = self.smile_loc[:,0]
        ypos = self.smile_loc[:,1]
        zpos = self.smile_loc[:,2]
        
        ax9.plot(xpos, zpos, c='k')
        ax10.plot(xpos, ypos, c='k')
        ax11.plot(ypos, zpos, c='k') 
        
        get_earth.make_earth(ax9, rotation=-90)
        get_earth.make_earth(ax10, rotation=-90)
        circle = Circle((0,0), 1, facecolor='w', edgecolor='navy')
        ax11.add_patch(circle)
        
        ax9.set(xlabel=r'X$_{GSE}$', ylabel=r'Z$_{GSE}$', xlim=(-15,15), ylim=(-5,25))
        ax10.set(xlabel=r'X$_{GSE}$', ylabel=r'Y$_{GSE}$', xlim=(-15,15), ylim=(-15,15))
        ax11.set(xlabel=r'Y$_{GSE}$', ylabel=r'Z$_{GSE}$', xlim=(-15,15), ylim=(-5,25))
        
        ax9.grid()
        ax10.grid()
        ax11.grid()
        
        if add_qf: 
            subplot_label.add_label(ax9, text=f'({letters[9]})', backcolour='w', textcolour='k', width=0.2, height=0.2, zorder=2, fontsize=10, alpha=0.8) 
            subplot_label.add_label(ax10, text=f'({letters[10]})', backcolour='w', textcolour='k', width=0.2, height=0.2, zorder=2, fontsize=10, alpha=0.8) 
            subplot_label.add_label(ax11, text=f'({letters[11]})', backcolour='w', textcolour='k', width=0.2, height=0.2, zorder=2, fontsize=10, alpha=0.8) 
        else:
            subplot_label.add_label(ax9, text=f'({letters[8]})', backcolour='w', textcolour='k', width=0.2, height=0.2, zorder=2, fontsize=10, alpha=0.8) 
            subplot_label.add_label(ax10, text=f'({letters[9]})', backcolour='w', textcolour='k', width=0.2, height=0.2, zorder=2, fontsize=10, alpha=0.8) 
            subplot_label.add_label(ax11, text=f'({letters[10]})', backcolour='w', textcolour='k', width=0.2, height=0.2, zorder=2, fontsize=10, alpha=0.8) 
            
            
        ax9.set_aspect('equal')
        ax10.set_aspect('equal')
        ax11.set_aspect('equal')
        
        title=f'GAC Results with CMEM2: {self.tag.capitalize()}'
        fig.text(0.5, 0.95, title, ha='center', fontsize=12) 
        
        if save: 
            qftag = '_qf_dec' if add_qf else ''
            cusptag = '_cusp' if add_cusp else '' 
            figname = f'plots/GAC_results_CMEM2_{self.tag}{qftag}{cusptag}.png'
            print ('Saving... ', figname)
            fig.savefig(figname)
    

