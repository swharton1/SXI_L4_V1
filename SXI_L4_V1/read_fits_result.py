#This will read in the fitted files produced by fit_sxi_sim so that analysis can be done afterwards. 
#There is also another class at the bottom that does animations. 

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable 
import pdb
import matplotlib.animation as animation 
import glob 

#from .SXI_Core import read_fits_cube
from . import make_image_axes
from . import overlays 
from . import mpause_overlay
from . import read_cmap
from . import quality_flag as qf
from . import bayane_cusp 



class results():
    '''This class contains all the information in a useable format for the fitted data FITS files. Adds it all to an object to be used by other analysis programs.''' 
    
    def __init__(self, filename='SMILE_SXI_L3_SCIM60x60-SCI-CXF_20260318T1910-20260318T1915_V01_fitted_cxfov_cmem2g_mpfit_BATSRUS_hybrid1.fits', path='/data/smile/shared/sims/GAC_V2_Sam/GAC_V2.2_1x1_300s/', counts_file=True, get_ecube=False, filetype='PPMLR', sw_info=False, bayane=True):
        '''This read in the FITS file and opens it.
        
        Parameters
        ----------
        filename - FITS file of the fitted image. 
        path - Path to the FITS file. 
        counts_file - Boolean. If true, FITS file is a full sxi counts map fit produced by fit_image_sxi. Else, FITS file is just an intensity fit produced by fit_image.py. 
        get_ecube - boolean to read in the emissivity cube used to simulate the original image. def = True. 
        filetype - def = 'PPMLR'. Needed when reading in the emissivity cube. 
        sw_info - boolean to read in solar wind info. 
        bayane - boolean to run Bayane's cusp ID method. 
        
        '''
        
        self.path = path 
        self.filename = filename 
        self.counts_file = counts_file 
        
        #Open the FITS file. 
        # Check file exists. If it does, open it. 
        try:
            print ('Read {}...'.format(self.filename)) 
            # Open the FITS file with this command. 
            with fits.open(self.path+self.filename) as hdul:
                self.hdul = hdul 
                
                #Get all extension names. 
                ext_names = [hdul[i].name for i in range(len(hdul))]
                
                # Add the headers and data to the object. 
                self.primary_header = self.hdul['PRIMARY'].header
                
                #Get the fitted image, original image, and model parameters. 
                self.fitted = self.hdul['PRIMARY'].data 
                
                if 'CXFOV' in ext_names:
                    self.original = self.hdul['CXFOV'].data 
                elif 'CTSMAP' in ext_names:
                    self.original = self.hdul['CTSMAP'].data
                else:
                    self.original = self.hdul['ORIGINAL'].data 
                
                if 'BKGMAP' in ext_names:
                    self.bkgmap = self.hdul['BKGMAP'].data 
                    
                self.param_table = self.hdul['PARAMETERS'].data
                
                #This table contains all the information on how the 
                #parameters and cost value varied throughout the fitting process. 
                self.p_variation = self.hdul['P_VARIATION'].data
                
                #Extract covariance matrix. Won't work for nmd results. 

                try:
                    self.covar = self.hdul['COVAR'].data
                except KeyError:
                    print ('No covariance found') 
                
                #Extract important parameter info from the table. 
                self.param_names = self.param_table['Parameter']
                self.params0 = self.param_table['Initial']
                self.params_best = self.param_table['Optimum']
                
                #Extra extensions for fit_image_sxi.py 
                if counts_file:
                    #self.los_intensity = self.hdul['PRIMARY'].data 
                    self.error = self.hdul['ERROR'].data
                
                    #Get exposure time. 
                    self.exposure = self.primary_header['EXPOS']
                    self.params_error = self.param_table['Error']
                    
                #Extract all information from the header. 
                self.m_pixels = self.primary_header['NAXIS1']
                self.n_pixels = self.primary_header['NAXIS2']
                self.image_filename = self.primary_header['FILENAME']
                self.smile_loc = np.array([self.primary_header['POS_X'], self.primary_header['POS_Y'], self.primary_header['POS_Z']])
                self.target_loc = np.array([self.primary_header['AIM_X'], self.primary_header['AIM_Y'], self.primary_header['AIM_Z']])
                
                #Pixel widths. 
                self.xdeg_sep = self.primary_header['CDELT1']
                self.ydeg_sep = self.primary_header['CDELT2']
                
                #Lower bounds. 
                self.xdeg_min = self.primary_header['CRVAL1']
                self.ydeg_min = self.primary_header['CRVAL2']
                
                #Calculate FOV. 
                self.phi_fov = -2*self.xdeg_min 
                self.theta_fov = -2*self.ydeg_min   
                
                #Solar wind information. 
                #Get information from the primary header about the solar wind. 
                if sw_info: 
                    self.density = self.primary_header['DENS_SW']
                    self.vx = self.primary_header['VX']
                    self.vy = self.primary_header['VY']
                    self.vz = self.primary_header['VZ']
                    self.bx = self.primary_header['BX']
                    self.by = self.primary_header['BY']
                    self.bz = self.primary_header['BZ'] 
                    self.temp = self.primary_header['TEMP_SW']
                    self.flux = abs(self.primary_header['FLUX']) 
                    self.pdyn = self.primary_header['PDYN']
                    self.pmag = self.primary_header['PMAG']
                self.dipole = self.primary_header['DIPOLE']
                
                if hasattr(self.primary_header, 'HISTORY'):
                    self.history = self.primary_header['HISTORY']
                self.model = self.primary_header['MODEL']
                self.fit_method = self.primary_header['FIT']
                
                #Get rcmem and the fitted parameters.
                self.rcmem = self.primary_header['RCMEM'] 
                self.rcmem_err = self.primary_header['R_ERR']
                
                #Get fitting info. 
                self.minimum_cost = self.primary_header['MIN_COST']
                self.iterations = self.primary_header['ITER']
                self.opt_time = self.primary_header['OPT_TIME']    
                
                
                #Get time. 
                self.time_string = self.primary_header['DATE-OBS']
                #self.time_string = '' #Temporary. 
               
                self.p_spacing = self.primary_header['P_SPACE']                
                self.p0_method = self.primary_header['P0_METH']
                self.total_count_rate = self.primary_header['CNT_RATE'] 
                
                #Calculate bkg count rate. 
                self.bkg_count_rate = self.bkgmap.sum()/self.exposure 
                #Get any smoothed information. 
                #if counts_file:
                    #self.smooth = self.primary_header['SMOOTH']
                    #self.smooth_sigma = self.primary_header['SM_SIGMA']
                    #if self.smooth:
                    #    self.smoothed_image = self.hdul['SMOOTHED'].data
                
                #if get_ecube:   
                    #Get max Lx from ecube. 
                    #Get ecube filename. 
                #    self.get_ecube_name() 
        
                    #Read in the ecube. 
                #    print ('Read... ', self.ecube_fname)
                #    ecube = read_fits_cube.read_fits_cube(filename=self.ecube_fname, filetype=filetype)  
                
                #    self.maxIx = ecube.maxIx
                #    self.f25 = ecube.f
                #    self.maxdIx = ecube.maxdIx
                #else:
                self.maxIx = None 
                
                #Get the image unit vectors for magnetopause projection. 
                self.get_fov()    
                
                try: 
                    self.qf = self.primary_header['QF'] 
                except KeyError: 
                    try: 
                        self.qf = qf.calc_quality_flag(self.target_loc, self.smile_loc, self.exposure, self.hdul['CXFOV'].data)
                    except KeyError: 
                        #If cxfov is not available, estimate it with this calculation. Don't have the vignetting info unfortunately. 
                        self.qf = qf.calc_quality_flag(self.target_loc, self.smile_loc, self.exposure, self.hdul['CTSMAP'].data - self.hdul['BKGMAP'].data)
                    #raise ValueError('No QF flag')     
                try: 
                    self.qf2 = self.primary_header['QF2'] 
                except KeyError:       
                    self.qf2 = qf.calc_quality_flag_2(self.target_loc, self.smile_loc, self.hdul['CXFOV'].data, self.hdul['BKGMAP'].data)     
                
                #Bayane's cusp routine. 
                self.bayane = bayane 
                if bayane: 
                    self.cusp = bayane_cusp.test2_is_cusp(self.hdul['CXFOV'].data) 
                    
        except (FileNotFoundError, IOError):
            print ("Filename not found: {}".format(self.path+self.filename))
        
        
            
            
    def get_ecube_name(self):
        '''This uses the history and does some checks to make sure we have the right FITS file.''' 
        
        #Get the ecube from the history. 
        #Make sure you just have the filename without the path. 
        self.ecube_fname = str(self.history).split('/')[-1] 
        
        #Check for .dat and replace with .fits if needed. 
        split_name = self.ecube_fname.split('.')
        if split_name[-1] == 'dat':
            self.ecube_fname = split_name[0]+'.fits'
            
    def get_fov(self):
        '''This will use the aim and camera information to calculate 
        the unit vectors for the image. You would need this for the magnetopause projections.  
        '''
        
        #Get magnitude of spacecraft vector - its radial position. 
        self.smag = np.sqrt(self.smile_loc[0]**2 + self.smile_loc[1]**2 + self.smile_loc[2]**2)
        
        #Get alpha angle. 
        self.alpha_angle = np.arctan2(self.smile_loc[0],np.sqrt(self.smile_loc[1]**2 + self.smile_loc[2]**2))
        
        #Get earth angle. 
        self.r_angle = np.arcsin(1/self.smag) 

        #Get Look vector. 
        self.L = self.target_loc - self.smile_loc 
        
        #Get magnitude of Look vector. 
        self.Lmag = np.sqrt(self.L[0]**2 + self.L[1]**2 + self.L[2]**2)
        
        #Get unit vector for look directions. 
        self.L_unit = self.L/self.Lmag
        
        #Get angle between look vector and smile vector. This should be 20.3 + r_angle.  
        self.limb_c = np.arccos(np.dot(-self.smile_loc, self.L)/(self.smag*self.Lmag))

        #Get perpendicular b vector. 
        self.b = self.smile_loc + (self.smag*np.cos(self.limb_c)*self.L)/self.Lmag
        self.bmag = np.sqrt(self.b[0]**2 + self.b[1]**2 + self.b[2]**2) 
        self.b_unit = self.b/self.bmag 

        #Get unit vectors for the image. 
        self.xi_unit = self.b_unit

        self.yi = np.cross(self.b, self.L) 
        self.yi_unit = self.yi/(self.yi[0]**2 + self.yi[1]**2 + self.yi[2]**2)**0.5 
        
        
                   
    def plot_image(self, header='original', cmap='lundi', vmin=0, vmax=None, add_overlays=True):
        '''This will plot a simple image of the chosen extension for you.
        
        Parameters
        ----------
        header - 'original' for original image. 
                'primary' for fitted image. 
                'intensity' for fitted intensity map. 
                
        '''
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()
            
        #Create the figure. 
        fig = plt.figure(figsize=(6,6))
        fig.subplots_adjust(top=0.8, left=0.20, right=0.85)
        ax = fig.add_subplot(111) 
        
        
        #Select correct array. 
        if header.lower() == 'original':
            array = self.original
            cbar_label = "Counts/pixel"
        elif header.lower() == 'primary':
            array = self.fitted
            cbar_label = "Counts/pixel"
        elif header.lower() == 'intensity':
            array = self.los_intensity
            cbar_label = "Intensity [keV cm"+r"$^{-2}$"+" s"+r"$^{-1}$"+" sr"+r"$^{-1}$"+"]"
        elif header.lower() == 'smoothed':
            array = self.smoothed_image 
            cbar_label = "Counts/pixel" 
        else:
            raise ValueError("Unsure what you wish to plot. Current options are 'original', 'primary', 'intensity' or 'smoothed'.")
        
        #Make axis one. 
        make_image_axes.make_image_axes(ax, array, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title=cbar_label)
            
        #Try adding coordinate axes at aim point. 
        if add_overlays:
            overlays.add_axes_at_aim_point(ax, self.xi_unit, self.yi_unit, -self.L_unit, self.smile_loc, self.target_loc)  
            
            #Try adding a magnetopause overlay. 
            #mpause_overlay.mpause(ax, self.xi_unit, self.yi_unit, -self.L_unit, self.smile_loc, self.params_best)
            
        #Sort out a title. 
        #coords = 'GSE' if self.use_GSE else 'GSM'
        pos = 'SMILE_{}: ({:.2f},{:.2f},{:.2f})'.format('GSE', *self.smile_loc) 
        aim = 'Aim: {:.2f}'.format(self.target_loc[0])  
        
        if (header.upper() == 'ORIGINAL') or (header.upper() == 'PRIMARY'):
            exposure_title = '\nExposure = {} s'.format(self.exposure)
        else:
            exposure_title = ''    
        #Title
        if self.time_string == '':
            fig.text(0.5, 0.95, self.filename[0:60]+'\n'+self.filename[60:]+'\n'+pos+' '+aim+exposure_title+'\n\n', fontsize=10)
        else:
            
            fig.text(0.5, 0.95, self.filename[0:60]+'\n'+self.filename[60:]+'\n'+pos+' '+aim+'\n{} {} {}\n\n'.format(*self.time_string.split('T'), exposure_title), fontsize=10, va='top', ha='center')
         
        if self.bayane: 
            fig.text(0.95, 0.05, f'Cusp = {self.cusp}', ha='right')    
            
    def plot_image_comparison(self, cmap='lundi', vmin=0, vmax=None, add_overlays=True, save=True, ylims=(-13.5,13.5), add_histograms=True):
        '''This will plot the fitted image and original image together with the fit parameters.
        
        
        Parameters
        ----------
        cmap - colourmap of image. def = 'lundi' 
        vmin - minimum colour value. def = 0 
        vmax - maximum colour value. 
        add_overlays - adds projected magnetopause and axes. 
        save - boolean to save the plot. 
        ylims - limits on y axes (needed for reduced images). 
        add_histograms - boolean to add the histograms to the top and right of the plot. 
                
        '''
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()
            
        #Create the figure. 
        fig = plt.figure(figsize=(8,6))
        fig.subplots_adjust(top=0.9, left=0.10, right=0.90, wspace=0.5, bottom=0.05)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122) 
        
       
        if self.counts_file:
            cbar_label = "Counts/pixel"
        else:
            cbar_label = "Intensity [keV cm"+r"$^{-2}$"+" s"+r"$^{-1}$"+" sr"+r"$^{-1}$"+"]"
        
        
        #Make axis one. 
        make_image_axes.make_image_axes(ax1, self.original, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title=cbar_label, add_histograms=add_histograms)
        
        make_image_axes.make_image_axes(ax2, self.fitted, self.xdeg_min, self.ydeg_min, self.n_pixels, self.m_pixels, cmap=cmap, vmin=vmin, vmax=vmax, cbar_title=cbar_label, add_histograms=add_histograms)
        

        ax1.set_title('Original\n\n', fontsize=10)
        if self.rcmem_err == 9999:
            ax2.set_title('Fitted\nr'+r'$_{CMEM}$'+r'$ = {:.2f}$ RE'.format(self.rcmem)+'\n', fontsize=10)
        else:
            ax2.set_title('Fitted\nr'+r'$_{CMEM}$'+r'$ = {:.2f} \pm {:.2f}$ RE'.format(self.rcmem, self.rcmem_err)+'\n', fontsize=10) 
        
        ax1.set_ylim(ylims)
        ax2.set_ylim(ylims) 
        
        #Try adding coordinate axes at aim point. 
        if add_overlays:
            #For axis one. 
            overlays.add_axes_at_aim_point(ax1, self.xi_unit, self.yi_unit, -self.L_unit, self.smile_loc, self.target_loc)  
            
            #Try adding a magnetopause overlay. 
            #mpause_overlay.mpause(ax1, self.xi_unit, self.yi_unit, -self.L_unit, self.smile_loc, self.params_best)
            
            #For axis two. 
            overlays.add_axes_at_aim_point(ax2, self.xi_unit, self.yi_unit, -self.L_unit, self.smile_loc, self.target_loc)  
            
            #Try adding a magnetopause overlay. 
            mpause_overlay.mpause(ax2, self.xi_unit, self.yi_unit, -self.L_unit, self.smile_loc, self.params_best)
        
        #Sort out a title. 
        #coords = 'GSE' if self.use_GSE else 'GSM'
        pos = 'SMILE_{}: ({:.2f},{:.2f},{:.2f})'.format('GSE', *self.smile_loc) 
        aim = 'Aim: {:.2f}'.format(self.target_loc[0])  
        
        if self.counts_file:
            exposure_title = '\nExposure = {} s'.format(self.exposure)
        else:
            exposure_title = ''  
         
        #Title
        if self.time_string == '':
            fig.text(0.5, 0.95, self.filename[0:60]+'\n'+self.filename[60:]+'\n'+pos+' '+aim+exposure_title+'\n\n', fontsize=10)
        else:
            
            fig.text(0.5, 0.95, self.filename[0:60]+'\n'+self.filename[60:]+'\n'+pos+' '+aim+'\n{} {} {}\n\n'.format(*self.time_string.split('T'), exposure_title), fontsize=10, va='top', ha='center')
        
        if self.bayane: 
            fig.text(0.95, 0.05, f'Cusp = {self.cusp}', ha='right')    
            
        #Save the plot. 
        if save: 
            figname = self.path+self.filename[0:-5]+'.png'
            fig.savefig(figname)     
    
    def plot_p_variation(self):
        '''This will plot the data in the parameter variation table.''' 
        
        #Make figure. 
        plt.close("all")
        fig = plt.figure(figsize=(6,8))
        fig.subplots_adjust(wspace=0.3, top=0.8)
        
        #For each parameter. 
        for i, param in enumerate(self.param_names):
            
            ax = fig.add_subplot(int(np.ceil(len(self.param_names)/2)),2,i+1) 
            
            pdata = self.p_variation[param] 
            x = np.arange(len(pdata))
            
            ax.plot(x, pdata) 
            ax.set_ylabel(param)
            ax.set_xlabel('Function Evaluations') 
            ax.grid()
            print (pdata)
        
        pos = 'SMILE_{}: ({:.2f},{:.2f},{:.2f})'.format('GSE', *self.smile_loc) 
        aim = 'Aim: {:.2f}'.format(self.target_loc[0])  
        
        if self.counts_file:
            exposure_title = '\nExposure = {} s'.format(self.exposure)
        else:
            exposure_title = ''  
         
             
        #Title
        if self.time_string == '':
            fig.text(0.5, 0.95, self.filename[0:60]+'\n'+self.filename[60:]+'\n'+pos+' '+aim+exposure_title+'\n\n', fontsize=10)
        else:
            
            fig.text(0.5, 0.95, self.filename[0:60]+'\n'+self.filename[60:]+'\n'+pos+' '+aim+'\n{} {} {}\n\n'.format(*self.time_string.split('T'), exposure_title), fontsize=10, va='top', ha='center')
            
            
class animate():
    '''This class contains all the information in a useable format for the fitted data FITS files. Adds it all to an object to be used by other analysis programs.''' 
    
    def __init__(self, filestring='SMILE_SXI_L3_SCIM60x60-SCI-CXF_*_V01_fitted_cxfov_cmem2g_mpfit_BATSRUS_hybrid1.fits', path='/data/smile/shared/sims/GAC_V2_Sam/GAC_V2.2_1x1_300s/'):
        '''This read in the FITS file and opens it.
        
        Parameters
        ----------
        filestring - Unix file string that will list all files with a certain naming convention.  
        path - Path to the FITS file. 
        
        
        '''
        
        self.path = path 
        self.filestring = filestring
        
        #Get all the filenames. 
        self.filenames = sorted(glob.glob1(self.path, self.filestring)) 
        
        #Read all those filenames. 
        #self.results = [results(filename=f, path=self.path, sw_info=False) for f in self.filenames] 
        
    
    def animate_results(self, cmap='lundi', vmin=0, vmax=150, fps=10, add_overlays=True):
        '''This will make an animation cycling through all the images. 
        CURRENTLY DOES NOT ADD THE HISTOGRAMS! 
        
        '''
        
        #Make initial image and extract first result. 
        result = results(filename=self.filenames[0], path=self.path, sw_info=False) 
        self.result0 = result 
        
        #Get custom lundi colormap.
        if cmap == 'lundi':
            cmap = read_cmap.txt2matplotlib()
            
        #Create the figure.
        plt.close("all") 
        fig = plt.figure(figsize=(8,6))
        fig.subplots_adjust(top=0.8, left=0.10, right=0.90, wspace=0.5, bottom=0.1)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122) 
        
        
        #Get 1D pixel arrays for plotting. The edges of the pixels.  
        xarray_los = np.linspace(result.xdeg_min, -result.xdeg_min, result.m_pixels+1)
        yarray_los = np.linspace(result.ydeg_min, -result.ydeg_min, result.n_pixels+1)
        
        #Make 2D arrays for x and y. 
        X_los, Y_los = np.meshgrid(xarray_los, yarray_los)
        
        #Make the images 
        mesh_data = ax1.pcolormesh(X_los, Y_los, result.original, cmap=cmap, vmin=vmin, vmax=vmax)
        mesh_model = ax2.pcolormesh(X_los, Y_los, result.fitted, cmap=cmap, vmin=vmin, vmax=vmax)
        
        #Add colourbars. 
        cbar1 = plt.colorbar(mesh_data, ax=ax1, shrink=0.5, pad=0.1)
        cbar1.set_label("Counts per Pixel")
        cbar2 = plt.colorbar(mesh_model, ax=ax2, shrink=0.5, pad=0.1)
        cbar2.set_label("Counts per Pixel")
        
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        
        title_data = ax1.set_title('Original', fontsize=10)
        title_model = ax2.set_title('Fitted\nr'+r'$_{CMEM}$'+r'$ = {:.2f} \pm {:.2f}$ RE'.format(result.rcmem, result.rcmem_err), fontsize=10) 
        
        pos = 'SMILE_{}: ({:.2f},{:.2f},{:.2f})'.format('GSE', *result.smile_loc) 
        aim = 'Aim: {:.2f}'.format(result.target_loc[0])  
        exposure_title = '\nExposure = {} s'.format(result.exposure)
        
        title_fig = fig.text(0.5, 0.95, self.filenames[0][0:62]+'\n'+self.filenames[0][62:]+'\n'+pos+' '+aim+'\n{} {} {}\n\n'.format(*result.time_string.split('T'), exposure_title), fontsize=10, va='top', ha='center')
        
        if add_overlays:
            
            #Work out projections of axes. 
            alphax, betax, alphay, betay, alphaz, betaz = overlays.add_axes_at_aim_point(ax1, result.xi_unit, result.yi_unit, -result.L_unit, result.smile_loc, result.target_loc, add_to_plot=False, add_labels=False)
            xdata, = ax1.plot(alphax, betax, color='b', lw=2)
            ydata, = ax1.plot(alphay, betay, color='r', lw=2)
            zdata, = ax1.plot(alphaz, betaz, color='lime', lw=2)
            xmodel, = ax2.plot(alphax, betax, color='b', lw=2)
            ymodel, = ax2.plot(alphay, betay, color='r', lw=2)
            zmodel, = ax2.plot(alphaz, betaz, color='lime', lw=2)
            
            #Add labels to clarify what each axis line is. 
            ax1.text(1.11, 0.99, 'x', color='b', ha='left', va='top', transform=ax1.transAxes)
            ax1.text(1.11, 0.94, 'y', color='r', ha='left', va='top', transform=ax1.transAxes)
            ax1.text(1.11, 0.89, 'z', color='lime', ha='left', va='top', transform=ax1.transAxes)
            ax2.text(1.11, 0.99, 'x', color='b', ha='left', va='top', transform=ax2.transAxes)
            ax2.text(1.11, 0.94, 'y', color='r', ha='left', va='top', transform=ax2.transAxes)
            ax2.text(1.11, 0.89, 'z', color='lime', ha='left', va='top', transform=ax2.transAxes)
             
        def update(frame):
            print (frame)
            #Get the next data file. 
            result = results(filename=self.filenames[frame], path=self.path, sw_info=False) 
            
            #Update the meshes. 
            mesh_data.set_array(result.original)
            mesh_model.set_array(result.fitted) 
            
            #Update model title. 
            title_model.set_text('Fitted\nr'+r'$_{CMEM}$'+r'$ = {:.2f} \pm {:.2f}$ RE'.format(result.rcmem, result.rcmem_err))
            
            #Update figure title. 
            pos = 'SMILE_{}: ({:.2f},{:.2f},{:.2f})'.format('GSE', *result.smile_loc) 
            aim = 'Aim: {:.2f}'.format(result.target_loc[0])  
            exposure_title = '\nExposure = {} s'.format(result.exposure)
            
            title_fig.set_text(self.filenames[frame][0:62]+'\n'+self.filenames[frame][62:]+'\n'+pos+' '+aim+'\n{} {} {}\n\n'.format(*result.time_string.split('T'), exposure_title))
            
            #Work out projections of axes. 
            if add_overlays:
                alphax, betax, alphay, betay, alphaz, betaz = overlays.add_axes_at_aim_point(ax1, result.xi_unit, result.yi_unit, -result.L_unit, result.smile_loc, result.target_loc, add_to_plot=False, add_labels=False)
                
                xdata.set_xdata(alphax)
                xdata.set_ydata(betax)
                ydata.set_xdata(alphay)
                ydata.set_ydata(betay)
                zdata.set_xdata(alphaz)
                zdata.set_xdata(betaz) 
                
                xmodel.set_xdata(alphax)
                xmodel.set_ydata(betax)
                ymodel.set_xdata(alphay)
                ymodel.set_ydata(betay)
                zmodel.set_xdata(alphaz)
                zmodel.set_xdata(betaz) 
                
                return  (mesh_data, mesh_model, title_model, title_fig, xdata, ydata, zdata, xmodel, ymodel, zmodel)    
            else:
                return  (mesh_data, mesh_model, title_model, title_fig)
            
        
        #Run the animation. 
        animname = f'{self.filestring[0:-5]}.mp4'
        ani = animation.FuncAnimation(fig=fig, func=update,  frames=len(self.filenames), interval=int(1000/fps), repeat=False) 
        FFwriter = animation.FFMpegWriter(fps=fps)
        ani.save(self.path+animname, writer=FFwriter) 
        print ('Saved: ', self.path+animname)   
