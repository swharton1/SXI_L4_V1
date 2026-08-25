#This calculates a set of lines on the magnetopause that can then be overlaid on the image. 

import numpy as np 
from . import boundary_emissivity_functions as bef
from . import overlays 
from . import coord_conv as cconv 

class mpause():
    '''This class will work out a magnetopause surface as a series of lines that can be overplotted.''' 
    
    def __init__(self, ax, xi, yi, zi, sxi_loc, params, model='cmem2g', dipole=0, pdyn=None, pmag=None, bz=None, color='b', lw=1, marker=None):
        '''This takes in the parameters of the magnetopause. 
        
        Parameters
        ----------
        ax - axis object 
        xi - image unit vector x
        yi - image unit vector y 
        zi - image unit vector z 
        sxi_loc - location of smile 
        params - model parameters 
        model - 'cmem2g' is default
        dipole
        pdyn
        pmag
        bz
        
        '''
        

        
        
        #Get coefficients and extract the right parameters. 
        if model in ['jorg']:
            
            #Get p0, p1, p2 and p3
            r0 = params[0]
            ay = params[7]
            az = params[8]
            
        
        elif model in ["cmem", "acmem", "cmem2a", "cmem2b"]:
            #Get Lin coefficients. 
            self.lin_coeffs = bef.get_lin_coeffs(dipole, pdyn, pmag, bz)
            
            p0 = params[0]
            p1 = params[7]
            p2 = params[8]
            p3 = params[9]
            
        elif model in ["cmem2c", "cmem2d", "cmem2e", "cmem2f", "cmem2g"]:
            #Get Lin coefficients. 
            self.lin_coeffs = bef.get_lin_coeffs_simple(dipole) 
            
            #Get p0, p1, p2 and p3
            if model == 'cmem2f':
                #Get subsolar magnetopause. 
                p0 = params[0]
                p1 = params[5]
                p2 = params[6]
                p3 = params[7]
            elif model == 'cmem2g':
                p0 = params[0]
                p1 = params[4]
                p2 = params[5]
                p3 = params[6] 
            else:
                #Get subsolar magnetopause. 
                p0 = params[0]
                p1 = params[7]
                p2 = params[8]
                p3 = params[9]
        else:
            raise ValueError("Not picked a cmem model")
            
        #Create the arrays of lines you will need. 
        theta = np.deg2rad(np.linspace(10,90,9))
        phi = np.deg2rad(np.linspace(0,360,361))
        
        #For each line, calculate r, then convert to cartesian, then add to the image. 
        for t in theta:
            t = np.ones(phi.size)*t 
            
            #Calculate r. 
            if model in ['jorg']:
                r = bef.shue_func(t, phi, r0, ay, az)
            
            elif model in ["cmem", "acmem", "cmem2a", "cmem2b"]:
                r = bef.lin_scaled_func(t, phi, *self.lin_coeffs, p0=p0, p1=p1, p2=p2, p3=p3)
            elif model in ["cmem2c", "cmem2d", "cmem2e", "cmem2f", "cmem2g"]:  
                r = bef.lin_scaled_func_simple(t, phi, *self.lin_coeffs, p0=p0, p1=p1, p2=p2, p3=p3)

            else:
                raise ValueError("Not picked a cmem model")    
            
            
            #Now convert to Cartesian. 
            x, y, z = cconv.convert_shue_to_xyz_coords(r, t, phi)   
            
            #Now add to the image. 
            overlays.transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color=color, lw=lw, marker=marker)   
            
        #Now do the other direction. 
        theta = np.deg2rad(np.linspace(0,90,91))
        phi = np.deg2rad(np.linspace(0, 330, 12))   
        
        #For each line, calculate r, then convert to cartesian, then add to the image. 
        for p in phi:
            p = np.ones(theta.size)*p 
            
            #Calculate r. 
            if model in ["cmem", "acmem", "cmem2a", "cmem2b"]:
                r = bef.lin_scaled_func(theta, p, *self.lin_coeffs, p0=p0, p1=p1, p2=p2, p3=p3)
            elif model in ["cmem2c", "cmem2d", "cmem2e", "cmem2f", "cmem2g"]:  
                r = bef.lin_scaled_func_simple(theta, p, *self.lin_coeffs, p0=p0, p1=p1, p2=p2, p3=p3)

            else:
                raise ValueError("Not picked a cmem model")    
            
            
            #Now convert to Cartesian. 
            x, y, z = cconv.convert_shue_to_xyz_coords(r, theta, p)        
            
            #Now add to the image. 
            overlays.transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color=color, lw=lw, marker=marker) 
