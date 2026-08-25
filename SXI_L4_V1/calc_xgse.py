#This will do a calculation of the range of the X GSE axis in the image. 
#This is a simplication of the original as it does not use a file. Just give it the position. 

import numpy as np 
import os 
#from . import read_fits_image

class get_image():
    '''This will take in a fits file image, and use the header information to work out the X-GSE range. Or you can give it the smile_loc directly. ''' 
    
    def __init__(self, smile_loc=None):
        '''Takes in the spacecraft position and uses the guidance law to work out the bounds. ''' 
        
        #if smile_loc is None: 
        
        #    #Open file. 
        #    self.rf = read_fits_image.read_fits(filename=filename, path=path, header_type='PRIMARY', plot_image=False)
        
            #Extract SMILE position and aim. 
        #    self.smile_loc = self.rf.smile_loc
            #self.target_loc = self.rf.target_loc 
        #else:
        self.smile_loc = smile_loc
            
        #Get magnitude of spacecraft vector - its radial position. 
        self.smag = np.sqrt(self.smile_loc[0]**2 + self.smile_loc[1]**2 + self.smile_loc[2]**2)
        
        #Get alpha angle. 
        self.get_alpha_angle()
        
        #Get earth angle. 
        self.r_angle = np.arcsin(1/self.smag) 
        
        #Get combined limb and r angle. 
        #Added correction - 20.3 deg to centre of FOV, not edge. 
        self.limb = 20.3*(np.pi/180)
        self.limb_c = (self.limb + self.r_angle)
        
        #Set up array of angles in the FOV. 
        self.dtheta = np.linspace(-8,8,17)*(np.pi/180) 
        
        #Get the X-GSE distance to all these points. 
        self.x_gse = self.get_xgse(self.dtheta) 
        
    def get_alpha_angle(self, dtheta=0):
        '''This will calculate alpha, the angle between the spacecraft vector and the perpendicular to the x axis. '''
        self.alpha_angle = np.arctan2(self.smile_loc[0],np.sqrt(self.smile_loc[1]**2 + self.smile_loc[2]**2))   
    
    def get_xgse(self, dtheta=0):
        '''This will calculate the x-gse distance for a given look angle. lx + sx ''' 
        
        return np.sqrt(self.smile_loc[1]**2 + self.smile_loc[2]**2)*np.tan(self.limb_c-self.alpha_angle+dtheta) + self.smile_loc[0]
       

