#THIS WILL CONTAINS FUNCTIONS TO CONVERT BETWEEN 3D COORDINATE SYSTEMS. 
import numpy as np
#Cartesian-SHUE coordinates. 

def convert_xyz_to_shue_coords(x, y, z):
    '''This will convert the x,y,z coordinates to those used in the Shue model 
     of the magnetopause and bowshock. 

    Parameters
    ---------
    x, y, z - now 3D. Must be entered as numpy arrays.  

    Returns
    -------
    r, theta (rad) and phi (rad)
    '''

   
    # r 
    r = (x**2 + y**2 + z**2)**0.5
       
    # theta - only calc. where coordinate singularities won't occur. 
    theta = np.zeros(r.shape)
    i = np.where(r != 0)
    theta[i] =  np.arccos(x[i]/r[i])

    # phi - only calc. where coordinate singularities won't occur. 
    phi = np.zeros(r.shape)
    phi[i] = np.arctan2(z[i], y[i]) 
     
    return r, theta, phi
    
def convert_shue_to_xyz_coords(r, theta, phi):
    '''This will convert the Shue coordinates back to xyz coordinates. 
        
    Parameters
    ----------
    r, theta (rad), phi (rad)
        
    Returns
    -------
    x,y,z
    '''

    x = r*np.cos(theta)
    y = r*np.sin(theta)*np.cos(phi)
    z = r*np.sin(theta)*np.sin(phi)

    return x,y,z 

#Cartesian-Spherical coords. 
def convert_xyz_to_spherical_coords(x, y, z):
    '''This will convert Cartesian coordinates to spherical coordinates.
    
    Parameters
    ----------
    x, y, z 
    
    Returns
    -------
    r, alpha (rad), beta (rad) 
    '''
    
    # r 
    r = (x**2 + y**2 + z**2)**0.5
    
    # alpha - only calc. where coordinate singularities won't occur. 
    alpha = np.zeros(r.shape)
    i = np.where(r != 0)
    alpha[i] =  np.arccos(z[i]/r[i])

    # beta - only calc. where coordinate singularities won't occur. 
    beta = np.zeros(r.shape)
    beta[i] = np.arctan2(y[i], x[i]) 
     
    return r, alpha, beta 

def convert_spherical_to_xyz_coords(r, alpha, beta):
    '''This will convert spherical to Cartesian coordinates.
    
    Parameters
    ----------
    r, alpha (rad), beta (rad) 
    
    Returns
    -------
    x, y, z
    '''
    
    x = r*np.sin(alpha)*np.cos(beta)
    y = r*np.sin(alpha)*np.sin(beta)
    z = r*np.cos(alpha) 
    
    return x, y, z 
    
#Spherical-Shue coords. 
def convert_shue_to_spherical_coords(r, theta, phi):
    '''This will convert the Shue coordinates to spherical coordinaes. 
    
    Parameters
    ----------
    r, theta (rad), phi (rad)
    
    Returns
    -------
    r, alpha (rad), beta (rad)
    '''
    
    alpha = np.zeros(r.shape)
    i = np.where(r != 0) 
    alpha[i] = np.arccos(np.sin(phi[i])*np.sin(theta[i]))
    
    beta = np.zeros(r.shape)
    beta[i] = np.arctan2(np.cos(phi[i])*np.sin(theta[i]), np.cos(theta[i])) 
    
    
    return r, alpha, beta 
    
def convert_spherical_to_shue_coords(r, alpha, beta):
    '''This will convert the spherical coordinates to Shue coordinates.
    
    Parameters
    ----------
    r, alpha (rad), beta (rad)
    
    Returns
    -------
    r, theta (rad), phi (rad) 
    '''
    
    theta = np.zeros(r.shape)
    i = np.where(r != 0)
    theta[i] = np.arccos(np.cos(beta[i])*np.sin(alpha[i]))
    
    phi = np.zeros(r.shape)
    phi[i] = np.arctan2(np.cos(alpha[i]), np.sin(alpha[i])*np.sin(beta[i]))
    return r, theta, phi 
    





#Cartesian-Aberrated coordinates. 

def convert_xyz_to_aberrated_shue_coords(x, y, z, dgamma, ddelta): 
    '''This will convert the x,y,z coordinates to those used in the Shue model 
     of the magnetopause and bowshock. 

    Parameters
    ---------
    x, y, z - now 3D. Must be entered as numpy arrays.  
    dgamma - rotation in the horizontal plane (rad)
    ddelta - rotation in the vertical plane (rad) 
    
    Returns
    -------
    r, theta (rad) and phi (rad)
    '''

    #Horizontal rotation. 
    x2 = x*np.cos(-dgamma) - y*np.sin(-dgamma)
    y2 = x*np.sin(-dgamma) + y*np.cos(-dgamma) 
    z2 = z 
    
    #Vertical rotation. 
    x3 = x2*np.cos(-ddelta) - z2*np.sin(-ddelta)
    y3 = y2 
    z3 = x2*np.sin(-ddelta) + z2*np.cos(-ddelta) 


    # r 
    r = (x3**2 + y3**2 + z3**2)**0.5
       
    # theta - only calc. where coordinate singularities won't occur. 
    theta = np.zeros(r.shape)
    i = np.where(r != 0)
    theta[i] =  np.arccos(x3[i]/r[i])

    # phi - only calc. where coordinate singularities won't occur. 
    phi = np.zeros(r.shape)
    phi[i] = np.arctan2(z3[i], y3[i]) 
     
    return r, theta, phi


