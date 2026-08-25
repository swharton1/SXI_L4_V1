#This is based on the overlays code from JXI. This will allow me to put axes, boundaries and other objects into the FOV. 

import numpy as np 
from . import transformations as trans 
from . import coastlines 

def add_planetary_outline(ax, xi, yi, zi, sxi_loc, limb, lw=0.5, colour='o', marker=None, alpha_colour=0.5, fill=False):
    '''This adds a circular outline for a planet or moon. This is now using an updated method over the JXI code. It actually calculates the circle in image coordinates and then converts to angles instead of just drawing a circle in the angular coordinates. 
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    limb - limb angle to planet/moon in radians. 
    colour - colour of outline
    lw - linewidth of outline 
    alpha_colour - transparency if filling it in. 
    fill - False
    
    Return 
    ------
    alpha - alpha coordinate
    beta - beta coordinate
    
    '''
    
    print ('Add planetary outline...')
    
    x0 = np.array([0])
    y0 = np.array([0])
    z0 = np.array([0])
        
    #Package into a single vector.    
    positions = np.zeros((3,len(x0)))
    positions[0] = x0
    positions[1] = y0
    positions[2] = z0

    #Define toy FOV object. 
    fov = fov_class(xi, yi, zi) 
    
    
    #Convert to image coordinates. This is the centre of the circle in image coords. 
    pos_i = trans.convert_cartesian_to_image_coords_2(fov, positions.T, sxi_loc)
    
    
    #Do a test here to test the reverse transformation. This test works.  
    #pos_c = trans.convert_image_to_cartesian_coords_2(fov, pos_i.T, sxi_loc) 
    
    
    #Get the centre of the circle in image coordinates. 
    #This is also the normal vector as the spacecraft is at the origin. 
    c = np.squeeze(pos_i) 

    #Get radius in image coordinates (same as GSE). 
    r=1
    
    #Now get two vectors orthogonal to the normal vector. 
    u = np.array([-c[1], c[0], 0]) 
    v = np.cross(c, u) 
    
    #Get them as unit vectors. 
    uhat = u/np.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    vhat = v/np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    
    #Now create the locus of points defining the circle. 
    t_angles = np.linspace(0,2*np.pi,361)
    locus = np.zeros((3,361))

    for t, tval in enumerate(t_angles):
        locus[0,t] = c[0] + r*np.cos(tval)*uhat[0] + r*np.sin(tval)*vhat[0]
        locus[1,t] = c[1] + r*np.cos(tval)*uhat[1] + r*np.sin(tval)*vhat[1]
        locus[2,t] = c[2] + r*np.cos(tval)*uhat[2] + r*np.sin(tval)*vhat[2]
 
    
    #Convert to alpha/beta coordinates. 
    alpha_locus, beta_locus = trans.get_angular_position(locus)
    
    alpha_locus = np.rad2deg(alpha_locus)
    beta_locus = np.rad2deg(beta_locus) 
    
    #Add circle to image.    
    ax.fill(alpha_locus, beta_locus, fill=fill, color=colour, lw=lw, alpha=alpha_colour, ec=None)


def add_terminator(ax, xi, yi, zi, sxi_loc, lw=2, marker=None):
    '''This will add on the terminator.  
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    time - dt.datetime object of the time of day in UT. 
    
    '''
    
    print ("Terminator...") 
    
    #Define in GSE. 
    t = np.linspace(0, 2*np.pi, 361) 
    y0 = np.cos(t)
    z0 = np.sin(t)
    x0 = np.zeros(len(t)) 

    #New one for SXI. 
    transform_to_image(x0, y0, z0, ax, xi, yi, zi, sxi_loc, color='r', lw=lw, marker=marker, apply_shadow=True)
    
   
        
        
def add_nightside(ax, xi, yi, zi, sxi_loc, lw=2, marker=None):
    '''This will shade in between the terminator and the planetary outline on the nightside.   
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    time - dt.datetime object of the time of day in UT. 
    
    '''
    
    print ("Nightside...")   
    
    #Define terminator in GSE first. 
    t = np.linspace(0, 2*np.pi, 361) 
    y0 = np.cos(t)
    z0 = np.sin(t)
    x0 = np.zeros(len(t)) 
    
    #Get alpha and beta values of the visible terminator.  
    #alpha_term, beta_term = transform_to_image(x0, y0, z0, ax, xi, yi, zi, sxi_loc, color='r', lw=lw, marker=marker, apply_shadow=False, add_to_plot=False, return_angles=True)
    
    #Package into a single vector.    
    positions = np.zeros((3,len(x0)))
    positions[0] = x0
    positions[1] = y0
    positions[2] = z0

    #Define toy FOV object. 
    fov = fov_class(xi, yi, zi) 
    
    
    #Convert to image coordinates. 
    pos_i = trans.convert_cartesian_to_image_coords_2(fov, positions.T, sxi_loc)
           
    #Convert to alpha/beta coordinates. 
    alpha_term, beta_term = trans.get_angular_position(pos_i)
    
    alpha_term = np.rad2deg(alpha_term)
    beta_term = np.rad2deg(beta_term) 
    
    #Sort out the shadowing and reordering here. 
    #For filter, you need the radial distance to each point on the line. 
    mer_r = (pos_i[0]**2 + pos_i[1]**2 + pos_i[2]**2)**0.5 
        
    #Calculate range to Jupiter again, as you need it to decide which meridians to show. 
    sxi_r = (sxi_loc[0]**2 + sxi_loc[1]**2 + sxi_loc[2]**2)**0.5 
         
    #Filter by those on this side of the limb. 
    #near = pos_i[0] < sxi_r
    near = mer_r < (sxi_r**2 - 1)**0.5 
    
    #Check if array needs reordering. 
    idarray = make_idarray(near)
    near, alpha_term, beta_term = reorder_array_by_id(idarray, near, alpha_term, beta_term) 
    

    alpha_term = alpha_term[near]
    beta_term = beta_term[near]
    
    #Define the planetary outline in image coordinates next. 
    x0 = np.array([0])
    y0 = np.array([0])
    z0 = np.array([0])
        
    #Package into a single vector.    
    positions = np.zeros((3,len(x0)))
    positions[0] = x0
    positions[1] = y0
    positions[2] = z0

    #Define toy FOV object. 
    fov = fov_class(xi, yi, zi) 
    
    
    #Convert to image coordinates. This is the centre of the circle in image coords. 
    pos_i = trans.convert_cartesian_to_image_coords_2(fov, positions.T, sxi_loc)
    
    #Get the centre of the circle in image coordinates. 
    #This is also the normal vector as the spacecraft is at the origin. 
    c = np.squeeze(pos_i) 

    #Get radius in image coordinates (same as GSE). 
    r=1
    
    #Now get two vectors orthogonal to the normal vector. 
    u = np.array([-c[1], c[0], 0]) 
    v = np.cross(c, u) 
    
    #Get them as unit vectors. 
    uhat = u/np.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    vhat = v/np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    
    #Now create the locus of points defining the circle. 
    t_angles = np.linspace(0,2*np.pi,361)
    locus = np.zeros((3,361))

    for t, tval in enumerate(t_angles):
        locus[0,t] = c[0] + r*np.cos(tval)*uhat[0] + r*np.sin(tval)*vhat[0]
        locus[1,t] = c[1] + r*np.cos(tval)*uhat[1] + r*np.sin(tval)*vhat[1]
        locus[2,t] = c[2] + r*np.cos(tval)*uhat[2] + r*np.sin(tval)*vhat[2]
    
    #Convert to alpha/beta coordinates. 
    alpha_locus, beta_locus = trans.get_angular_position(locus)
    
    #These are the alpha/beta coordinates of the full locus. 
    alpha_locus = np.rad2deg(alpha_locus)
    beta_locus = np.rad2deg(beta_locus) 

        
    #Convert the locus of points forming the outline to GSE so you can filter them.   
    locus_c = trans.convert_image_to_cartesian_coords_2(fov, locus.T, sxi_loc) 
       

    #Filter by x value. 
    night = locus_c[0] <= 0 
    
    #Check if the array needs reordering. 
    idarray = make_idarray(night)
    night, alpha_locus, beta_locus = reorder_array_by_id(idarray, night, alpha_locus, beta_locus)
     
    #Filter the array. 
    alpha_night = alpha_locus[night]
    beta_night = beta_locus[night] 
    
    #Concatenate the terminator and the nightside outline. 
    alpha_comb = np.concatenate([alpha_night, alpha_term, [alpha_night[0]]])
    beta_comb = np.concatenate([beta_night, beta_term, [beta_night[0]]]) 
    
    #Use this filter on alpha and beta. 
    #ax.plot(alpha_night, beta_night, color='g', lw=2, marker=marker)  
    #ax.plot(alpha_term, beta_term, color='r', lw=2, marker=marker) 
    #ax.plot(alpha_comb, beta_comb, color='r', lw=2, marker=marker)   
    ax.fill(alpha_comb, beta_comb, fill=True, color='grey', lw=lw, alpha=0.5, ec=None)
    
    
    
def add_coastlines(ax, xi, yi, zi, sxi_loc, time, lw=0.5, marker=None):
    '''This will add on the Earth's coastlines. 
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    time - dt.datetime object of the time of day in UT. 
    
    '''
    
    print ("Earth's coastlines...") 
    
    #Read in the coastlines data and extract the x, y and z data. 
    coast = coastlines.coastlines(time=time) 

    
    unique = np.unique(coast.indices)
    
    #Loop through each island. 
    for u in unique: 
        i = coast.indices == u 
        
        x0 = coast.x[i] 
        y0 = coast.y[i]
        z0 = coast.z[i] 
    
        #New one for SXI. 
        transform_to_image(x0, y0, z0, ax, xi, yi, zi, sxi_loc, color='k', lw=lw, marker=marker, apply_shadow=True)
        
def add_parallels(ax, xi, yi, zi, sxi_loc, time, color='k', lw=0.5, marker=None):
    '''This will add the parallels on a planet.
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    time - dt.datetime object of the time of day in UT. 
    
    '''
    
    print ('Parallels...') 
    
    #Read in the coastlines data and extract the x, y and z data. 
    coast = coastlines.coastlines(time=time) 

    #Do all five parallels. 
    for i in range(5):
        
        x0 = coast.para_x[361*i:361*(i+1)]
        y0 = coast.para_y[361*i:361*(i+1)]
        z0 = coast.para_z[361*i:361*(i+1)]
    
        #New one for SXI. 
        transform_to_image(x0, y0, z0, ax, xi, yi, zi, sxi_loc, color=color, lw=lw, marker=marker, apply_shadow=True)
        
def add_meridians(ax, xi, yi, zi, sxi_loc, time, color='k', lw=0.5, marker=None):
    '''This will add the meridians on a planet.
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    time - dt.datetime object of the time of day in UT. 
    
    '''
    
    print ('Meridians...') 
    
    #Read in the coastlines data and extract the x, y and z data. 
    coast = coastlines.coastlines(time=time) 

    #Do all five parallels. 
    for i in range(12):
        
        x0 = coast.meri_x[181*i:181*(i+1)]
        y0 = coast.meri_y[181*i:181*(i+1)]
        z0 = coast.meri_z[181*i:181*(i+1)]
    
        #New one for SXI. 
        transform_to_image(x0, y0, z0, ax, xi, yi, zi, sxi_loc, color=color, lw=lw, marker=marker, apply_shadow=True)
  
    
       
            
def add_ecliptic_plane(ax, xi, yi, zi, sxi_loc, color='b', lw=0.5, marker=None):
    '''This will add a projection of the ecliptic plane as a 1 RE grid.
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    color - colour of overlay
    lw - linewidth of overlay 
    
    '''
    
    print ('Add ecliptic plane') 
    
    #Add lines parallel to x first. 
    x = np.linspace(-10,20,31)
    z = np.zeros(x.size)
    
    for y in np.linspace(-20,20,41):
        y = np.ones(x.size)*y 
        transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color=color, lw=lw, marker=marker)
    
    #Add lines parallel to y next. 
    y = np.linspace(-20,20,41)
    z = np.zeros(y.size)
    
    for x in np.linspace(-10,20,31):
        x = np.ones(y.size)*x 
        transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color=color, lw=lw, marker=marker)    
    
def add_axes_at_aim_point(ax, xi, yi, zi, sxi_loc, target_loc, color='b', lw=2, marker=None, add_to_plot=True, add_labels=True):
    '''This will add y and z directions on top of the GSE x axis at the aim point. 
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    target_loc - (x,y,z) vector of the aim point. 
    color - colour of overlay
    lw - linewidth of overlay 
    add_to_plot - boolean to add lines to plot. def = True. 
    add_labels - boolean to add labels to plot. def = True. 
    
    Returns
    -------
    alphax - alpha positions for x axis. 
    betax - beta positions for x axis. 
    alphay - alpha positions for y axis. 
    betay - beta positions for y axis. 
    alphaz - alpha positions for z axis. 
    betaz - beta positions for z axis. 
    
    '''
    
    print ('Add axes at tangent point...') 
    
    #x axis first. 
    x = np.linspace(target_loc[0], target_loc[0]+1, 10)
    y = np.zeros(x.size)
    z = np.zeros(x.size) 
    alphax, betax = transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color='b', lw=lw, marker=marker, return_angles=True, add_to_plot=add_to_plot)
    
    #Add a small y axis. 
    x = np.ones(10)*target_loc[0]
    y = np.linspace(0,1,10)
    z = np.zeros(10) 
    alphay, betay = transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color='c', lw=lw, marker=marker, return_angles=True, add_to_plot=add_to_plot)
    
    #Add a small y axis. 
    x = np.ones(10)*target_loc[0]
    y = np.zeros(10)
    z = np.linspace(0,1,10) 
    alphaz, betaz = transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color='k', lw=lw, marker=marker, return_angles=True, add_to_plot=add_to_plot)
    
    #Add little labels so people know what each coloured line is. 
    if add_labels:
        ax.text(1.11, 0.99, 'x', color='b', ha='left', va='top', transform=ax.transAxes)
        ax.text(1.11, 0.94, 'y', color='c', ha='left', va='top', transform=ax.transAxes)
        ax.text(1.11, 0.89, 'z', color='k', ha='left', va='top', transform=ax.transAxes)
    
    return alphax, betax, alphay, betay, alphaz, betaz 
    
def add_x_axis(ax, xi, yi, zi, sxi_loc, color='b', lw=2, marker='+'):
    '''This will add the GSE x axis. 
    
    Parameters
    ----------
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    color - colour of overlay
    lw - linewidth of overlay 
    
    Returns 
    -------
    x, y, z - Cartesian positions of X-axis. 
    
    '''
    
    print ('GSE x-axis...') 
    
    #Define the x axis in GSE space. 
    x = np.arange(20)
    y = np.zeros(x.size)
    z = np.zeros(x.size) 
    
    transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color=color, lw=lw, marker=marker)
    
    return x, y, z 
    
def transform_to_image(x, y, z, ax, xi, yi, zi, sxi_loc, color='b', lw=2, marker=None, apply_shadow=False, add_to_plot=True, return_angles=False):
    '''This bit actually does the transformation of the structure to the image. 
    
    Parameters
    ----------
    x - x values of structure as an array. 
    y - y values of structure as an array. 
    z - z values of structure as an array. 
    
    ax - axis for the image. 
    xi - image unit vector x
    yi - image unit vector y 
    zi - image unit vector z 
    sxi_loc - (x,y,z) vector of SXI.
    color - colour of overlay
    lw - linewidth of overlay
    apply_shadow - boolean to apply shadow so you can't see what is behind the planet. 
    add_to_plot - boolean to actually add the final structure to the image. def = True. 
    return_angles - boolean to return the alpha and beta values in the plot. def = False. 
    '''
    
    #Package into a single vector.    
    positions = np.zeros((3,len(x)))
    positions[0] = x
    positions[1] = y
    positions[2] = z

    #Define toy FOV object. 
    fov = fov_class(xi, yi, zi) 
    
    
    #Convert to image coordinates. 
    pos_i = trans.convert_cartesian_to_image_coords_2(fov, positions.T, sxi_loc)
           
    #Convert to alpha/beta coordinates. 
    alpha, beta = trans.get_angular_position(pos_i)
    
    alpha = np.rad2deg(alpha)
    beta = np.rad2deg(beta) 

    if apply_shadow: 
        #For filter, you need the radial distance to each point on the line. 
        mer_r = (pos_i[0]**2 + pos_i[1]**2 + pos_i[2]**2)**0.5 
        
        #Calculate range to Jupiter again, as you need it to decide which meridians to show. 
        sxi_r = (sxi_loc[0]**2 + sxi_loc[1]**2 + sxi_loc[2]**2)**0.5 
        
        #Plot. 
        #Filter by those on this side of the limb. 
        #near = pos_i[0] < sxi_r
        near = mer_r < (sxi_r**2 - 1)**0.5 
        
        
        #You will need to plot in bits. 
        #Get ids for each section. zero is the far side so ignore. 
        idarray = make_idarray(near)
        id_unique = np.unique(idarray) 
        
        #Get maximum in idarray. 
        idmax = int(idarray.max())
        
        
        if (idmax > 0) & (add_to_plot):
            for idx in range(1,idmax+1): 
                i = idarray == idx 
                ax.plot(alpha[i], beta[i], color=color, lw=lw, marker=marker)  
        
        
        #ax.plot(alpha[near], beta[near], color=color, lw=lw, marker=marker) 
        if return_angles:
            return alpha[near], beta[near] 
    else:
        if add_to_plot:
            ax.plot(alpha, beta, color=color, lw=lw, marker=marker)
        if return_angles:
            return alpha, beta


def make_idarray(bool_array):
    '''This will identify and label all the sections in the array which are true, and give each section a number in an array called idarray.
    
    Parameters
    ----------
    bool_array - Array of booleans. 
    
    Returns
    -------
    idarray - Array with numerical labels identifying sections of True values. False values get zero. E.g. [1,1,1,0,0,2,2,2,0,0,3,0,4,4,4,4]
    
    '''
    
    #You will need to plot in bits. 
    idarray = np.zeros(len(bool_array))
    cnt = 1 

    for n, nval in enumerate(bool_array):
            
        if nval == True: 
            idarray[n] = cnt 
        elif (nval == False) & (bool_array[n-1] == True):
            cnt += 1 
        else:
            pass 

    return idarray 

def reorder_array_by_id(id_array, bool_array, alpha, beta):
    '''This will reorder the array if the True values in bool_array are split.
    
    Parameters
    ----------
    id_array - array of ids labeling the blocks of True values. Made by make_idarray()
    bool_array - Array of booleans. 
    alpha - Array of alpha values. 
    beta - Array of beta values. 
    
    Returns
    -------
    bool_array - Array of booleans but reordered.
    alpha - Array of alpha values but reordered.  
    beta - Array of beta values but reordered. 
    
    ''' 

    #Get maximum in idarray. 
    idmax = int(id_array.max())
    
    #If there are two blocks (assumed), stick the front on the end. 
    #Otherwise, do nothing.  
    if idmax > 1:
        #Identify the two blocks and swap them over. 
        block1 = np.where(id_array == 1) 
        block2 = np.where(id_array != 1) 
    
        #Swap over arrays if the max in idarray is more than 1. 
        bool_array = np.concatenate([bool_array[block2], bool_array[block1]]) 
        alpha = np.concatenate([alpha[block2], alpha[block1]]) 
        beta = np.concatenate([beta[block2], beta[block1]]) 
    else:
        pass 
        
    return bool_array, alpha, beta   
         
class fov_class():
    '''This class just replicates the original fov class, and just has the information added to it that transformations needs, i.e. the Cartesian and image unit vectors.''' 
    
    def __init__(self, xi, yi, zi):
        '''Parameters
        -------------
        xi - 3-element array describing the image x unit vector. 
        yi - 3-element array describing the image y unit vector. 
        zi - 3-element array describing the image z unit vector/LOS vector. 
        
        '''
        
        self.x0 = np.array([1,0,0])
        self.y0 = np.array([0,1,0])
        self.z0 = np.array([0,0,1]) 
        
        self.xi = xi
        self.yi = yi
        self.zi = zi  
