#This contains functions to add the Earth to plots. 

import numpy as np 
from matplotlib.patches import Wedge, Polygon, Circle

def make_earth(ax, rotation=0):
    '''This will add a little plot of the Earth on top for reference. This if for 2D   pictures. 
    
    Parameters
    ----------
    ax - axis object to add the earth too. 
    rotation - rotation of the Earth. def = 0. 
    '''

    # Add white circle first. 
    r=1
    circle = Circle((0,0), r, facecolor='w', edgecolor='navy')
    ax.add_patch(circle)

    # Add nightside. 
    theta2 = np.arange(181)-180+rotation
    xval2 = np.append(r*np.cos(theta2*(np.pi/180)),0)
    yval2 = np.append(r*np.sin(theta2*(np.pi/180)),0)
    verts2 = [[xval2[i],yval2[i]] for i in range(len(xval2))]

    polygon2 = Polygon(verts2, closed=True, edgecolor='navy', facecolor='navy', alpha=1) 
    ax.add_patch(polygon2)
    
def make_earth_3d(ax, color='k'):
    '''This will add a sphere for the Earth in a 3D plot. '''
        
    #Create a spherical surface. 
    radius = 1
    u = np.linspace(0, 2*np.pi, 100) 
    v = np.linspace(0, np.pi, 100) 
    x = radius* np.outer(np.cos(u), np.sin(v))
    y = radius* np.outer(np.sin(u), np.sin(v))
    z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, color=color, lw=0, alpha=1)
    
def make_earth_3d_2(ax):
    '''This will add a sphere for the Earth. It is black one side, white the other.'''
        
    #Create a spherical surface. 
    radius = 1
    u = np.linspace(np.pi/2, 1.5*np.pi, 100) 
    v = np.linspace(0, np.pi, 100) 
    x = radius* np.outer(np.cos(u), np.sin(v)) 
    y = radius* np.outer(np.sin(u), np.sin(v))
    z = radius* np.outer(np.ones(np.size(u)), np.cos(v))
        
    ax.plot_surface(x, y, z, color='k', lw=0, alpha=1)
        
    u = np.linspace(-np.pi/2, np.pi/2, 100) 
    v = np.linspace(0, np.pi, 100) 
    x = radius* np.outer(np.cos(u), np.sin(v)) 
    y = radius* np.outer(np.sin(u), np.sin(v))
    z = radius* np.outer(np.ones(np.size(u)), np.cos(v))
        
    ax.plot_surface(x, y, z, color='w', lw=0, alpha=1) 
