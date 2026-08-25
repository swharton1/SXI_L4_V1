#This runs the code to read in the coastlines file. 

import geopandas 
import matplotlib.pyplot as plt 
import numpy as np 
import warnings 
warnings.simplefilter('ignore', UserWarning)
from spacepy import coordinates as coord
from spacepy.time import Ticktock 
import datetime as dt 

class coastlines():
    '''This class will contain functions to read in the coastlines data and manipulate it.''' 
    
    def __init__(self, filename='/home/s/sw682/Code/geomaps/ne_110m_coastline.shp', time=dt.datetime(2020,1,1)):
        '''This reads in the coastline data.''' 
        
        #Read in the data. 
        self.filename = filename 
        print ('Read in Coastlines file...') 
        self.world = geopandas.read_file(self.filename) 
        self.time = time
        
        #Extract the latitutes and longitudes.
        print ('Extract coastline longitudes and latitudes...')  
        self.get_coordinates()
        
        #Convert coastlines to GSE. 
        print ('Convert coastlines to GSE...') 
        self.x, self.y, self.z = self.convert_geo_to_gse(self.r, self.lat, self.lon) 
        
        #Now get the parallels and convert them to GSE too.
        print ("Get the Earth's parallels...") 
        self.get_parallels() 
        self.para_x, self.para_y, self.para_z = self.convert_geo_to_gse(self.para_r, self.para_lat, self.para_lon) 
        
        #Now get the meridians and convert them to GSE. 
        print ("Get the Earth's meridians...") 
        self.get_meridians()
        self.meri_x, self.meri_y, self.meri_z = self.convert_geo_to_gse(self.meri_r, self.meri_lat, self.meri_lon) 
        
    def plot_world_map(self, color='lightgray', edgecolor='k', alpha=0.5):
        '''This plots a basic world map with geopandas.''' 
        
        fig = plt.figure()
        ax = fig.add_subplot(111) 
        
        self.world.plot(ax=ax, color=color, edgecolor=edgecolor, alpha=alpha) 
    
    
    def plot_world_map_2(self):
        '''This plots a world map but using my extraction of the coordinates, instead of the built in geopandas function.''' 
        
        plt.style.use('classic')
        fig = plt.figure()
        ax = fig.add_subplot(111) 
        
        #How to plot coastlines. 
        unique = np.unique(self.indices)
        for u in unique: 
            i = self.indices == u 
            ax.plot(self.lon[i], self.lat[i], color='grey') 
        
        #How to plot parallels. 
        n=5 
        for i in range(n):
            ax.plot(self.para_lon[361*i:361*(i+1)], self.para_lat[361*i:361*(i+1)], color='k') 
        
        #How to plot meridians. 
        n=12
        for i in range(n):
            ax.plot(self.meri_lon[181*i:181*(i+1)], self.meri_lat[181*i:181*(i+1)], color='k') 
        
        ax.set(xlim=(-180,180), ylim=(-90,90))      
        ax.set_aspect('equal') 
    
    def plot_earth_3d(self):
        '''This will plot the Earth in 3D in the GSE coordinate system.''' 
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d') 
        
        #Create a spherical surface. 
        radius = 0.99
        u = np.linspace(0, 2*np.pi, 100) 
        v = np.linspace(0, np.pi, 100) 
        x = radius* np.outer(np.cos(u), np.sin(v))
        y = radius* np.outer(np.sin(u), np.sin(v))
        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color=None, lw=0, alpha=0.1)
        
        unique = np.unique(self.indices) 
        for u in unique: 
            i = self.indices == u 
            ax.plot(self.x[i], self.y[i], self.z[i], color='k', lw=0.5) 
    
        #Labels. 
        ax.set(xlabel=r'X$_{GSE}$', ylabel=r'Y$_{GSE}$', zlabel=r'Z$_{GSE}$') 
        ax.set_aspect('equal') 
        
        
        
        
        
    def get_coordinates(self):
        '''This gets the longitudes and latitudes of all the points in the coastlines map.''' 
        
        #This returns the coordinates into a pandas core data frame. 
        #Longitude in deg = column 'x', Latitude in deg = column 'y'. 
        coords = self.world.get_coordinates() 
        
        #Separate this into the different coast objects. 
        self.indices = np.array(coords.index) 
        self.lon = np.array(coords['x'])
        self.lat = np.array(coords['y'])
        self.r = np.ones(len(coords['x']))
    
    def get_parallels(self):
        '''This gets some of the parallels.''' 
        
        #Define the angles for the Earth's parallels. 
        para_angle = np.array([-60, -30, 0, 30, 60])
        phi = np.linspace(-180, 180, 361)  
        
        PHI, PARA = np.meshgrid(phi, para_angle) 
        
        #Flatten to get 1D arrays of lons and lats.
        #It switches parallel every 361 values.  
        self.para_lon = PHI.flatten()
        self.para_lat = PARA.flatten() 
        self.para_r = np.ones(len(self.para_lat))     
    
    
    def get_meridians(self):
        '''This gets some of the meridians.''' 
        
        #Define the angles for the Earth's meridians. 
        meri_angle = np.array([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]) 
        lat = np.linspace(-90,90,181)
        
        LAT, MERI = np.meshgrid(lat, meri_angle) 
        
        #Flatten to 1D arrays of lons and lats. 
        #It switches meridian every 181 values. 
        self.meri_lon = MERI.flatten()
        self.meri_lat = LAT.flatten()
        self.meri_r = np.ones(len(self.meri_lat))   
         
    def convert_geo_to_gse(self, r, lat, lon):
        '''This uses spacepy to convert from geographic coordinates to GSE.'''
        
        #Package up the geographic coordinates correctly. 
        coords_start = np.zeros((len(r), 3))
        coords_start[:,0] = r
        coords_start[:,1] = lat
        coords_start[:,2] = lon
        
        #Make the coordinate object. 
        coord_obj = coord.Coords(coords_start, 'GEO', 'sph') 
        
        #Add the time. 
        time_array = [self.time for i in range(len(r))]
        coord_obj.ticks = Ticktock(time_array, 'UTC') 
        
        #Convert to GSE. 
        coords_end = coord_obj.convert('GSE', 'car') 
        
        #Extract the x, y and z positions in GSE. 
        return coords_end.x, coords_end.y, coords_end.z

