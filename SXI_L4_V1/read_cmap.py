from matplotlib.colors import ListedColormap
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) 
import os

def txt2matplotlib(cbpath=os.path.dirname(__file__), cbname='smile_colt_lundi.txt'):
    '''
    Example:
    cbpath = '/data/smile/sxi_coltab/'
    cbname = 'smile_colt_lundi.txt'
    txt2matplotlib(cbpath, cbname)

    '''
    # Read the file and load it as a DataFrame
    print (cbpath)
    cb = pd.read_csv(os.path.join(cbpath +'/'+ cbname), header=None, sep=r"\s+")
    
    # Create a new column by converting RGB values to tuples
    cb['color'] = cb.apply(lambda row: tuple(row[:3] / 255), axis=1)

    # Generate the custom color map
    custom_cmap = ListedColormap(cb['color'].tolist())
    return custom_cmap
