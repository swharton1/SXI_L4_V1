#This is an example script for how to run the summary plot. 
#You can either run this script or copy this code into ipython. 

import SXI_L4_V1

#Specify the path where your fitted L4 files are stored. 
#It will read in all of the ones that fit wildcards! 
path='/data/smile/shared/sims/GAC_V2_Sam/GAC_V2.2_1x1_300s/'

#Type of fit. Use 'original' for method as described in Wharton et al. (2026). 
tag='original' 

#Load in the L4 files and extract the data. 
analyse = SXI_L4_V1.analyse_GAC.analyse_gac(path=path, tag=tag) 

#Make the summary plot. Set figname = None to use default name. Plot will appear in plots directory.  
analyse.plot_results(save=True, add_qf=True, add_cusp=False, figname=None) 
