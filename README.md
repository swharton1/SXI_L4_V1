# SXI_L4_V1

S. J. Wharton 
August 2026 

This is my code to plot the L4 output files from CMEM. It is V1 as it is based on the fits to the GAC files. 

The code can be used to produce one of my summary plots. These display: 
Six example images 
A main subplot of the magnetopause position and uncertainty with the SXI FOV. 
A subplot of the foreground and background count rates. 
Three subplots of the spacecraft position in each Cartesian plane. 

Optional subplot to show the quality flag decomposition. 
Optional subplot to show the output of the Bayane cusp ID. 

All the required functions are present. You don't need any of my other libraries. It just requires standard numpy, matplotlib, glob, os, datetime and string libraries which come by default. 

Users should run the example.py script to see how to run the code. E.g.:
python example.py 

Users can make their own adapted scripts to run it if they wish. 

