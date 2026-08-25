from matplotlib.patches import Rectangle

def add_label(ax, text='', corner='topleft', width = 0.1, height = 0.1, 
	backcolour='k', textcolour='w',fontsize=12, zorder=1, alpha=1):
	'''This will add a figure label to a subplot. 
	
	Parameters
	----------
	ax - axis to add the label too. 
	text - text to put in the box. def = ''
	corner - 'topright' (def), 'topleft', 'bottomleft', 'bottomright'
	width - label width as a fraction of the plot, def = 0.1
	height - label width as a fraction of the plot, def = 0.1
	backcolour - label colour, def = 'k'
	textcolour - text colour, def = 'w'
	fontsize - size of the text, def = 12
	zorder - order in which to lay down objects onto the axis. def = 3, 
		which should put it on top of other stuff. 
	alpha - transparency
	 
	'''

	
	if corner == 'topright':
		cnr = (1,1)
		width = -width
		height = -height
		va = 'top'
		ha = 'right'
		#text_loc = (cnr[0]-0.01, cnr[1]-0.01)
	elif corner == 'topleft':
		cnr = (0,1)
		width = width
		height = -height
		va = 'top'
		ha = 'left'
		#text_loc = (cnr[0]+0.01, cnr[1]-0.01)
	elif corner == 'bottomleft':
		cnr = (0,0)
		width = width
		height = height
		va = 'bottom'
		ha = 'left'
		#text_loc = (cnr[0]+0.01, cnr[1]+0.01)
	elif corner == 'bottomright':
		cnr = (1,0)
		width = -width
		height = height
		va = 'bottom'
		ha = 'right' 
		#text_loc = (cnr[0]-0.01, cnr[1]+0.01)
		
	
	#textcentre=(cnr[0]+(width/2.0),cnr[1]+(height/2.0)) 	
	
	ax.add_patch(Rectangle(cnr, width, height, facecolor=backcolour, edgecolor=backcolour,
		transform=ax.transAxes, zorder=zorder, alpha=alpha))
	ax.text(*cnr, text, color=textcolour, va=va,
		ha=ha, transform=ax.transAxes, fontsize=fontsize, zorder=zorder+1) 
