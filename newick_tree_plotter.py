#!/usr/bin/env python3

"""
Module for parsing newick phylogentic tree files and outputting matplotlib-based phylogeny figure.
Primarily tested for ~100 or fewer terminal nodes (i.e. strains, leaves, etc.)
Makes heavy use of Biopython's Phylo module.
"""

#Imported modules
from Bio import Phylo
from matplotlib import rc, patches
import matplotlib.pyplot as plt
from collections import defaultdict

	
# def inspect_class_object(class_object):
# 	"""
# 	Legacy code that was useful in understanding the tree and clade class objects.
# 	Leaving it here for future reference.
# 	"""
# 	import inspect
# 	print(inspect.getmembers(class_object, lambda a:not(inspect.isroutine(a))))


def none_label(leaf):
	"""
	Just takes label names from plot function and returns None for each to suppress Phylo output.
	"""
	return None


def get_XY_positions(tree):
	"""
	Parse tree object to extract the exact x,y coords for each node/leaf.
	https://github.com/biopython/biopython/blob/a09f38b468cbedaf25aec6cfed1475b0849011e6/Bio/Phylo/_utils.py#L320
	"""

	#Depths is equal to x coords i.e. sum of branch lengths leading to node
	depths = tree.depths()
	# If there are no branch lengths, assume unit branch lengths
	if not max(depths.values()):
		depths = tree.depths(unit_branch_lengths=True)


	#Heights is equal to Y coords. Should be integers for leaves.
	maxheight = tree.count_terminals()
	# Rows are defined by the tips
	heights = {
		tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
	}

	# Internal nodes: place at midpoint of children
	def calc_row(clade):
		for subclade in clade:
			if subclade not in heights:
				calc_row(subclade)
		# Closure over heights
		heights[clade] = (
			heights[clade.clades[0]] + heights[clade.clades[-1]]
		) / 2.0

	if tree.root.clades:
		calc_row(tree.root)
	
	#Combine X and Y coordinate dictionaries above into one dictionary with X,Y as list.
	combined_dict = defaultdict(list)
	dict_list = [depths, heights]
	for dictionary in dict_list:
		for key, value in dictionary.items():
			combined_dict[key].append(value)
	
	
	return(combined_dict)





def plot_tree(newick_tree_file, output_file, sample_name, plot_title, low_boot=75, high_boot=95):
	"""
	Primary function to parse tree data and plot elements.
	With specific logic to apply text colors, bootstrap markers, etc.
	"""
	
	class TreeException(Exception):
		pass

	
	if low_boot >= high_boot or high_boot > 99:
		raise TreeException("Bootstrap cutoffs are wrong")
	
	
	#Import tree file
	tree = Phylo.read(newick_tree_file, "newick")
	#Arrange from low to high branch depths (overall tree is cleaner and easier to read)
	tree.ladderize()
	
	#Set global plot formats
	rc('font', size=6)
	plt.rcParams["font.sans-serif"] = "Arial"
	
	#Original idea was to plot boxes and labels all on seperate axes, but I don't think that will be necessary.
	#I'm leaving the code as a 1x1 subplot though in case I cange my mind later.
	fig, (ax0) = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(10, 20))

	#Axes customizations
	#Pretty much just hiding all the unused stuff.
	#See plot notes below for more detail on how to overcome certain Phylo defaults.
	ax0.set_title(plot_title, fontsize=16)
	
	ax0.tick_params(axis='x', labelsize=10)
	ax0.spines['top'].set_visible(False)
	ax0.spines['right'].set_visible(False)
	ax0.spines['left'].set_visible(False)
	ax0.axes.get_yaxis().set_visible(False)
	
	#I originally tried to make the axis clean by forcing ticks. But for highly related alignments, this was way to large.
	#I'm letting the plotter handle axes scaling for now, but leaving this if I need to revisit.
# 	ax0.spines.bottom.set_bounds((0, 1))
# 	ax0.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
	
	
	#Parse node info (Pylo calls them clades)
	clades_dict = get_XY_positions(tree)
	max_y = 0
	for key in clades_dict:
		clade_name = key.name
		clade_branch = key.branch_length
		clade_confidence = key.confidence
		x_value = clades_dict[key][0]
		y_value = clades_dict[key][1]
		
		#Get total vertical size of tree plot to place legend info
		#Y axis is inverted for some reason in Phylo (zero at top)
		if y_value > max_y:
			max_y = y_value

		#This labels the leaves (internal nodes have name 'None')		
		if key.name != None:
			if clade_name == sample_name:
				ax0.text(x_value+0.005, y_value, s=clade_name, zorder=20, color='red', fontsize=8, ha='left', va='center')
			else:
				ax0.text(x_value+0.005, y_value, s=clade_name.replace('_', ' '), zorder=20, color='black', fontsize=6, ha='left', va='center')
		else:
			pass
		
		#Logic to mark nodes by bootstrap values
		if clade_confidence != None:
			if low_boot <= clade_confidence < high_boot:
				ax0.plot(x_value, y_value, marker="o", markersize=4, zorder=50, color='black')
			elif high_boot <= clade_confidence:
				ax0.plot(x_value, y_value, marker="o", markersize=4, zorder=50, markerfacecolor='white', markeredgecolor='black', markeredgewidth=1)

		#Hack to fill in gaps where Phylo draws lines that intersect at the middle, leaving an ugly notch at each node.
		#Circle marker 'o' looks a little nicer than square, but has a little artifact at high zoom.
		if clade_branch != None:
			ax0.plot(x_value-clade_branch, y_value, marker="s", markersize=.5, color="black", zorder=1)
	
	
	#Create legend in bottom left (should be space there due to ladderize)
	ax0.plot(0, max_y-1, marker="o", markersize=8, zorder=50, color='black')
	ax0.text(0.01, max_y-1, s=f"Bootstrap support >{low_boot}", zorder=20, color='black', fontsize=8, ha='left', va='center')
	
	ax0.plot(0, max_y, marker="o", markersize=8, zorder=50, markerfacecolor='white', markeredgecolor='black', markeredgewidth=1)
	ax0.text(0.01, max_y, s=f"Bootstrap support >{high_boot}", zorder=20, color='black', fontsize=8, ha='left', va='center')
	
# 	plt.Rectangle((0, max_y), width=.03, height=max_y-1, linewidth=1, color='black', zorder=50)
# 	ax0.add_patch(patches.Rectangle((-0.006, max_y+0.6), 0.083, -2.25, facecolor='None', edgecolor='black', linewidth=1))
	
	#Draw tree. Phylo hard codes a few things including the axes labels which is stupid and took a while to sort out.
	#You can overwrite with kwargs, but it uses a tricky syntax.
	#https://github.com/biopython/biopython/blob/a09f38b468cbedaf25aec6cfed1475b0849011e6/Bio/Phylo/_utils.py#L474
	Phylo.draw(tree, 
		axes=ax0, 
		do_show=False, 
		label_func=none_label, 
		branch_labels=none_label, 
		xlabel=({'xlabel':'Branch Distance\n(nucleotide substitutions per site)', 'fontsize':14}))
		
	fig.savefig(output_file,bbox_inches='tight')

	return(True)




if __name__ == "__main__":
	import sys
	
	newick_file = sys.argv[1]
	plot_tree(f"{newick_file}", f"{newick_file}.pdf", sample_name=None, plot_title='Title text')


