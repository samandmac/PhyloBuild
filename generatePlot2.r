#!/usr/bin/env Rscript
#Because we added in the user input options, we need to change file locations to
#the user choices by adding in parameters.
args = commandArgs(trailingOnly = TRUE)

#Installation of treedataverse, which contains all the packages we need
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("remotes", quietly = TRUE))
    BiocManager::install("remotes")   
if (!require("YuLab-SMU/treedataverse", quietly = TRUE))
    BiocManager::install("YuLab-SMU/treedataverse")
if (!require("devEMF", quietly = TRUE))
    install.packages("devEMF")


#Loading said packages

library(treedataverse)
library(devEMF)

save_plot<- function(ggp,plot_height,plot_width,plot_path)
{
  emf(plot_path, height=plot_height/90, width=plot_width/90)
  print(ggp)
  dev.off()
  #Files must be saved as EMF file
  # clears all devices - as a safety measure
  while (dev.cur()>1) dev.off()
}

#####CIRCLE AND HEATMAP#####
#last examples used p as a loaded in tree, obviously in past examples we had the colour of the tip
#names, which was because we had added in the grouping as data in the tree.tbl. This is an attempt
#to have both the colour and the heatmap.
#This is all the previous code, just with a new bit about adding a heatmap onto what we had already
#Loading data
tree = read.newick(paste(args[1],"/phylipFor.phy_phyml_tree.txt", sep = ""),
                   node.label = "support")
get.tree(tree)
#Converting to tibble for adding in phylogroup data
tree.tbl = as_tibble(tree)
#Removing .blastn at end of name so that the phylogroup and tree use same label names
tree.tbl$label = gsub(".blastn", "", tree.tbl$label)
#Loading phylogrouping from file, the gene matches, gene names, and genomes added
phylogroup = read.table(paste(args[1],"/phylogeny_list.txt", sep = ""),
                        col.names = c("group","label"))
phylogroup$group = as.factor(phylogroup$group)
phylogroup$group = factor(phylogroup$group, levels = c("A", "B1","B2","C","D",
                                                       "E","F","G","cladeI","cladeII",
                                                       "cladeIII","cladeIV","cladeV",
                                                       "albertii","fergusonii"))

geneMatches=read.table(paste(args[1],"/geneMatches.txt", sep=""),
                       col.names = c("gene","label"))
genes = read.table(paste(args[1],"/genesOI.txt", sep=""))
geneOrder = as.vector(unlist(genes))
geneMatches$gene = as.factor(geneMatches$gene)
geneMatches$gene = factor(geneMatches$gene, levels = geneOrder)

levels(geneMatches$gene)
levels(phylogroup$group)

genomes_added_in = read.table(paste(args[1],"/genomes.added.txt", sep=""), 
                              col.names = c("genomes"))
#Joining phylogroup into tree.tbl
tree.tbl = full_join(tree.tbl, phylogroup, by = 'label')
#tree.tbl = full_join(tree.tbl, geneMatches, by = 'label')
#Making the tree.tbl back into a treedata object
tree = as.treedata(tree.tbl)
#Manipulating new dataframe, just the phylogroup information and row names as labels
phylogroup1 = phylogroup[1]
row.names(phylogroup1) = phylogroup[,2]
#Plotting the tree, circular form, SAVED AS TREE.PLOT
# tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE)+
#   geom_tiplab(size=3, aes(colour=group), offset=4)+
#   scale_color_discrete(breaks = c("A", "B1","B2","C","D","E","F","G","cladeI","cladeII",
#                                   "cladeIII","cladeIV","cladeV","albertii","fergusonii"))+
#   guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))#+
# #Adding the heatmap onto the tree plot
# plot_heatmap = gheatmap(tree.plot, phylogroup1, offset=0, width=0.1, font.size=3, colnames = FALSE)+
#   scale_fill_discrete(breaks=c("A", "B1","B2","C","D","E","F","G","cladeI","cladeII",
#                                "cladeIII","cladeIV","cladeV","albertii","fergusonii"), 
#                       name="Phylogrouping")

#HEATMAP FOR GENE CARRIAGE
datalist = list()
for (i in seq_along(levels(geneMatches$gene))) {
  if(i <= length(levels(geneMatches$gene))) {
    pgroup = levels(geneMatches$gene)[i]
    subs = subset(geneMatches, gene == pgroup)
    subs = distinct(subs)
    #NOTE we use phylogroup as the join and not the original dataset, as we need the 
    #label to NOT repeat itself!
    new = full_join(subs, phylogroup, by="label")
    row.names(new) = new[,2]
    new = new[1]
    names(new) = pgroup
    datalist[[i]] = new
  }
}

heatmap = Reduce(merge, lapply(datalist, function(x) data.frame(x, rn = row.names(x))))
row.names(heatmap) = heatmap[,1]
heatmap = heatmap[,c(2:24)]

#This allows us to select the nodes for the genomes added in, for altering the plot!
merge_of_genomes_added_in = merge(tree.tbl, genomes_added_in, by.x = 4, by.y = 1)
nodes_of_above =unlist(merge_of_genomes_added_in[,3])
nodes_of_above

#Plotting the tree, circular form, SAVED AS TREE.PLOT
#Added in colour=group to ggtree to change the branch colour to group 
tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, aes(colour=group), show.legend=F)+ #show.legend =F to get rid of line in legend
  #Changed geom_tiplab to be only the genomes added in
  #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_of_above))+
  #Alternative, highlight branches of genomes added in, and give all the genomes names as tip labels?
  geom_tiplab(size=3, aes(colour=group), offset=26)+
  geom_hilight(mapping=aes(subset = node %in% nodes_of_above, fill = "red"))+
  scale_color_discrete("Group",breaks = c("A", "B1","B2","C","D","E","F","G","cladeI","cladeII",
                                          "cladeIII","cladeIV","cladeV","albertii","fergusonii"))+
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))#+

#This adds the heatmap onto the plot
plot_heatmap = gheatmap(tree.plot, heatmap, offset=0, width=1, font.size=3, colnames = FALSE)+
  scale_fill_discrete(breaks=geneOrder, 
                      name="Gene Carriage")

save_plot(plot_heatmap, 1000, 1000, paste(args[1],"/finalPlot.EMF", sep =""))
#If using multiple heatmaps (which is possible, so could do it for each gene of interest?) you can
#either have multiple columns in the dataframe you are using or add another heatmap on after using
#new_scale_fill() e.g.:
#p2 <- p1 + new_scale_fill()
# gheatmap(p2, df2, offset=15, width=.3,
#          colnames_angle=90, colnames_offset_y = .25) +
#   scale_fill_viridis_c(option="A", name="continuous\nvalue")

#You can also visualise a multiple sequence alignment, but maybe not necessary here, I reckon that
#would be too cumbersome a plot - considering it already gets a bit heavy

#Visualisation with images is also possible - it's pointless for our stuff! :)

#Can use other tree-like objects, not so relevant for what we do

#ggextra allows for extra features on a circular plot, but none are really relevant to what we do -
#I think what we have at the moment could be fine as is?

#####Example feature list#####
#Features specifically added for these packages
# geom_balance	highlights the two direct descendant clades of an internal node
# geom_cladelab	annotate a clade with bar and text label (or image)
# geom_facet	plot associated data in a specific panel (facet) and align the plot with the tree
# geom_hilight	highlight selected clade with rectangular or round shape
# geom_inset	add insets (subplots) to tree nodes
# geom_label2	the modified version of geom_label, with subset aesthetic supported
# geom_nodepoint	annotate internal nodes with symbolic points
# geom_point2	the modified version of geom_point, with subset aesthetic supported
# geom_range	bar layer to present uncertainty of evolutionary inference
# geom_rootpoint	annotate root node with symbolic point
# geom_rootedge	add root edge to a tree
# geom_segment2	the modified version of geom_segment, with subset aesthetic supported
# geom_strip	annotate associated taxa with bar and (optional) text label
# geom_taxalink	Linking related taxa
# geom_text2	the modified version of geom_text, with subset aesthetic supported
# geom_tiplab	the layer of tip labels
# geom_tippoint	annotate external nodes with symbolic points
# geom_tree	tree structure layer, with multiple layouts supported
# geom_treescale	tree branch scale legend
