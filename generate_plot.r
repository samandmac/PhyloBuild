#!/usr/bin/env Rscript

#Because we added in the user input options, we need to change file locations to
#the user choices by adding in parameters.
args = commandArgs(trailingOnly = TRUE) #We just run Rscript with a path to the output files which contains files for plots.

#Installation of treedataverse, which contains all the packages we need
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#if (!require("remotes", quietly = TRUE))
#    BiocManager::install("remotes")   
#if (!require("YuLab-SMU/treedataverse", quietly = TRUE))
#    BiocManager::install("YuLab-SMU/treedataverse")
#if (!require("devEMF", quietly = TRUE))
#    install.packages("devEMF")
#if (!require("data.table", quietly = TRUE))
#  install.packages("data.table")

#install.packages("data.table")
#Loading said packages

library(treedataverse)
library(devEMF)
library(data.table)

#This is the function for saving plots as an EMF image.
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
#This is where we make the plot itself, containing the tree, labelled tips, a heatmap of carriage,
#labels on the nodes to indicate phylogroup groups.

#The phylipFor.phy_phyml_tree.txt is in a newick file format, so we can easily load it in as below.
tree = read.newick(paste(args[1],"/phylipFor.phy_phyml_tree.txt", sep = ""),
                   node.label = "support")

#This just gives some information about the tree
get.tree(tree)

#Converting to tibble for adding in phylogroup data - easier than adjusting tree object.
tree.tbl = as_tibble(tree)

#Removing .blastn at end of name so that the phylogroup and tree use same label names
tree.tbl$label = gsub(".blastn", "", tree.tbl$label)

groups_yn = args[3]

if (groups_yn == "yes"){
  #Loading phylogrouping from file, the gene matches, gene names, and genomes added
  phylogroup = read.table(paste(args[1],"/../group_list.txt", sep = ""),
                          col.names = c("group","label"))
  phylogroup$group = as.factor(phylogroup$group) #Changing the grouping to be a factor
}

#Because we want phylogroup to become a generic group for indicating strain group, we have to take this out - leave commented out so user can look at script and decide factor level?
#phylogroup$group = factor(phylogroup$group, levels = c("A", "B1","B2","C","D",
#                                                      "E","F","G","I","II",
#                                                     "III","IV","V",
#                                                    "ALB","FER")) #Changing the order of factors

gene_names = args[4]

if (gene_names == "yes"){
  new_gene_names = read.table(paste(args[1],"/../new_gene_names.txt", sep = ""), sep = "\t")
  names(new_gene_names) = c("old_names","new_names")
}

geneMatches=read.table(paste(args[1],"/geneMatches.txt", sep=""),
                       col.names = c("gene","label"))
genes = read.table(paste(args[1],"/genesOI.txt", sep=""))
geneOrder = as.vector(unlist(genes))
geneMatches$gene = as.factor(geneMatches$gene)
geneMatches$gene = factor(geneMatches$gene, levels = geneOrder)

#levels(geneMatches$gene)
#levels(phylogroup$group)

genomes_added_in = read.table(paste(args[1],"/genomes.added.txt", sep=""), 
                              col.names = c("genomes"))

renaming_genomes = args[2]

#This part just finds if the user replied yes or no and gives a difFERent phylogroup list which is used later in a loop.
if (renaming_genomes == "yes") {
  #This little bit of script can be used to change names if required. If you do, place the "new_names" file into the output folder.
  new_genome_names = read.table(paste(args[1],"../new_genome_names.txt",sep=""),sep = "\t")
  names(new_genome_names) = c("label","new_names")
  lb = get.tree(tree)$tip.label
  d = data.frame(label = lb, label2=lb)
  d = merge(d,new_genome_names, by="label", all.x = TRUE)
  for (x in 1:length(d$label)){
    print(x)
    if (is.na(d$new_names[x]) == FALSE){
      d$label2[x] = d$new_names[x]
    }
  }
  d = d[,-3]
} else {
  lb = get.tree(tree)$tip.label
  d = data.frame(label = lb, label2=lb)
}

if (groups_yn == "yes"){
  #Joining phylogroup into tree.tbl
  tree.tbl = full_join(tree.tbl, phylogroup, by = 'label')
  #tree.tbl = full_join(tree.tbl, geneMatches, by = 'label') # Not needed.
  
  #Making the tree.tbl back into a treedata object
  tree = as.treedata(tree.tbl)
  
  #Manipulating new dataframe, just the phylogroup information and row names as labels
  phylogroup1 = phylogroup[1]
  row.names(phylogroup1) = phylogroup[,2]
}

#HEATMAP FOR GENE CARRIAGE
#Below makes the dataset for the heatmap for the tree plot

#If grouping 
if (groups_yn == "yes"){
  datalist = list() #Makes empty list for filling in loop
  for (i in seq_along(levels(geneMatches$gene))) { #For each gene count i
    
    if(i <= length(levels(geneMatches$gene))) { #If i is less than the length of the genes
      
      gene_name = levels(geneMatches$gene)[i] #gene_name = the name of the gene match
      subs = subset(geneMatches, gene == gene_name) #subs = subset of geneMatches where the gene matches the gene_name
      subs = distinct(subs) #subs = only distinct subs, so no repeats
      new = full_join(subs, phylogroup, by="label") #new = a join of subs and phylogroup to obtain a dataset containing phylogrouping
      #and NA's where there is no gene match present.
      row.names(new) = new[,2] #Changing the row names to the genome names
      new = new[1] #new = only the column giving the gene name that matched genomes - NA's present if there was no match.
      names(new) = gene_name #changes the column name to gene_name
      datalist[[i]] = new #We add this to the data list that was empty - so for each gene tested for carriage this repeats. Gives a
      #column of gene names (in a row for each genome) if there was a match, NA's if there was no match, and 
      #saves results in the datalist for use later.
    }
  }
  #Next we combine all the data saved in the datalist into one data table.
  heatmap = Reduce(merge, lapply(datalist, function(x) data.frame(x, rn = row.names(x))))
  row.names(heatmap) = heatmap[,1] #CHange the row names to the genome names
  heatmap = heatmap[,-1] #Remove the genome name as a column.
}

#If not grouping
if (groups_yn == "no"){
  datalist = list() #Makes empty list for filling in loop
  for (i in seq_along(levels(geneMatches$gene))) { #For each gene count i
    
    if(i <= length(levels(geneMatches$gene))) { #If i is less than the length of the genes
      
      gene_name = levels(geneMatches$gene)[i] #gene_name = the name of the gene match
      subs = subset(geneMatches, gene == gene_name) #subs = subset of geneMatches where the gene matches the gene_name
      subs = distinct(subs) #subs = only distinct subs, so no repeats
      new = full_join(subs, tree.tbl, by="label") #new = a join of subs and the tree.tbl to obtain a dataset containing labels and NA's where there is no gene match present.
      new = new[,-c(3:6)] #We need to remove columns generated from the join that aren't in the phylogrouping version
      new = distinct(new) #Filter to only unique labels
      new = head(new,-1) #Remove the node that doesn't contain any data
      row.names(new) = new[,2] #Changing the row names to the genome names
      new = new[1] #new = only the column giving the gene name that matched genomes - NA's present if there was no match.
      names(new) = gene_name #changes the column name to gene_name
      datalist[[i]] = new #We add this to the data list that was empty - so for each gene tested for carriage this repeats. Gives a
      #column of gene names (in a row for each genome) if there was a match, NA's if there was no match, and 
      #saves results in the datalist for use later.
    }
  }
  #Next we combine all the data saved in the datalist into one data table.
  heatmap = Reduce(merge, lapply(datalist, function(x) data.frame(x, rn = row.names(x))))
  row.names(heatmap) = heatmap[,1] #CHange the row names to the genome names
  heatmap = heatmap[,-1] #Remove the genome name as a column.
}

#This allows us to select the nodes for the genomes added in, for altering the plot!
merge_of_genomes_added_in = merge(tree.tbl, genomes_added_in, by.x = 4, by.y = 1) #We previously made genomes added in, so we can just
#merge them by their label name
nodes_added_in =unlist(merge_of_genomes_added_in[,3]) #Saving the nodes of the genomes_added_in
nodes_added_in

if (groups_yn == "yes"){
  if (gene_names == "no"){
    #Below is the script for when pgroup labels aren't required.
    #Back to making the actual plot now, and saving that.
    #Added in colour=group to ggtree to change the branch colour to group 
    tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, aes(colour=group), show.legend=F) %<+% d + #show.legend =F to get rid of line in legend
      #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_added_in))+ #Changed geom_tiplab to be only the genomes added in
      #Alternative plot below, highlight branches of genomes added in, and give all the genomes names as tip labels?
      geom_tiplab(size=3, aes(colour=group, label=label2), offset=18)+ #changed offset to 28 from 26 for Clermont genes
      geom_hilight(mapping=aes(subset = node %in% nodes_added_in, fill = "red"))+ #Highlights nodes of genomes that were added in
      scale_color_discrete("Group")+ #This is just another way of getting the legend in order
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    #This adds the heatmap onto the plot
    plot_heatmap = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
      scale_fill_discrete(breaks=geneOrder, 
                          name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))
    
    #And finally saving that plot.
    save_plot(plot_heatmap, 1500, 1500, paste(args[1],"/finalPlot.EMF", sep =""))
  }
  
  if (gene_names == "yes"){
    #Below is the script for when pgroup labels aren't required.
    #Back to making the actual plot now, and saving that.
    #Added in colour=group to ggtree to change the branch colour to group 
    tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, aes(colour=group), show.legend=F) %<+% d + #show.legend =F to get rid of line in legend
      #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_added_in))+ #Changed geom_tiplab to be only the genomes added in
      #Alternative plot below, highlight branches of genomes added in, and give all the genomes names as tip labels?
      geom_tiplab(size=3, aes(colour=group, label=label2), offset=18)+ #changed offset to 28 from 26 for Clermont genes
      geom_hilight(mapping=aes(subset = node %in% nodes_added_in, fill = "red"))+ #Highlights nodes of genomes that were added in
      scale_color_discrete("Group")+ #This is just another way of getting the legend in order
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    #This adds the heatmap onto the plot
    plot_heatmap = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
      scale_fill_discrete(breaks=geneOrder, labels = c(new_gene_names$new_names),
                          name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))
    #And finally saving that plot.
    save_plot(plot_heatmap, 1500, 1500, paste(args[1],"/finalPlot.EMF", sep =""))
  }
}

if (groups_yn == "no"){
  if (gene_names == "no"){
    #Below is the script for when pgroup labels aren't required.
    #Back to making the actual plot now, and saving that.
    #Added in colour=group to ggtree to change the branch colour to group 
    tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, show.legend=F) %<+% d + #show.legend =F to get rid of line in legend
      #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_added_in))+ #Changed geom_tiplab to be only the genomes added in
      #Alternative plot below, highlight branches of genomes added in, and give all the genomes names as tip labels?
      geom_tiplab(size=3, aes(label=label2), offset=18)+ #changed offset to 28 from 26 for Clermont genes
      geom_hilight(mapping=aes(subset = node %in% nodes_added_in, fill = "red"))+ #Highlights nodes of genomes that were added in
      scale_color_discrete("Group")+
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    #This adds the heatmap onto the plot
    plot_heatmap = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
      scale_fill_discrete(breaks=geneOrder, 
                          name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))
    
    #And finally saving that plot.
    save_plot(plot_heatmap, 1500, 1500, paste(args[1],"/finalPlot.EMF", sep =""))
  }
  if (gene_names == "yes"){
    #Below is the script for when pgroup labels aren't required.
    #Back to making the actual plot now, and saving that.
    #Added in colour=group to ggtree to change the branch colour to group 
    tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, show.legend=F) %<+% d + #show.legend =F to get rid of line in legend
      #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_added_in))+ #Changed geom_tiplab to be only the genomes added in
      #Alternative plot below, highlight branches of genomes added in, and give all the genomes names as tip labels?
      geom_tiplab(size=3, aes(label=label2), offset=18)+ #changed offset to 28 from 26 for Clermont genes
      geom_hilight(mapping=aes(subset = node %in% nodes_added_in, fill = "red"))+ #Highlights nodes of genomes that were added in
      scale_color_discrete("Group")+
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    #This adds the heatmap onto the plot
    plot_heatmap = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
      scale_fill_discrete(breaks=geneOrder, labels = c(new_gene_names$new_names),
                          name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))
    
    #And finally saving that plot.
    save_plot(plot_heatmap, 1500, 1500, paste(args[1],"/finalPlot.EMF", sep =""))
  }
}

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
# geom_range	bar layer to present uncertainty of evolutionary inFERence
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
