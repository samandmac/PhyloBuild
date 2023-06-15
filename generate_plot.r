#!/usr/bin/env Rscript

#Because we added in the user input options, we need to change file locations to
#the user choices by adding in parameters.
args = commandArgs(trailingOnly = TRUE) #We just run Rscript with a path to the output files which contains files for plots.

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

#The treeDataFile is in a newick file format, so we can easily load it in as below.
tree = read.newick(paste(args[1],"/treeDataFile.txt", sep = ""),
                   node.label = "support")

#This just gives some information about the tree
get.tree(tree)

#Converting to tibble for adding in phylogroup data - easier than adjusting tree object.
tree.tbl = as_tibble(tree)

#Removing .blastn at end of name so that the phylogroup and tree use same label names
tree.tbl$label = gsub(".blastn", "", tree.tbl$label)

sample_sheet = read.table("sample_sheet.txt", header = T)

groups_yn = args[3]

if (groups_yn == "yes"){
  #Loading phylogrouping from file, the gene matches, gene names, and genomes added
  phylogroup = sample_sheet[,c(1,3)]
  phylogroup$group = as.factor(phylogroup$group) #Changing the grouping to be a factor
  names(phylogroup) = c("label", "group")
}

gene_names = args[4]

if (gene_names == "yes"){
  new_gene_names = read.table(paste(args[1],"/../new_gene_names.txt", sep = ""), sep = "\t")
  names(new_gene_names) = c("old_names","new_names")
}

genes_oi = args[5]

if (genes_oi == "yes"){
  geneMatches=read.table(paste(args[1],"/geneMatches.txt", sep=""),
                         col.names = c("gene","label"))
  genes = read.table(paste(args[1],"/genesOI.txt", sep=""))
  geneOrder = as.vector(unlist(genes))
  geneMatches$gene = as.factor(geneMatches$gene)
  geneMatches$gene = factor(geneMatches$gene, levels = geneOrder)
}

#levels(geneMatches$gene)
#levels(phylogroup$group)


genomes_oi = args[6]
if (genomes_oi == "yes"){ 
  genomes_added_in = vector()
  
  for (x in 1:length(row.names(sample_sheet))){
    if (sample_sheet$of_interest[x] == "yes"){
      genomes_added_in = append(genomes_added_in, x)
    }
  }
  
  genomes_added_in = sample_sheet$strain[genomes_added_in]
}

renaming_genomes = args[2]

#This part just finds if the user replied yes or no and gives a difFERent phylogroup list which is used later in a loop.
if (renaming_genomes == "yes") {
  #This little bit of script can be used to change names if required. If you do, place the "new_names" file into the output folder.
  new_genome_names = sample_sheet[,c(1,2)]
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
  phylogroup1 = phylogroup[2]
  row.names(phylogroup1) = phylogroup[,1]
}

#HEATMAP FOR GENE CARRIAGE
#Below makes the dataset for the heatmap for the tree plot

if (genes_oi == "yes"){
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
    heatmap = heatmap[,-1)] #Remove the genome name as a column.
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
}


#Let's also generate a heatmap of just the gene carriage alongside the tree plot.
#We can do one for the genomes of interest, and one for all genes.
if (genes_oi == "yes"){
  if (gene_names == "no"){
    if(renaming_genomes == "no"){
      gene_and_genome_only = subset(geneMatches, label %in% genomes_added_in)
      map = ggplot(gene_and_genome_only, aes(x = gene, y = label, fill = gene))+
        geom_tile(colour="white", size=0.25)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none", plot.background=element_blank(), panel.border=element_blank())+
        #ylab(label = "Strain")+
        #xlab(label = "Gene")
        labs(x="", y="")
    }
    else {
      gene_and_genome_only = subset(geneMatches, label %in% genomes_added_in)
      
      gene_and_genome_only$gene = as.character(gene_and_genome_only$gene)
      gene_and_genome_only$label = as.character(gene_and_genome_only$label)
      
      
      for (x in 1:length(new_genome_names$label)){
        replace = grep(new_genome_names$label[x],gene_and_genome_only$label)
        gene_and_genome_only$label[replace] = new_genome_names$new_names[x]
      }
      
      gene_and_genome_only$gene = as.factor(gene_and_genome_only$gene)
      gene_and_genome_only$label = as.factor(gene_and_genome_only$label)
      
      map = ggplot(gene_and_genome_only, aes(x = gene, y = label, fill = gene))+
        geom_tile(colour="white", size=0.25)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none", plot.background=element_blank(), panel.border=element_blank())+
        labs(x="", y="")
        #ylab(label = "Strain")+
        #xlab(label = "Gene")
    }
    save_plot(map, 500, 500, paste(args[1],"/heatmap_carriage.emf", sep =""))
  }
  
  if (gene_names == "yes"){
    if (renaming_genomes == "no"){
      gene_and_genome_only = subset(geneMatches, label %in% genomes_added_in)
      gene_and_genome_only$gene = as.character(gene_and_genome_only$gene)
      gene_and_genome_only$label = as.character(gene_and_genome_only$label)
      
      
      for (x in 1:length(new_gene_names$old_names)){
        replace = grep(new_gene_names$old_names[x],gene_and_genome_only$gene)
        gene_and_genome_only$gene[replace] = new_gene_names$new_names[x]
      }
      
      gene_and_genome_only$gene = as.factor(gene_and_genome_only$gene)
      gene_and_genome_only$label = as.factor(gene_and_genome_only$label)
      
      map = ggplot(gene_and_genome_only, aes(x = gene, y = label, fill = gene))+
        geom_tile(colour="white", size=0.25)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.background=element_blank(), panel.border=element_blank())+
        scale_fill_discrete(breaks=geneOrder, labels = c(new_gene_names$new_names),
                            name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))+
        #ylab(label = "Strain")+
        #xlab(label = "Gene")+
        #labs(fill = "Genes")
        labs(x="", y="")
        }
    else{
      gene_and_genome_only = subset(geneMatches, label %in% genomes_added_in)
      gene_and_genome_only$gene = as.character(gene_and_genome_only$gene)
      gene_and_genome_only$label = as.character(gene_and_genome_only$label)
      
      
      for (x in 1:length(new_gene_names$old_names)){
        replace = grep(new_gene_names$old_names[x],gene_and_genome_only$gene)
        gene_and_genome_only$gene[replace] = new_gene_names$new_names[x]
      }
      
      for (x in 1:length(new_genome_names$label)){
        replace = grep(new_genome_names$label[x],gene_and_genome_only$label)
        gene_and_genome_only$label[replace] = new_genome_names$new_names[x]
      }
      
      gene_and_genome_only$gene = as.factor(gene_and_genome_only$gene)
      gene_and_genome_only$label = as.factor(gene_and_genome_only$label)
      
      map = ggplot(gene_and_genome_only, aes(x = gene, y = label, fill = gene))+
        geom_tile(colour="white", size=0.25)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.background=element_blank(), panel.border=element_blank())+
        scale_fill_discrete(breaks=geneOrder, labels = c(new_gene_names$new_names),
                            name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))+
        #ylab(label = "Strain")+
        #xlab(label = "Gene")+
        #labs(fill = "Genes")
        labs(x="", y="")
    }
    save_plot(map, 500, 500, paste(args[1],"/heatmap_carriage.emf", sep =""))
  }
}

if (groups_yn == "yes"){
  if (gene_names == "no"){
    #Below is the script for when pgroup labels aren't required.
    #Back to making the actual plot now, and saving that.
    #Added in colour=group to ggtree to change the branch colour to group 
    tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, aes(colour=group), show.legend=F) %<+% d + #show.legend =F to get rid of line in legend
      #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_added_in))+ #Changed geom_tiplab to be only the genomes added in
      #Alternative plot below, highlight branches of genomes added in, and give all the genomes names as tip labels?
      geom_tiplab(size=3, aes(colour=group, label=label2), offset=15)+ #changed offset to 28 from 26 for Clermont genes
      scale_color_discrete("Group")+ #This is just another way of getting the legend in order
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    if (genes_oi == "yes") {
      #This adds the heatmap onto the plot
      tree.plot = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
        scale_fill_discrete(breaks=geneOrder, 
                            name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))+
        theme(legend.box.spacing = unit(0, "pt"))
    }
    #And finally saving that plot.
    save_plot(tree.plot, 1300, 1300, paste(args[1],"/finalPlot.EMF", sep =""))
  }
  
  if (gene_names == "yes"){
    #Below is the script for when pgroup labels aren't required.
    #Back to making the actual plot now, and saving that.
    #Added in colour=group to ggtree to change the branch colour to group 
    tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, aes(colour=group), show.legend=F) %<+% d + #show.legend =F to get rid of line in legend
      #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_added_in))+ #Changed geom_tiplab to be only the genomes added in
      #Alternative plot below, highlight branches of genomes added in, and give all the genomes names as tip labels?
      geom_tiplab(size=3, aes(colour=group, label=label2), offset=15)+ #changed offset to 28 from 26 for Clermont genes
      scale_color_discrete("Group")+ #This is just another way of getting the legend in order
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    if (genes_oi == "yes"){
      #This adds the heatmap onto the plot
      tree.plot = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
        scale_fill_discrete(breaks=geneOrder, labels = c(new_gene_names$new_names),
                            name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))+
        theme(legend.box.spacing = unit(0, "pt"))
    }
      #And finally saving that plot.
    save_plot(tree.plot, 1300, 1300, paste(args[1],"/finalPlot.EMF", sep =""))
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
      geom_tiplab(size=3, aes(label=label2), offset=15)+ #changed offset to 28 from 26 for Clermont genes
      scale_color_discrete("Group")+
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    if (genes_oi == "yes"){
      #This adds the heatmap onto the plot
      tree.plot = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
        scale_fill_discrete(breaks=geneOrder, 
                            name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))+
        theme(legend.box.spacing = unit(0, "pt"))
    }
    
    #And finally saving that plot.
    save_plot(tree.plot, 1300, 1300, paste(args[1],"/finalPlot.EMF", sep =""))
  }
  if (gene_names == "yes"){
    #Below is the script for when pgroup labels aren't required.
    #Back to making the actual plot now, and saving that.
    #Added in colour=group to ggtree to change the branch colour to group 
    tree.plot = ggtree(tree, layout='circular',branch.length = 'none', ladderize = FALSE, show.legend=F) %<+% d + #show.legend =F to get rid of line in legend
      #geom_tiplab(size=3, aes(colour=group), offset=23, data = td_filter(isTip & node %in% nodes_added_in))+ #Changed geom_tiplab to be only the genomes added in
      #Alternative plot below, highlight branches of genomes added in, and give all the genomes names as tip labels?
      geom_tiplab(size=3, aes(label=label2), offset=15)+ #changed offset to 28 from 26 for Clermont genes
      scale_color_discrete("Group")+
      guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 3)))
    
    if (genes_oi == "yes"){
      #This adds the heatmap onto the plot
      tree.plot = gheatmap(tree.plot, heatmap, offset=0, width=0.5, font.size=3, colnames = FALSE)+
        scale_fill_discrete(breaks=geneOrder, labels = c(new_gene_names$new_names),
                            name="Gene Carriage", guide = guide_legend(label.theme = element_text(face = "italic")))+
        theme(legend.box.spacing = unit(0, "pt"))
    }
    #And finally saving that plot.
    save_plot(tree.plot, 1300, 1300, paste(args[1],"/finalPlot.EMF", sep =""))
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
        #genomes_oi = args[6]
        #if (genomes_oi == "yes"){
          #This allows us to select the nodes for the genomes added in, for altering the plot!
          #merge_of_genomes_added_in = merge(tree.tbl, genomes_added_in, by.x = 4, by.y = 1) #We previously made genomes added in, so we can just merge them by their label name
        #nodes_added_in =unlist(merge_of_genomes_added_in[,3]) #Saving the nodes of the genomes_added_in
        #nodes_added_in
      #}
      #And layer for plot
      # geom_hilight(mapping=aes(subset = node %in% nodes_added_in, fill = "red"))+ #Highlights nodes of genomes that were added in
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

#To remove strains in your tree e.g.
# to_drop = c("R20291_reference","FM2.5_PaLoc")
# 
# tree = drop.tip(tree, to_drop)
# 
# tree2 = root(tree, outgroup = "R20291_denovo", edgelabel = T)
# 
# tree_plot = ggtree(tree2, branch.length = "none") +
#   geom_tiplab()+
#   ggplot2::xlim(0, 4)
# 
# save_plot(tree_plot, 500,500,"/media/sam/Expansion/tree_plot_no_PALoc_no_ref.emf")
