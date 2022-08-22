[Usage](#usage) 

[Steps in pipeline](#steps-in-pipeline) 

[Output](#output)

[Required dependencies](#required-dependencies-and-their-associated-dependencies) 

[Example plot](#example-plot-generated)


# PhyloTree - Pipeline
This is a pipeline that can be used to generate a phylogenetic tree of species of interest, including a heatmap showing carriage of specific genes, the assigned grouping of strains, and the names of each strain. Required for this to run is a list of genome sequences of interest, a list of gene sequences of interest, genomes that can be used to help expand the tree (as template genomes), and a list of genes that can split the species used into groups. Throughout this explanation, we'll refer to *Escherichia coli* (*E.coli*) as an example species with multiple strains that can be used to generate our Phylogenetic tree.


# Usage 
`bash PhyloTree.sh [OPTIONAL PARAMETERS]`

   **[--genes, -p]** : Path to location of genes of interest directory, which should contain the genes of interest. Default is set to your working directory, and the file "Genes". E.g. your_working_directory/Genes. Inside this file should be the Template_Genes file, which contains a FASTA format list of genes used to separate strains for the tree 

   **[--genomes, -h]** : Path to location of genomes of interest directory, which should contain the genomes of interest. Default is set to your working directory, and the file "Genomes". E.g. your_working_directory/Genomes. Inside this folder should be a file called Template_Genomes - containing a number of genomes to help tree building - these can be separate from the genomes of interest, which should just be stored in the Genomes file. 

   **[--output, -y]** : File name after this should indicate the path to the output file. Default is set to your working directory, in a file called "output". E.g. your_working_directory/output. 

   **[--rename_genomes, -l]** : Indicates whether user has "new_names.txt" file in working directory to rename old genomes into new names. Default is "no". If used, the file should contain tab separated old names and new names, line by line. 

   **[--phylogroup, -o]** : Indicates whether the user has used the PhyloGroup.sh script, which auto-generates a group_list.txt which is used to group *E. coli* strains ONLY. If this parameter is set to "yes", then the grouping parameter is automatically set to "yes", and an output folder is created which contains a list of genomes and a group_list.txt file. Default is "no". 
   
   **[--grouping, -t]** : This indicates whether the user has a group_list.txt file which contains a grouping factor tab separated from the name of the strain being used. Default is "no". The .txt file should be located in the working directory. 

# Steps in pipeline

This is the script where we want to examine the carriage of specific genes of interest over a number of E.coli (strains) genomes.
To do so, we generate a tree and use this to plot the carriage. The steps of this script are below:

   1. We need a fasta file of the genes of interest and independent files of the genomes of interest. You can put these in the Genes and Genomes directories respectively - not the template file though, those are to make the tree. Do not use symbol names in file names or in gene names inside the geneIO file, these might complicate things later.
	
   2. We run a blastn on each genome (with genes gathered by the ISME method used as queries) to determine which genes are present in the genome and their sequence in the respective genome. The results of this are stored in the BlastResults folder for each genome.
	
   3. We also require a list of genes of interest as an input file (concatenate the genes of interest into one single txt file, following a specific formatting i.e:
       
        \>geneName  
        ATCTAATTATATACATACATATAT etc.
   
      Stick to uncomplicated gene names (no symbols). On this, we will run a blastn to find genomes that have a match to the genes of interest in the list to identify carriage in the tree plot, and then save the list containing both the name of the genes and which genomes have matched with them. The result is a list of genes (per genome) that have been matched to the specific genome. This is saved as a specific output (see output section).
	
   4. We wish to perform alignment on the blastn results. I use part of Thomas' script for this, which is written in PERL, so I've tried to explain what goes on when we run it. Essentially we put all the results from the blastn (from step 3, step 4 blastn is finished and does not used get further except for the R plot) into one file per gene type, so if a gene is present in the genome file, it gets added to one file of that specific gene type, and then in the next genome if that gene is present the sequence is added to the same file - and so on for each gene in each genome. Essentially, you get a file for each gene and it contains the gene sequences from each genome that had a match. These are saved in the working directory, while they then undergo MAFFT alignment, which generates the files in the FINALTEST.Alignments folder (or whatever you decide to name the directory in the arguments for the perl script) which are, as the name suggests, the aligned gene sequences (of each genome match) for each gene. After this step, all the alignments for all the genes are added to the same file: Final.FINALTEST.aln. This generates a concatenated file of the alignments (note that alignment is before concatenation [of all the genes, not of individual genes though!] which is not exactly what Thomas recommended, but it still works!)
	
   5. This generates a concatenated file of the alignments, which is great but we aren't done yet. We need to make the tree with PhyML, however we cannot use a .aln file to generate this. So, we need to convert it to a phylip file format, using ALTER. This simply generates a .phy file, which I've designated phylipfor.phy (inventive, I know).
	
   6. Once we have a .phy file, we can generate a tree using phyml. This generates a whole bunch of files, but the important one is saved to the output folder (it's the one needed for R plots).
	
   7. Once that's done, we run the R script and generate a plot. Hooray! 

# Output
We generate a number of files in the output, they should be:
   
   **geneMatches.txt** (a tab separated file where one column is genes of interest and the other is genomes that had a match to those genes)
  
   **genesOI.txt** (just a column containing the genes of interest that were input by the user)
   
   **genomesAdded.txt** (a column containing the genomes of interest that were added by the user)
   
   **List.genomes.txt** (a column of the entire list of genomes used in the analysis - will be adjusted if the clades are removed or not)
   
   **phylipFor.phy_phyml_tree.txt** (the Newick file format containing the tree information, used in R to produce the tree)
   
   **phylogeny_list.txt** (tab-separated file, one column containing assigned phylogroup and the other being the matching genome)
   
   **newRscript.r** (the R script, but with the args[x] swapped with their actual parameter e.g. file path in the R script, allows user to use Rscript right after running with all data needed - they can adjust the R script as necessary)
   
   **The plot! As an EMF file. Shows heatmap indicating gene carriage for each strain, highlights genomes that were added in, distinct colours based on group and carriage.**  

If the user needs to change the names of some of the genomes (for instance, those that were added in) simply make a txt file called new_names.txt and have the old genome names on the left, tab, then the new names on the right. You can place that in the working directory. If required after making the original plot, simply take genomesAdded.txt from the output directory, add new names on the same line after tabbing, and run that part of the script in the newRscript.r output in the output directory.

# Required dependencies (and their associated dependencies)
**ALTER** https://github.com/sing-group/ALTER  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**JAVA** https://java.com/en/download/  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Maven Tool** https://maven.apache.org/download.cgi
   
**NCBI BLAST** https://www.ncbi.nlm.nih.gov/books/NBK569861/

**R** https://cran.rstudio.com/  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(possibly **RSTUDIO** as well) https://www.rstudio.com/products/rstudio/download/#download  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With other packages (install manually):  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**treeDataVerse** 
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**BiocManager** 
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**remotes**  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**devEMF** 
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**data.table**
      
**ClermonTyping** https://github.com/A-BN/ClermonTyping  
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**BioPython** from Python3 https://www.python.org/downloads/  
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**pandoc** https://pandoc.org/installing.html  
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;And other R packages     
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**readr**  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**dplyr**  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**tidyr**  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**stringr**  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**knitr**  
   
**PhyML** https://github.com/stephaneguindon/phyml  
&nbsp;&nbsp;&nbsp;If on Mac, you'll need to download the binary from the PhyML website, and move the UNZIPPED file to the same working directory as the rest of the script: http://www.atgc-montpellier.fr/phyml/download.php 
&nbsp;&nbsp;&nbsp;In this case, swap line 254 in PhyloTree.sh from `phyml -i phylipFor.phy -b 100` to `PhyML-3.1/PhyML-3.1_macOS-MountainLion -i phylipFor.phy -b 100`

**MAFFT** https://mafft.cbrc.jp/alignment/software/source.html

# Example plot generated
![finalPlot jpeg](https://user-images.githubusercontent.com/100131598/167879830-0587c396-07ad-456a-92fa-24dccf75653a.jpg)

