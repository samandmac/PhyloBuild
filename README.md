# Ecoli-Phylogentic-Tree-Pipeline
This is a pipeline that can be used to generate a phylogenetic tree of E.coli, including a heatmap showing carriage of specific genes, the assigned phylogroup of each strain, and the names of each strain.

# Usage 
`bash ecoli_phylogroup_tree.sh [OPTIONAL PARAMETERS]`

   **[--genes, -s]** : Path to location of genes of interest directory, which should contain the genes of interest. Default is set to your_working_directory/GeneInterest

   **[--genomes, -a]** : Path to location of genomes of interest directory, which should contain the genomes of interest. Default is set to your_working_directory/GenomeInterest

   **[--output, -m]** : File name after this should indicate the path to the output file. Default is set to your_working_directory/output.

   **[--clades, -u]** : Indicates whether clades {I,II,III,IV,V,VI} should be included in the plot. "yes" or "no" selections only. Default is "no".

   **[--label, -e]** : Indicates whether the labels should be included in the end plot. Default is set to "no".

# Steps in pipeline

This is the script where we want to examine the carriage of specific genes of interest over a number of E.coli (strains) genomes.
To do so, we generate a tree and use this to plot the carriage. The steps of this script are below:

   1. We need a fasta file of the genes of interest (we've taken this from 1 strain as a reference). For this, I initially used the Clermont scheme genes (i.e. looked at what genes the primers are supposed to select for, and took those). I switched that to a method similar to 		the ISME paper, which iteratively looks for genes present in the genomes and only uses those that are kept in the sequences. The intial genes are all those present in a reference strain. Therefore, the script initially determines the shared genes in all genome sequences being queried (and those that make up the reference tree!). For the moment, this is only set to the top 10 results (no metric, just the first 10 in the list), but this can be changed as necessary.
		
   2. We have to determine the phylogrouping of each genome that is being queried (and the references), so therefore the script also runs ClermonTyping and the results of which are saved (see output section).
	
   3. We run a blastn on each genome (with genes gathered by the ISME method used as queries) to determine which genes are present in the genome and their sequence in the respective genome. The results of this are stored in the BlastResults folder for each genome.
	
   4. We also require a list of genes of interest as an input file (concatenate the genes of interest into one single txt file, following a specfific formatting i.e:
       
        \>geneName  
        ATCTAATTATATACATACATATAT etc.
   
      Preferably, stick to uncomplicated gene names (no symbols). On this, we will run a blastn to find genomes that have a match to the genes of interest in the list to identify carriage in the tree plot, and then save the list containing both the name of the genes and which genomes have matched with them. The result is a list of genes (per genome) that have been matched to the specific genome. This is saved as a specific output (see output section).
	
   5. We wish to perform alignment on the blastn results. I use part of Thomas' script for this, which is written in PERL, so I've tried to explain what goes on when we run it. Essentially we put all the results from the blastn (from step 3, step 4 blastn is finished and does not used get further except for the R plot) into one file per gene type, so if a gene is present in the genome file, it gets added to one file of that specific gene type, and then in the next genome if that gene is present the sequence is added to the same file - and so on for each gene in each genome. Essentially, you get a file for each gene and it contains the gene sequences from each genome that had a match. These are saved in the working directory, while they then undergo MAFFT alignment, which generates the files in the FINALTEST.Alignments folder (or whatever you decide to name the directory in the arguments for the perl script) which are, as the name suggests, the aligned gene sequences (of each genome match) for each gene. After this step, all the alignments for all the genes are added to the same file: Final.FINALTEST.aln. This generates a concatenated file of the alignments (note that alignment is before concatenation [of all the genes, not of individual genes though!] which is not exactly what Thomas recommended, but it still works!)
	
   6. This generates a concatenated file of the alignments, which is great but we aren't done yet. We need to make the tree with PhyML, however we cannot use a .aln file to generate this. So, we need to convert it to a phylip file format, using ALTER. This simply generates a .phy file, which I've designated phylipfor.phy (inventive, I know).
	
   7. Once we have a .phy file, we can generate a tree using phyml. This generates a whole bunch of files, but the important one is saved to the output folder (it's the one needed for R plots).
	
   8. Once that's done, we run the R script and generate a plot. Hooray! 

# Output
We generate a number of files in the output, they should be:
   
   **geneMatches.txt** (a tab separated file where one column is genes of interest and the other is genomes that had a match to those genes)
  
   **genesOI.txt** (just a column containing the genes of interest that were input by the user)
   
   **genomesAdded.txt** (a column containing the genomes of interest that were added by the user)
   
   **List.genomes.txt** (a column of the entire list of genomes used in the analysis - will be adjusted if the clades are removed or not)
   
   **phylipFor.phy_phyml_tree.txt** (the Newick file format containing the tree information, used in R to produce the tree)
   
   **phylogeny_list.txt** (tab-separated file, one column containing assigned phylogroup and the other being the matching genome)
   
   **newRscript.r** (the R script, but with the args[x] swapped with their actual parameter e.g. file path in the R script, allows user to peruse it at their leisure, and change things if they wish)
   
   **The plot! As an EMF file. Shows heatmap indicating gene carriage for each strain, highlights genomes that were added in, distinct colours based on group and carriage.**
 
# Required dependencies (and their associated dependencies)
**ALTER** https://github.com/sing-group/ALTER  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**JAVA** https://java.com/en/download/  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Maven Tool** https://maven.apache.org/download.cgi
   
**NCBI BLAST** https://www.ncbi.nlm.nih.gov/books/NBK569861/

**R** https://cran.rstudio.com/  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(possibly **RSTUDIO** as well) https://www.rstudio.com/products/rstudio/download/#download  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;With other packages:  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**treeDataVerse** (installed in script)  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**BiocManager** (installed in script)  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**remotes** (installed in script)  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**devEMF** (installed in script)  
      
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
&nbsp;&nbsp;&nbsp;In this case, swap line 254 in ecoli_phylogroup_tree.sh from `phyml -i phylipFor.phy -b 100` to `PhyML-3.1/PhyML-3.1_macOS-MountainLion -i phylipFor.phy -b 100`

**MAFFT** https://mafft.cbrc.jp/alignment/software/source.html

# Example plot generated
![finalPlot jpeg](https://user-images.githubusercontent.com/100131598/167879830-0587c396-07ad-456a-92fa-24dccf75653a.jpg)

