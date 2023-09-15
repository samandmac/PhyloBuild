[Usage](#usage) 

[Steps in pipeline](#steps-in-pipeline) 

[Output](#output)

[Required dependencies](#required-dependencies-and-their-associated-dependencies) 

[Example plot](#example-plot-generated)


# PhyloBuild - Pipeline
Included in this software are a number of scripts that can be used to generate a phylogenetic tree of particular genomes/strains, investigate carriage of specific genes of interest in genomes/strains, or use the software to visualise a tree previously obtained through PhyloBuild or other means. A number of different parameters allow for user customisation of options for their run, and with visualisation. These are further explained below.

# Prior to PhyloBuild.sh
Some tools may be used before PhyloBuild.sh is ran, which will be described here.  

   **PhyloGenes.sh** allows user to auto-generate a "rough-and-ready" geneList.txt file (which is used to split the genomes/strains on the tree by their evolutionary similarity) containing orthologues present in all genomes/strains used. *The user should consider generating their own list of orthologues, or regions of similarity, instead of relying on this method. However, it exists for now as a quick way to run PhyloBuild.*  By inputting a reference genome full set of gene sequences. This can easily be downloaded from NCBI, although a caveat of this step is that the more genomes used the better (around 131 different strains in our tests gives off good separation in the tree) - as some orthologs may have little variation in fewer genomes/strains. Alternatively, the user can collect their own orthologs or genes known to be present in a subset of the species to split the genomes/strains by their evolutionary similarity.
   
   **PhyloGroup.sh** allows the user to group E.coli strains by their phylogroup by generating a group_list.txt file (this is simply used to colour genome/strain names on the final plot, and is run optionally), designated through the Clermont Scheme sequences. Although a lot of testing has been performed utilising E.coli datasets, the software is not hard-coded to one particular species, and this should be considered an extra feature. User can indicate their own grouping scheme for strains simply by making their own group_list.txt file (tab separated, grouping then the strain name).

These can be ran simply as `bash Phyloxxxx.sh`, which makes the assumption you are working from the working directory containing the files as in the gitHub default. The scripts do have default parameter options, so do ensure that you make sure you are selecting the parameters specific to what you want to happen in the run.

# Usage 
`bash PhyloBuild.sh [OPTIONAL PARAMETERS]`

   **[--genes_interest, -a]** : Path to location of genes of interest directory, which should contain the genes of interest. Default is set to your working directory, and the file "Genes". E.g. your_working_directory/Genes. Inside this file should be the Template_Genes file, which contains a FASTA format list of genes used to separate strains for the tree 

   **[--genomes_interest, -b]** : Path to location of genomes of interest directory, which should contain the genomes of interest. Default is set to your working directory, and the file "Genomes". E.g. your_working_directory/Genomes. Inside this folder should be a file called Template_Genomes - containing a number of genomes to help tree building - these can be separate from the genomes of interest, which should just be stored in the Genomes file. 

   **[--output, -c]** : File name after this should indicate the path to the output file. Default is set to your working directory, in a file called "output". E.g. your_working_directory/output. 
 
   **[--tree_genes, -d]** : Indicates user directory for genes used to help build the tree. Default is set to working directory, in a filed called Tree_Genes.

   **[--tree_genomes, -e]** : Indicates user directory for genomes used to help build the tree. Default is set to working directory, in a file called Tree_Genomes.
   
   **[--grouping, -f]** : This indicates whether the user has a "group_list.txt" file in the working directory which contains a grouping factor tab separated from the name of the strain being used. Default is "no". If using PhyloGroup first, a group_list.txt is automatically generated.  
   
   **[--rename_genomes, -g]** : Indicates whether user has "new_genome_names.txt" file in working directory to rename old genomes into new names. Default is "no". If used, the file should contain tab separated old names and new names, line by line. 

   **[--phylogroup, -h]** : Indicates whether the user has used the PhyloGroup.sh script, which auto-generates a group_list.txt which is used to group *E. coli* strains ONLY. If this parameter is set to "yes", then the grouping parameter is automatically set to "yes", and an output folder is created which contains a list of genomes and a group_list.txt file. Default is "no".
   
   **[--mac, -i]** : This indicates whether the user has a mac or not (used for PHYML, check below to see steps for downloading PHYML on MAC). Default is "no".  
   
   **[--rename_genes, -j]** : This indicates whether the user has a file called "new_gene_names.txt" in their working directory to rename old gene names into new gene names on the plot. Default is set to "no". Should be a tab separated file, containing the old name, then the new name.  
   
   **[--phylogenes, -k]** : Indicates whether the user previously utilised the PhyloGenes.sh command, which auto-generates a "geneList.txt" file in the Tree_Genes folder. Default is set to "no". See instructions for tool in geneList1.txt file in PhyloGenes directory, and above.  
   
   **[--genomes_oi, -l]** : Indicates whether the user has genomes of interest, these should be located in the Genomes_Interest folder and will be run alongside the other genomes that are not considered of interest but are used to help with populating the tree. Default is set to "no".  
   
   **[--genes_oi, -m]** : Indicates whether the user has genes of interest, these should be located in the Carriage_Genes folder and so the genomes used will be subjected to a blastn to determine which genes are carried in which genomes - this will generate a geneMatches.txt file int he output directory. Default is set to "no".  
   
   **[--sample_sheet, -n]** : Indicates whether the user has already generated a sample sheet that contains information pertaining to the visualisation step. Should this be set to "no" and the tree is generated, the sample sheet will be automatically generated to indicate parameters chosen by the user. Default is set to "no".  
   
   **[--run_phyloplot, -o]** : Indicates whether the user wants to run the visualisation step after generating a tree or not. Default is set to "no".  

   **[--viral_strains, -p]** : Indicates whether the user is using viral strains or not, if they are then the parameters for the blastn are slightly altered due to the smaller size in viral genomes. Default is set to "no".
   
   **[--make_tree, -q]** : Indicates whether the user wants to actually make a tree or not, if this is set to "no" and genes_oi is set to "yes" then the process will continue up to checking for gene carriage and then quit. Default is set to "yes".  
   
# Steps in pipeline

We could use PhyloBuild.sh to generate a tree and use this to plot the carriage. The typical steps of this script are below, but of course depend on what parameters the user has chosen:

   1. We need a fasta file containing sequences of genes that will be used to separate the tree by evolutionary similarity, and we need a number of different fasta files of genomes/strains that we can use for the tree. Optionally, we could also include a fasta file containing the gene sequences for investigating carriage, and also populate the Genomes_Interest folder with genomes/strains that we are specifically interested in highlighting.
	
   2. We run a blastn on each genome using the genes to split the tree, and optionally to determine which genes are carried in the genomes/strains, and extract their sequence in each respective genome. The results of this are stored in the BlastResults folder for each genome, which is deleted after the run is over.
	
   3. We wish to perform alignment on the blastn results for each extracted sequence. Essentially we put all the results from the blastn into one file per gene type, so if a gene is present in the genome file, it gets added to one file of that specific gene type, and then in the next genome if that gene is present the sequence is added to the same file - and so on for each gene in each genome. Essentially, you get a file for each gene and it contains the gene sequences from each genome that had a match. These are saved in the working directory, while they then undergo MAFFT alignment, which generates the files in the FINALTEST.Alignments folder (or whatever you decide to name the directory in the arguments for the perl script) which are, as the name suggests, the aligned gene sequences (of each genome match) for each gene. After this step, all the alignments for all the genes are added to the same file: Final.FINALTEST.aln. This generates a concatenated file of the alignments.
	
   5. This generates a concatenated file of the alignments, which is great but we aren't done yet. We need to make the tree with PhyML, however we cannot use a .aln file to generate this. So, we need to convert it to a phylip file format, using ALTER. This simply generates a .phy file, which I've designated phylipfor.phy.
	
   6. Once we have a .phy file, we can generate a tree using phyml. This generates a whole bunch of files, but the important one is saved to the output folder (it's the one needed for R plots).
	
   7. Once that's done, we can run the R script and generate a plot. 

# Output
We generate a number of files in the output, they could be:
   
   **geneMatches.txt** (a tab separated file where one column is genes of interest and the other is genomes that had a match to those genes)
  
   **genesOI.txt** (just a column containing the genes of interest that were input by the user)
   
   **genomesAdded.txt** (a column containing the genomes of interest that were added by the user)
   
   **ListGenomes.txt** (a column of the entire list of genomes used in the analysis - will be adjusted if the clades are removed or not)
   
   **TreeDataFile.txt** (the Newick file format containing the tree information, used in R to produce the tree)
   
   **phylogeny_list.txt** (tab-separated file, one column containing assigned phylogroup and the other being the matching genome)
   
   **newRscript.r** (the R script, but with the args[x] swapped with their actual parameter e.g. file path in the R script, allows user to use Rscript right after running with all data needed - they can adjust the R script as necessary)
   
   **The plot! As an EMF file. Shows heatmap indicating gene carriage for each strain, highlights genomes that were added in, distinct colours based on group and carriage.**  
   
   **A heatmap showing carriage only, in the strains that you specified as of interest (so not every single strain, unless they are located in the Genomes_Interest folder and you specified that this was a parameter), not including a tree. As an EMF file.**  
   
   **A tab-separated .txt file called sample_sheet.txt in your working directory, that will contain the information designated by the parameters specified by the user for the visualisation step - this is automatically generated once PhyloBuild.sh is run, unless the user has created their own sample sheet and specified this in the parameters.**  

If the user needs to change the names of some of the (carriage) genes or genomes, simply make a txt file called new_gene_names.txt or new_genome_names.txt and have the old genome names on the left, tab, then the new names on the right. You can place that in the working directory. If required after making the original plot, simply take genomesAdded.txt from the output directory, add new names on the same line after tabbing, and run that part of the script in the newRscript.r output in the output directory. 

# Required dependencies (and their associated dependencies)
**trimAl** https://github.com/sing-group/ALTER  
   
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
&nbsp;&nbsp;&nbsp;If using MAC, ensure "--mac yes" is used in user parameters.

**MAFFT** https://mafft.cbrc.jp/alignment/software/source.html

# Example plots generated
![finalPlot](https://user-images.githubusercontent.com/100131598/199327999-0e6ba88a-a502-4ada-ac83-a64e69d39491.jpg)
![heatmap](https://user-images.githubusercontent.com/100131598/228218314-bca69119-0426-4ea8-9d37-bf3a9a0acb37.jpg)
