# Ecoli-Phylogentic-Tree-Pipeline

This is a pipeline that can be used to generate a phylogenetic tree of E.coli, including a heatmap showing carriage of specific genes, the assigned phylogroup of each strain, and the names of each strain.

The user can download the zip file, which contains the script and necessary folders to run it.

Some other packages were used in the process, including:  
&nbsp;&nbsp;ALTER (requires java)  
&emsp;&emsp;&emsp;&emsp; citation - https://pubmed.ncbi.nlm.nih.gov/20439312/  
&nbsp;&nbsp;phyML (can be downloaded as command line tool)  
&emsp;&emsp;&emsp;&emsp;citation - https://academic.oup.com/sysbio/article/59/3/307/1702850  
&nbsp;&nbsp;R (requires treedataverse package - may be installed in the script if the user does not already have it installed)  
&emsp;&emsp;&emsp;&emsp;citation - https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.96  
&nbsp;&nbsp;BLAST (requires download of blast package)  
&emsp;&emsp;&emsp;&emsp;citation - https://pubmed.ncbi.nlm.nih.gov/2231712/  
&nbsp;&nbsp;ClermonTyping (requires download of the ClermonTyping command line tool)  
&emsp;&emsp;&emsp;&emsp;citation - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6113867/  
&nbsp;&nbsp;mafft alignment (requires download of mafft aligner)  
&emsp;&emsp;&emsp;&emsp;citation - https://academic.oup.com/bib/article/20/4/1160/4106928?login=false  
  
BE AWARE: The tools presented above are likely to have further dependencies needed to be downloaded for use.

Thank you to Thomas Otto for providing perl script to run mafft alignment and parse alignments into suitable format for ALTER.

The steps of this script are below:

&emsp;&emsp;1. We need a fasta file of the genes of interest (we've taken this from 1 strain as a reference). For this, I initially used the Clermont scheme genes (i.e. looked at what genes the primers are supposed to select for, and took those). I switched that to a method similar to an ISME paper (https://www.nature.com/articles/ismej2014242), which iteratively looks for genes present in the genomes and only uses those that are kept in the sequences. The intial genes are all those present in a reference strain. Therefore, the script initially determines the shared genes in all genome sequences being queried (and those that make up the reference tree!). For the moment, this is only set to the first 10 genes (for no reason other than my computer runs very slowly, not enough cores!).
		
&emsp;&emsp;2. We have to determine the phylogrouping of each genome that is being queried (and the references), so therefore the script also runs ClermonTyping and the results of which are saved (see output section).
	
&emsp;&emsp;3. We run a blastn on each genome (with genes gathered by the ISME method used as queries) to determine which genes are present in the genome and their sequence in the respective genome. The results of this are stored in the BlastResults folder for each genome.
	
&emsp;&emsp;4. We also require a list of genes of interest as an input file (concatenate the genes of interest into one single txt file, following a specfific formatting i.e:  
	&emsp;&emsp;&emsp;&emsp;	>geneName  
	&emsp;&emsp;&emsp;&emsp;	ATCTAATTATATACATACATATAT etc.  
Preferably, stick to uncomplicated gene names (no symbols - something like EC0_06_19646_7.fasta is perfect, but the script can handle altering EC0_06_19646_7#6.contigs_velvet.fa and converting it to EC0_06_19646_7.fasta. NOTE: outside of this format however, the script would need to be adjusted). On this, we will run a blastn to find genomes that have a match to the genes of interest in the list to identify carriage in the tree plot, and then save the list containing both the name of the genes and which genomes have matched with them. The result is a list of genes (per genome) that have been matched to the specific genome. This is saved as a specific output (see output section).
	
&emsp;&emsp;5. We wish to perform alignment on the blastn results. I use part of a script provided by Thomas Otto for this, which is written in PERL, so I've tried to explain what goes on when we run it. Essentially we put all the results from the blastn (from step 3, step 4 blastn is finished and does not used get further except for the R plot) into one file per gene type, so if a gene is present in the genome file, it gets added to one file of that specific gene type, and then in the next genome if that gene is present the sequence is added to the same file - and so on for each gene in each genome. Essentially, you get a file for each gene and it contains the gene sequences identified from each genome that had a match. These are saved in the working directory. They then undergo MAFFT alignment, which generates the files in the ALIGNMENTS.Alignments folder (or whatever you decide to name the directory in the arguments for the perl script) which are, as the name suggests, the aligned gene sequences (of each genome match) for each gene. After this step, all the alignments for all the genes are added to the same file: Final.ALIGNMENTS.aln. This generates a concatenated file of the alignments under each genome (now separated by genome rather than genetype) (note that alignment is before concatenation [of all the genes, not of individual genes though!] which is not exactly what Thomas recommended, but it still works!)
	
&emsp;&emsp;6. This generates a concatenated file of the alignments, which is great but we aren't done yet. We need to make the tree with PhyML, however we cannot use a .aln file to generate this. So, we need to convert it to a phylip file format, using ALTER. This simply generates a .phy file, which I've designated phylipfor.phy (phylip format . phy).
	
&emsp;&emsp;7. Once we have a .phy file, we can generate a tree using phyml. This generates a whole bunch of files, but the important one is saved to the output folder (it's the one needed for R plots).
	
&emsp;&emsp;8. Once that's done, we run the R script and generate a plot including a heatmap of gene carriage and phylogrouping. Red highlights indicate genome strains that were added in by the user. 

In addition to the above, the script can be run with some user input - the script can be run with --genomes, --genes, --output, or -s, -a, -m (respectively) to indicate the path to folders containing the genomes of interest, the genes of interest, and the output file (which can be made by the script or an existing file). 
The script will also ask if you wish to keep strains assigned clade phylogrouping (user's choice, apparently it is indicative of strains that are not E.coli - so if you are 100% sure you know what you are working with feel free to type "y" when prompted).

![plot](https://user-images.githubusercontent.com/100131598/163408185-3d2a6ead-d6b0-4d1a-ac76-5c455fefc67f.jpg)
