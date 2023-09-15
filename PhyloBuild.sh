#!/bin/bash

#This is the script where we can plot a rudimentary phylogenetic tree given some simple inputs of strain sequences, gene sequences, and a few other options.
#To do so, we generate a tree and use this to plot the carriage. The steps of this script are below:

#1. Format check of input files to make sure they're in FASTA format. Will return an error if not. Also some general parsing of files, creation of directories, and movement of files.
	
#2. We run a blastn on each strain (with genes either determined by the user or by the PhyloGenes method) to determine which genes are present in each strain and extract these from the strain. Each gene is tested, and the results of the blast are stored in the BlastResults folder for each genome - this will be cleared at the end of the script.
	
#3. Optionally, the user can input a list of genes of interest (concatenate the genes of interest into one single txt file, following FASTA format). Preferably, stick to uncomplicated gene names (no symbols). On this, we will run a blastn to find genomes that have a match to the genes of interest in the list to identify carriage in the tree plot, and then save the list containing both the name of the genes and which strains have matched with them. The result is a list of genes (per strain) that have been matched to the specific strain. They would do so to generate the heatmap around the tree, as seen in examples.
	
#4. We wish to perform alignment on the blastn results so that we can differentiate between similar or dissimiliar strains by their sequence matches from step 2. Essentially, we put all the results from the blastn (step 2 results) into one file per gene type, so if a gene is present in the strain, it gets added to a corresponding file labelled by gene type, and then in the next strain if that gene is present the sequence is added to the same file - and so on for each gene in each strain. Essentially, a FASTA file for each gene is created that contains the sequence matches from each strain matched the gene, separated by the original strain name as a FASTA header. These files then undergo MAFFT alignment, which generates the files in the Alignments folder which are, as the name suggests, the MAFFT aligned file generated from each of the strains matching sequences to a particular gene. Should a strain be missing a match, then this is reflected by adding in empty dashes matching the length of the aligned sequence of a particular gene. After this step, the sequences belonging to each respective strain are re-extracted and concatenated back together. For each strain, the concatenated alignments are then treated as the sequence belonging to the strain (in FASTA format). This generates a concatenated file of the alignments.
	
#5. This generates a concatenated file of the alignments, which is great but we aren't done yet. We need to make the tree with PhyML, however we cannot use an .aln file to generate this. So, we need to convert it to a phylip file format, using ALTER. This simply generates a .phy file, which I've designated phylipfor.phy (inventive, I know).
	
#6. Once we have a .phy file, we can generate a tree using phyml. This generates a whole bunch of files, but the important one is saved to the output folder (it's the one needed for R plots).
	
#7. Once that's done, we can then run the R script and generate a plot. This will be made in the output folder, alongside other data files used in the R script - including the list of gene carriage for the strains, a list of the genes of interest, a list of the genomes used, a list of the genomes that were added to the template folder (and removed), and the Newick format file which contains the tree information. I've since separated the tree making step from the visualisation of the plot, so they will need to be run separately - PhyloPlots.sh should be run after the PhyloBuild.sh (and can even be run straight after through the parameters).

#### START OF THE SCRIPT ####
#Here we begin timing to see how long processes are taking.
start=$SECONDS

#It's set up to take some arguements, -p (or --genes) for genes of interest path, -h (or --genomes) for genomes of interest path, -y (or --output) for the output path, -l (or --rename_groups) for indicating whether the working directory contains a file called new_names.txt which contains the old names of strains and then new names of strains, for name replacement, -o (or --phylogroup) which is a yes or no setting to indicate whether you used PhyloGroup beforehand, and -t (or --grouping) which indicates whether you included group_list.txt file to provide groups for the strains included in the run. You can run without parameters and by default it will use some defaults - working_directory/Genes, working_directory/Genomes, working_directory/output, "no", "no", and "no", respectively.

for arg in "$@"; do #This just sets it so that you can use the double --xxxx and it will convert those to -x.
  shift
  case "$arg" in
    "--genes_interest") set -- "$@" "-a" ;;
    "--genomes_interest") set -- "$@" "-b" ;;
    "--output")   set -- "$@" "-c" ;;
    "--tree_genes")   set -- "$@" "-d" ;;
    "--tree_genomes")   set -- "$@" "-e" ;;
    "--grouping")   set -- "$@" "-f" ;;
    "--rename_genomes")   set -- "$@" "-g" ;;
    "--phylogroup")   set -- "$@" "-h" ;;
    "--mac_yn")   set -- "$@" "-i" ;;
    "--rename_genes")   set -- "$@" "-j" ;;
    "--phylogenes")   set -- "$@" "-k" ;;
    "--genomes_oi")   set -- "$@" "-l" ;;
    "--genes_oi")   set -- "$@" "-m" ;;
    "--sample_sheet")   set -- "$@" "-n" ;;
    "--run_phyloplot")   set -- "$@" "-o" ;;
    "--viral_strains")   set -- "$@" "-p" ;;
    "--make_tree")   set -- "$@" "-q" ;;
    *)        set -- "$@" "$arg"
  esac
done

#This simply takes the arguement from the user and stores that in a variable.
while getopts "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:" opt
do
	case "$opt" in
		a)
			gene_dir="${OPTARG}" #User can decide gene directory
			;;
		b)
			genome_dir="${OPTARG}" #User can decide genome directory
			;;
		c)
			output_dir="${OPTARG}" #User can decide output directory
			;;
		d)
			tree_genes="${OPTARG}" 
			;;
		e)
			tree_genomes="${OPTARG}" 
			;;
		f)
			grouping_yn="${OPTARG}" 
			;;
		g)
			rename_genomes="${OPTARG}" #added this in too! Allows user to rename genomes they added in.
			;;
		h)
			phylogroup_yn="${OPTARG}"
			;;
		i)
			mac_yn="${OPTARG}" 
			;;
		j)
			rename_genes="${OPTARG}" 
			;;
		k)
			phylogenes="${OPTARG}" 
			;;
		l)
			genomes_oi="${OPTARG}" 
			;;
		m)
			genes_oi="${OPTARG}" 
			;;
		n)
			sample_sheet="${OPTARG}" 
			;;
		o)
			run_phyloplot="${OPTARG}" 
			;;
		p)
			viral_strains="${OPTARG}" 
			;;
		q)
			make_tree="${OPTARG}" 
			;;		
	esac
done

shift $((OPTIND-1))

#Initialising a few paths as variables, removed Documents directory variable, as script should work from working directory. Made Genomes directory containing genomes for building the tree - these are not user input ones but a template list.

working_directory=`pwd` 						#Working directory, hard coded to wherever you're working
genomeInterest=${genome_dir:-${working_directory}/Genomes_Interest} 	#The directory containing the genomes of interest
geneInterest=${gene_dir:-${working_directory}/Carriage_Genes} 	#Directory containing the genes of interest
genesForTree=${tree_genes:-${working_directory}/Tree_Genes} 		#Directory for the genes to make the tree
genomesForTree=${tree_genomes:-${working_directory}/Tree_Genomes}	#The directory for the genomes used to populate the tree
plots=${output_dir:-${working_directory}/output}			#The output directory (needs full path)
blastResults=$working_directory/BlastResults				#The blast results folder
rename_genomes=${rename_genomes:-"no"}					#Indication of whether the genomes need to be renamed, this should be a file in wd
phylogroup=${phylogroup_yn:-"no"}					#Indication of whether PhyloGroup has been used beforehand or not
grouping=${grouping_yn:-"no"}						#Indication of whether a group_list.txt file exists in wd
mac=${mac:-"no"}							#Indication of whether user is on Mac or not
rename_genes=${rename_genes:-"no"}					#Indication of whether user wishes to rename genes
phylogenes=${phylogenes:-"no"}						#Indication of whether the user used PhyloGenes beforehand or not
genes_oi=${genes_oi:-"no"}						#Indication of whether the user will have genes of interest in file or not
genomes_oi=${genomes_oi:-"no"}						#Indication of whether the user will have genomes of interest or not
sample_sheet=${sample_sheet:-"no"}					#Indication of whether the user has a sample sheet pre-PhyloBuild.
run_phyloplot=${run_phyloplot:-"no"}					#Indication of whether the user is needing to run visualisation straight after.
viral_strains=${viral_strains:-"no"}					#Indication of whether the user is running an analysis with viral strains 
make_tree=${make_tree:-"yes"}						#Indication of whether the user wants to generate a tree or just run the carriage genes step

#Indicates to user what options and parameters they've chosen and halts progress for 5 seconds, giving them time to confirm that they are correct
echo "##### Gathered user parameters - do these look right to you? #####"
echo "Working directory: ${working_directory}"
echo "Genome interest directory: ${genomeInterest}"
echo "Carriage genes directory: ${geneInterest}"
echo "Genes for tree file: ${genesForTree}"
echo "Genomes for tree file: ${genomesForTree}"
echo "Output directory: ${plots}"
echo "Rename genomes: ${rename_genomes}"
echo "Used PhyloGroup: ${phylogroup}"
echo "Using grouping for plot: ${grouping}"
echo "Using Mac: ${mac}"
echo "Rename genes: ${rename_genes}"
echo "Used PhyloGenes: $phylogenes"
echo "Using strains of interest?: $genomes_oi"
echo "Using genes of interest?: $genes_oi"
echo "Run PhyloPlot after generating Tree file?: $run_phyloplot"
echo "Are you using viral strains?: $viral_strains"
echo "Are you generating a tree (if no, then you are just running to check gene carriage)?: $make_tree"

sleep 10


echo "##### Step 1 - File parsing, moving, and format checks #####"
#Making some directories - if they exist already something may have gone wrong with clearing them last time - script may have been ended earlier. 
mkdir -p $plots
mkdir -p $blastResults #ABSOLUTELY REQUIRED TO BE REMOVED - if you change the number of genes between runs, then tree will use any old genes left behind

#I also set it so that if you have used PhyloGroup then you must have a group_list.txt, as this was one of the outputs - so if phylogroup is set to "yes" then we automatically set grouping to "yes" as well.
if [ $phylogroup == "yes" ]
then
	grouping="yes"
	cp $working_directory/listGenomes.txt $plots/listGenomes.txt
fi

#This checks whether input files are in the correct format or not, if files aren't in FASTA format an error will be returned.
python3 py/CheckFileFormat.py $genomesForTree $genesForTree $geneInterest $genomeInterest $genes_oi $genomes_oi 

exit_status=$?  # store the exit status for later use

#Check if exit status was returned, and stop the process
if [ $exit_status -ne 0 ]; then
    echo "An input file  is not in FASTA format, please reexamine your input files and make sure they are in FASTA format"
    exit $exit_status  # exit the bash script with the same status
fi
# Or continue
echo "All is good, end of format check"


#And gene match list - again this would be different with different genes/genomes.
echo "Removing gene match list from $plots - if it exists"
rm $plots/geneMatches.txt

#Check if the user wants to have genomes of interest included, then rename if using different extensions.
if [ $genomes_oi == "yes" ]
then
	cd $genomeInterest

#I make some name changes. For example, changes EC0_06_19646_7#6.contigs_velvet.fa to EC0_06_19646_7.fasta. We only change extension to .fasta because the genome folder uses these extensions, which makes it easier in the script for the loops.
#It's been made a little bit more general - if the gene names are given in .fa format then we switch these to .fasta for convenience.
	echo "Converting .fa to .fasta - this may not work if the file is already in .fasta - don't worry in this case. In the future, consider using .fasta from the beginning."
	for filename in *.fa; do 
	   mv "${filename}" "${filename%#*.fa}"
	   mv ${filename} "${filename%.fa}.fasta"
   	echo "Changed $filename to have .fasta extension"
	done
fi

#Change back to original directory
cd $working_directory

#If the user indicated that there's strains in the genomes of interest folder, then transfer these to the same folder as your other strains.
if [ $genomes_oi == "yes" ]
then
	echo "Moving genomes added in to template genomes folder..."
	cp -a $genomeInterest/*.fasta $genomesForTree

	#We make a list of the genomes we have added in -this is used to both remove the genomes from the default genome folder (in case new genomes are added, the old ones that were tested are removed) and for R plot.
	cd $genomeInterest 
	echo "Making list of genomes added in..."
	ls *.fasta | sed 's/.fasta//g' > $plots/genomesAdded.txt
fi

#Check if the user has included genes of interest as a parameter.If yes, we'll then make a txt list of the genes of interest.
if [ $genes_oi == "yes" ]
then
	cd $working_directory
	echo "Making list of genes of interest..."
	sed -n -e '/>/p' $geneInterest/genesOI.txt | sed 's/^.*>//' > $plots/genesOI.txt
fi

#Change to place where genomes are stored and make a list of genome file names - if phylogroup was run (phylogroup = yes) then this will already exist.
if [ $phylogroup == "no" ]
then
	cd $genomesForTree
	echo "Making list of genomes used to make the tree..."
	ls *.fasta | sed 's/.fasta//g' > $plots/listGenomes.txt
fi

cd $working_directory

#Below, we are running step 2. We run a blastn where the geneList.fasta (or txt if the user used PhyloGenes) genes are queries, the genomes are the subject, and use coverage of 80% and identity of 70%. For each strain, this generates a fasta format file where the gene name is given, then the matching sequence from the strain. If the user has indicated that they have genes of interest to run through, we run another blastn - this time using the genes of interest as a query, and still the strains are a subject. In this case, we just grab the gene names of those that have a match in the strain. We then take the lines in the file, add the genome name in the same line as the gene name, then put that file in the plot_extras directory.
echo "##### Step 2: blastn for tree (and optional blastn for genes of interest) #####"

cd $working_directory


if [ $make_tree == "no" ]
then
	for x in `cat $plots/listGenomes.txt`
	do
		if [ $genes_oi == "yes" ]
		then	
			if [ $viral_strains == "no" ]
			then
				#This is the blastn just to identify which genomes have a match in the genes that are added in by the user.
				blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
			
			else
				blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 20 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
			fi
			
			#We then add the genome name in every line of the file, to indicate genome matches, and save this in plot_extras
			for f in blastMatch.txt
			do 
				cat blastMatch.txt | sed "s/$/\t$x/" >> $plots/geneMatches.txt
				echo "Testing $x for carriage" 
			done	
		fi
	done
	
	echo "You opted not to make a tree, and so the process will end here"
	sleep 5
	
	exit
fi

#We need to be aware of which geneList file to use, which is determined by running PhyloGenes or not.
if [ $phylogenes == "no" ]
then 
	for x in `cat $plots/listGenomes.txt`
	do

		#We run the blastn for the genes to make the tree and save in BlastResults
		echo "Running blastn on $x to identify gene carriage and sequences for template genes (for the tree)"

		if [ $viral_strains == "no" ]
		then		
			blastn -query $genesForTree/geneList.fasta -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid sseq pident" | sed 's/^\(.\{0\}\)/\1>/' | tr '\t' '\n' >  $blastResults/$x.fasta 
		else
			blastn -query $genesForTree/geneList.fasta -subject $genomesForTree/$x.fasta -perc_identity 70 -qcov_hsp_perc 20 -outfmt "6 qseqid sseq pident" -task blastn | sed 's/^\(.\{0\}\)/\1>/' | tr '\t' '\n' >  $blastResults/$x.fasta 
		fi
		
		if [ $genes_oi == "yes" ]
		then	
			if [ $viral_strains == "no" ]
			then
				#This is the blastn just to identify which genomes have a match in the genes that are added in by the user.
				blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
			
			else
				blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 20 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
			fi
			
			#We then add the genome name in every line of the file, to indicate genome matches, and save this in plot_extras
			for f in blastMatch.txt
			do 
				cat blastMatch.txt | sed "s/$/\t$x/" >> $plots/geneMatches.txt 
			done	
		fi
	done
	#This is a new addition, because we often get multiple matches for one gene, I put a restriction on this. Essentially, if there are multiple gene matches for one gene in the tree - then the highest match gets kept. For now. Will further test performance with this criteria
	python3 py/blast_formatting_updated.py $blastResults 

else #Now I've removed the alternative IF statement and put the ELSE instead.
#This version just uses the different file name. 
	for x in `cat $plots/listGenomes.txt`
	do
		if [ $viral_strains == "no" ]
		then
			#We run the blastn for the genes to make the tree and save in BlastResults
			echo "Running blastn on $x to identify gene carriage and sequences for template genes (for the tree)"
			blastn -query $genesForTree/geneList.fasta -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid sseq pident" | sed 's/^\(.\{0\}\)/\1>/' | tr '\t' '\n' >  $blastResults/$x.fasta 
		else
			blastn -query $genesForTree/geneList.fasta -subject $genomesForTree/$x.fasta -perc_identity 70 -qcov_hsp_perc 20 -outfmt "6 qseqid sseq pident" -task blastn | sed 's/^\(.\{0\}\)/\1>/' | tr '\t' '\n' >  $blastResults/$x.fasta 
		fi
		
		if [ $genes_oi == "yes" ]
		then
			if [ $viral_strains == "no" ]
			then
				#This is the blastn just to identify which genomes have a match in the genes that are added in by the user.
				blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
			else
				blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 20 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
			fi
			#We then add the genome name in every line of the file, to indicate genome matches, and save this in plot_extras
			for f in blastMatch.txt
			do 
				cat blastMatch.txt | sed "s/$/\t$x/" >> $plots/geneMatches.txt 
			done	
		fi
	done
	python3 py/blast_formatting_updated.py $blastResults
fi

#This is just a clean up after running the python script
cd $blastResults
shopt -s extglob
rm !(new*.fasta)
cd $working_directory

for filename in $blastResults/*.fasta
do
	mv "$filename" "${filename//new_/}"
done

echo "BLAST has been run over all the genomes."

cd $working_directory

#This is step 3, where we run alignment on the gene sequences found in the first blastn. A more detailed description of the script can be found in the steps at the start of the script.

echo "##### Step 3: Running alignment #####"

perl align_blast.pl BlastResults ALIGNMENTS
find . -name 'ALIGNMENTS.*.fasta' -delete #removes unecessary .fasta files generated from script
#One slightly annoying thing about this perl script is that concatenation back to the final file for some reason seems to mess with the number of sequences, such that if we ended up using different sequence for the geneList file, then occasionally we get an error which says that the length of sequences differs - which fails the script. Therefore, in addition to the align_blast.pl, I have to run MAFFT on the final concatenated file.
#mafft --auto Final.ALIGNMENTS.aln > Final.Aligned.aln

#This is the python script used to run alignment using MAFFT
#cd $working_directory/py
#python3 Alignment_Script.py $working_directory $plots
#cd $working_directory



#REMEMBER TO SWITCH BACK TO Final.Aligned.aln IF SWITCHING BACK TO ABOVE VERSION 
#Below is step 4, where we convert the resulting alignment file to a phylip format for tree producing in phyml.
echo "##### Step 4: Changing .aln format to .phy #####"
#java -jar $working_directory/ALTER/alter-lib/target/ALTER-1.3.4-jar-with-dependencies.jar -cg -i final_alignment.tsv -ia -io Linux -o phylipFor.phy -of PHYLIP -oo Linux -op PhyML

#java -jar $working_directory/ALTER/alter-lib/target/ALTER-1.3.4-jar-with-dependencies.jar -cg -i Final.Aligned.aln -ia -io Linux -o phylipFor.phy -of PHYLIP -oo Linux -op PhyML

trimal -in Final.ALIGNMENTS.aln -phylip -out phylipFor.phy -gappyout -keepheader

#cp final_alignment.tsv $plots/final_alignment.tsv

#Step 5, where we run phyml. We use a bootstrap repition of 100.
echo "##### Step 5: Producing Phylogenetic Tree #####"

#And now we need to know whether the user is using a Mac or not - if you are, make sure you've downloaded Phyml and extracted to your working directory.
if [ $mac == "yes" ]
then
	./PhyML-3.1/PhyML-3.1_macOS-MountainLion -i phylipFor.phy -b 100 --quiet
fi

if [ $mac == "no" ]
then
	phyml -i phylipFor.phy -b 100 --quiet
fi

#We copy the resulting newick file format tree data into the plot extras for R.
cp phylipFor.phy_phyml_tree.txt $plots/treeDataFile.txt

#Step 6 is where the R script is used to generate the plot - in case the user uses the same output folder multiple times I've set it so that the name is changed to a unique number each time. The output will also include the R script with the args switched with the user path - so the user can immediately run the same R script after running the script in case they wish to make any adjustments to the script/plot. But, in this new version we've separated tree generation and visualisation into two different steps.

#Check if user has already pre-generated a sample sheet
if [ $sample_sheet == "no" ]
then
	echo "##### Making sample sheet #####"
	python3 py/Sample_Sheet_Generation.py $working_directory $grouping $genomes_oi $rename_genomes $plots
fi

#Check if they need to run visualisation right after or are going to do it later.
if [ $run_phyloplot == "yes" ]
then
	echo "##### Generating Plot through R #####"
	Rscript generate_plot.r $plots $rename_genomes $grouping $rename_genes $genes_oi $genomes_oi
	cp $plots/finalPlot.svg $plots/finalPlot.$$.svg
	rm $plots/finalPlot.svg
	echo "Rscript ran, output finalPlot.$$.EMF should be in $plots"

	#This is simply the python script which switches args with the user's paths, options, etc.
	python3 py/rScriptParameters.py $plots $rename_genomes $working_directory $grouping $rename_genes $genes_oi $genomes_oi

fi

#Save parameters to txt file in case they need to be used for visualisation
echo "$plots $rename_genomes $working_directory $grouping $rename_genes $genes_oi $genomes_oi" > r_script_parameters.txt

#Remove the tmp directory, used in step 5 and no longer needed.
#rmdir tmp

if [ $genomes_oi == "yes" ]
then
	cd $genomesForTree
	#Removing the genomes of interest from the default genome list in Genomes
	echo "Removing genomes of interest from general genome file"
	for f in `cat $plots/genomesAdded.txt`
	do 
    		rm ${f}.fasta
	done
fi

cd $working_directory

#Clear the blast results directory
echo "Clearing Blast Result directory in case of changes to gene numbers"

#rm -r $blastResults

#Clear the alignments file (should stay constant but in case of updates old genes used will still show up without this)
echo "Removing ALIGNMENTS.Alignments file (if it exists)."
rm -r ALIGNMENTS.Alignments
rm Final.ALIGNMENTS.aln
rm Final.Aligned.aln
rm -r Alignment_Files
rm -r BlastResults
rm final_alignment.tsv

cd $working_directory
#Deleting extra files generated by programs that are no longer necessary.
echo "Deleting extra files that are unecessary after the script..."
find . -maxdepth 1 -name 'phylipFor*' -delete

end=$SECONDS
runtime=$((end - $start))
mins=$((runtime / 60))
secs=$((runtime % 60))
hour=0
while (( $mins >= 60 ))
do
	hour=$((hour+1))
	mins=$((mins-60))
done

if (( $hour >= 1 ))
then
	echo "Script took ${hour} hour, ${mins} minutes, and ${secs} seconds" 
else
	echo "Script took ${mins} minutes and ${secs} seconds" 
fi
