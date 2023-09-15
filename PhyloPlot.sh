#!/bin/bash

#Here's where the R script is ran, the parameters are gathered after running PhyloTree regardless of visualisation being run or not
start=$SECONDS

plots=`awk '{print $1}' r_script_parameters.txt`
rename_genomes=`awk '{print $2}' r_script_parameters.txt`
working_directory=`awk '{print $3}' r_script_parameters.txt`
grouping=`awk '{print $4}' r_script_parameters.txt`
rename_genes=`awk '{print $5}' r_script_parameters.txt`
genes_oi=`awk '{print $6}' r_script_parameters.txt`
genomes_oi=`awk '{print $7}' r_script_parameters.txt`

echo "##### Generating Plot through R #####"
Rscript generate_plot.r $plots $rename_genomes $grouping $rename_genes $genes_oi $genomes_oi
cp $plots/finalPlot.svg $plots/finalPlot.$$.svg
rm $plots/finalPlot.svg
echo "Rscript ran, output finalPlot.$$.EMF should be in $plots"

#This is simply the python script which switches args with the user's paths, options, etc.
python3 py/rScriptParameters.py $plots $rename_genomes $working_directory $grouping $rename_genes $genes_oi $genomes_oi

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
