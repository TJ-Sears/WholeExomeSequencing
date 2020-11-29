#!/bin/bash

#Wrapper for consensus variant calling from input paired tumor/normal sorted deduped bam files
#required arguments 1-sample list input file 2-output directory #sample input file should be in the format:
#sample_name-T1	fullpath/sample_tumor_sorted_deduped.bam	fullpath/sample_normal_sorted_deduped.bam
#each sample name should be unique, use "-T2" etc. to designate multiple samples from the same patient
#prior to running - modify the Settings below

##########     Settings - please modify this section as desired    #################
CNV_kit_switch="on" #"on" will run CNVkit for copynumber
CNV_reference_switch="off" #on will generate pooled reference for all samples in directory
CNV_target_antitarget="off" #will create target/antitarget coverage

export selector_bed=/scratch/users/tsears/WholeExome/MouseRef/Twist_Mouse_Exome_Target_Rev1_7APR20.bed #enter full path of the capture selector bed; default is MedExome
export selector_padded_bed=/scratch/users/tsears/WholeExome/MouseRef/mm10_exons.bed #enter full path of the capture selector bed; default is MedExome
export cnv_selector=/scratch/users/tsears/WholeExome/MouseRef/mm10.padded.bed #modified capture .bed file with only gene name in 4th column; needed for better visualization on cnv diagram- if unavailable can use the selector bed file

#######   Should not require modification past this line    ##########################################

#required scripts
export scripts_dir="/scratch/users/tsears/WholeExome/exome_pipeline_scripts" #dir containing all associated scripts
	#MSK_process_variants.sh tag
	#MSK_process_annovar_output.sh annovar.multianno.txt annovar.coding.file tag
	#MSK_run_varscan.sh	
	#merge_VCFs.R
	#MSK_run_stats.sh
	#plot_coverage.R 
	#MSK_polysolver.sh

#Dependencies
export GATK_index_dir="/scratch/users/tsears/WholeExome/exome_pipeline_scripts/indexes"
export annovar_dir="/scratch/users/tsears/WholeExome/annovar_2015-03-22" #dir to annovar scripts
export samtools="$scripts_dir/dependencies/samtools-0.1.19/samtools" #samtools command
export samtools_dir="$scripts_dir/dependencies/samtools-0.1.19" #samtools dir
export varscan2_4_1_dir="$scripts_dir/dependencies/varscan-master"
export bam_readcount="$scripts_dir/dependencies/bam-readcount/bin/bam-readcount"
export strelka_dir="$scripts_dir/dependencies/strelka_workflow-1.0.14"
export polysolver_dir="/scratch/users/tsears/WholeExome/polysolver"
export polysolver="/scratch/users/tsears/WholeExome/exome_pipeline_scripts/MSK_polysolver.sh"  
export GATK_dir="$scripts_dir/dependencies/gatk-4.1.9.0/"
export strelka_exome_config="$strelka_dir/etc/strelka_config_bwa_exome.ini" #this file may need to be modified for non-exome variant calling
export mutect_jar="$scripts_dir/dependencies/mutect/mutect-1.1.7.jar"
export cnvkit_script="$scripts_dir/dependencies/cnvkit-0.9.6/cnvkit.py"
export factera_script="/scratch/users/tsears/WholeExome/cappseq/scripts-v3/vars-fusions-factera.pl"
export plot_coverage="$scripts_dir/plot_coverage.R"
#export java_dir="/usr/lib/jvm/java-openjdk/jre/bin"
export mutect_dir="$scripts_dir/dependencies/mutect"
export novoalign_dir="$scripts_dir/dependencies/novocraft"
export PATH=$PATH:$scripts_dir/dependencies/ensembl-vep/:$scripts_dir/SuppScripts/bin/
export vep="$scripts_dir/dependencies/ensembl-vep/vep"
export convert2bed_script="$scripts_dir/dependencies/bin/convert2bed"
export bam2freq_script="$scripts_dir/bam-snvfreq-single.py"
export annotate_vars_script="$scripts_dir/SuppScripts/jchabon_scripts/annotate_var_list.py"

#indices/references
export reference_fasta="/scratch/users/tsears/WholeExome/MouseRef/mm10.fa" 
export db="/scratch/users/tsears/WholeExome/annovar_2015-03-22/mousedb" #dir to annovar human db files
export dbsnp_mutect="$scripts_dir/indices/dbsnp_132_b37.leftAligned.vcf.gz" #for mutect
export cosmic_vcf="$scripts_dir/indices/b37_cosmic_v54_120711.vcf.gz" #for mutect
export exons="/scratch/users/tsears/WholeExome/RefSeq_Gencodev17_022314.allexons.bed" #for factera
export genome_2bit="/scratch/users/tsears/WholeExome/hg19.2bit" #for factera
export cnv_access_file="$scripts_dir/dependencies/cnvkit-0.9.6/data/access-5kb.mm10.bed" #for cnvkit
export GATK_indel_standards_file="$scripts_dir/indices/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
export GATK_known_sites_file="$scripts_dir/indices/dbsnp_138.hg19.vcf"
export pon="/scratch/users/tsears/WholeExome/MouseRef/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
export germline="/scratch/users/tsears/WholeExome/MouseRef/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
export germline2="$scripts_dir/indices/somatic-b37_af-only-gnomad.raw.sites.vcf"
export mm10dbsnp="/scratch/users/tsears/WholeExome/MouseRef/mm10dbsnp.vcf"
export refFlat="/scratch/users/tsears/WholeExome/MouseRef/refFlat.txt"

if [ "$#" -ne "2" ]; then
	echo "Incorrect number of arguments.  Please run as"
	echo "MSK_exome.sh sample_list_file.txt output_directory"
	input_mode="multi"
	exit
else
	export input_file=$1
	export starting_dir=$(echo $2 | sed 's/\/$//') #removes last slash if present in directory name
	if [ $2 == "." ]; then
		starting_dir=$PWD
	fi
	export working_dir="$starting_dir/exome_pipeline"
fi


#functions
create_dir_structure() { #make directories; Argument 1 is space delimited list of samples
mkdir -p "$working_dir"
mkdir -p "$working_dir/bam_files"
mkdir -p "$working_dir/variants"
export sample_table="$working_dir/sample_table.temp"
awk '{print $1,$1,$2,$3}' $input_file | awk '{gsub(/-T.*/,"",$2)}1' > $working_dir/sample_table.temp #sample list file for downstream scripts
export tumor_samples=$(awk '{print $1}' $input_file)
export normal_samples=$(awk '{print $1}' $input_file | sed 's/-T.*//')

for sample in $tumor_samples
	do
	mkdir -p "$working_dir/$sample"
done
export bam_dir="$working_dir/bam_files"
export HLA_dir="$working_dir/HLA_files"
export run_stat_dir="$working_dir/run_stats"
export variant_dir="$working_dir/variants"
export cnv_dir="$working_dir/CNVkit"
}

make_bam_links() {
while read -r line
	do
		tumor_tag=$(echo $line | cut -d" " -f1)
		normal_tag=$(echo $line | cut -d" " -f2)
		tumor_bam=$(echo $line | cut -d" " -f3)
		normal_bam=$(echo $line | cut -d" " -f4)
		
		ln -fs $tumor_bam $bam_dir/$tumor_tag.tumor.sorted.deduped.bam
		ln -fs $normal_bam $bam_dir/$normal_tag.normal.sorted.deduped.bam
	if [ -e $tumor_bam.bai ]
		then	
			ln -fs $tumor_bam.bai $bam_dir/$tumor_tag.tumor.sorted.deduped.bam.bai
		else
			$samtools index $tumor_bam > $bam_dir/$tumor_tag.tumor.sorted.deduped.bam.bai
	fi
	if [ -e $normal_bam.bai ]
		then
			ln -fs $normal_bam.bai $bam_dir/$normal_tag.normal.sorted.deduped.bam.bai
		else
			$samtools index $normal_bam > $bam_dir/$normal_tag.normal.sorted.deduped.bam.bai
	fi
done < $sample_table
}

run_mutect_caller() {

mutect_jar=$1
reference_fasta=$2
pon=$3
germline=$4
working_dir=$5
tag=$6

tumorbam=/scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/$tag/${tag}_L.sorted.readgroup.bam
normalbam=/scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/$tag/${tag}_T.sorted.readgroup.bam
echo $reference_fasta
echo $tumorbam
echo $normalbam
echo $6

cd $working_dir/$tag

#Run Mutect2
$GATK_dir/gatk Mutect2 -R $reference_fasta -I $tumorbam -I $normalbam -O $tag.mutect.raw.out -normal 20 -L $selector_bed
$GATK_dir/gatk FilterMutectCalls -V $tag.mutect.raw.out -O $tag.mutect.pass.raw.out -R $reference_fasta

#Filter for reads that PASS filters
bcftools view -f PASS $tag.mutect.pass.raw.out > $tag.passed.vcf
#Filter out sections that are within RepeatMasker defined repeat regions
bedtools intersect -loj -f 0.50 -a $tag.passed.vcf -b $scripts_dir/dependencies/RepeatMasker/Outs/mm10.masked.bed > $tag.mutect.passed.vcf

#Old commands for consensus runs
awk -F'[\t:]' '{print 100*$18}' $tag.mutect.passed.vcf > $tag.mutect.vaf.temp
paste $tag.mutect.passed.vcf $tag.mutect.vaf.temp > $tag.mutect.vaf.pre.temp
awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$NF, "mutect"}' $tag.mutect.vaf.pre.temp > $tag.mutect.convert2annovar.input

cd $working_dir
}
export -f run_mutect_caller

call_mutect_variants() {
for i in $tumor_samples
do
run_mutect_caller $mutect_jar $reference_fasta $pon $germline $working_dir $i
done
}

cnv_analysis() {
#create sample file with switch - for single vs multiplex; creates sample file rather than just variable for future scripts
mkdir -p "$starting_dir/exome_pipeline/CNVkit"
cd $cnv_dir

 pooled_reference=$scripts_dir/dependencies/Reference.cnn
 tumorbam=/scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/*/*_L.sorted.readgroup.bam
 normalbam=/scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/*/*_T.sorted.readgroup.bam
 allbams=/scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/*/*.sorted.readgroup.bam

 $cnvkit_script autobin $allbams -t $cnv_selector -g $cnv_access_file --annotate $refFlat --short-names
 
 autobins=${cnv_selector##*/}

            if [ $CNV_target_antitarget == "on" ] #turn on for first time run
 	        then
 		for i in $allbams
 		do	 
 		    samp=${i##*/}
 		    $cnvkit_script coverage $i ${autobins%.*}.target.bed -o ${samp%%.*}.targetcoverage.cnn
 		    $cnvkit_script coverage $i ${autobins%.*}.antitarget.bed -o ${samp%%.*}.antitargetcoverage.cnn
 		done
 		fi

 	    if [ $CNV_reference_switch == "on" ] #turn on for first time run
 		then
 		$cnvkit_script reference *_L*coverage.cnn -f $reference_fasta -o $scripts_dir/dependencies/Reference.cnn
 		fi

        pooled_reference=$scripts_dir/dependencies/Reference.cnn

 	for i in $tumorbam
 	do
 	        samp=${i##*/}
 	    	echo $i			
 		$cnvkit_script fix ${samp%%.*}.targetcoverage.cnn ${samp%%.*}.antitargetcoverage.cnn $pooled_reference -o ${samp%%.*}.tumor.cnr
 		$cnvkit_script segment ${samp%%.*}.tumor.cnr -o ${samp%%.*}.tumor.cns
 		$cnvkit_script scatter ${samp%%.*}.tumor.cnr -s ${samp%%.*}.tumor.cns -o ${samp%%.*}.tumor.scatter.pdf
 		$cnvkit_script diagram -s ${samp%%.*}.tumor.cns ${samp%%.*}.tumor.cnr 
 		$cnvkit_script genemetrics ${samp%%.*}.tumor.cnr > ${samp%%.*}.genemetrics_no_cns.txt #-s ${samp%%.*}.tumor.cns > ${samp%%.*}.genemetrics.txt
 		$cnvkit_script call ${samp%%.*}.tumor.cns -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o ${samp%%.*}.call.cns
 		$cnvkit_script breaks ${samp%%.*}.tumor.cnr ${samp%%.*}.tumor.cns > ${samp%%.*}.gene.breaks.txt
 		cp ${samp%%.*}.call.cns /scratch/users/tsears/WholeExome/CNV_11_17/
 		cp ${samp%%.*}.tumor.scatter.pdf /scratch/users/tsears/WholeExome/CNV_11_17/
 		cp ${samp%%.*}.tumor.cns /scratch/users/tsears/WholeExome/CNV_11_17/
 		cp ${samp%%.*}.tumor.cnr /scratch/users/tsears/WholeExome/CNV_11_17/
 		cp ${samp%%.*}.genemetrics_no_cns.txt /scratch/users/tsears/WholeExome/CNV_11_17/
 		cp ${samp%%.*}.gene.breaks.txt /scratch/users/tsears/WholeExome/CNV_11_17/
 	done

 	$cnvkit_script heatmap *.cnr -d -o All.tumors.pdf
 	$cnvkit_script heatmap *.cnr -d -c chr2 -o All.chr2.pdf
}

run_VEP() {
echo "Running VEP"
for tag in $tumor_samples
do
	cd $working_dir/$tag/
	variant_output=$working_dir/$tag/$tag.mutect.passed.vcf
	$vep -i $variant_output --vcf --species mus_musculus -o $tag.variant.annotated.vcf.temp --cache --force_overwrite --format vcf --pick
	sed -i -e 1,3d $tag.variant.annotated.vcf.temp
	awk -v OFS="\t" 'FNR==1{print $1,$2,$3,$4,$5,$6,$7,$8,"INFO2","INFO3","INFO4"}1' $tag.variant.annotated.vcf.temp > $working_dir/variants/$tag.variant.annotated.vcf
	sed -i -e 2d $working_dir/variants/$tag.variant.annotated.vcf
done
}


annotate_depth() {
for tag in $tumor_samples
do
    cd $working_dir/$tag/
    variant_output_for_bed=$working_dir/$tag/$tag.mutect.passed.vcf
    echo "Generating BED file for FREQ generation"
    $convert2bed_script -i vcf < $variant_output_for_bed > $tag.variant_file.bed
    echo "Done generating BED file"

    tumorbam=/scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/$tag/${tag}_L.sorted.readgroup.bam
    normalbam=/scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/$tag/${tag}_T.sorted.readgroup.bam
    echo "Starting FREQ Generation"

    #Generate Tumor Freq
    python $bam2freq_script $tumorbam $reference_fasta $tag.variant_file.bed 
    #Generate Normal Freq
    python $bam2freq_script $normalbam $reference_fasta $tag.variant_file.bed
    
    #Copying and renaming FREQ files because the annotate_vars script expects it in this format
    cp /scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/$tag/${tag}_T.sorted.readgroup.freq.paired.Q30.txt $working_dir/$tag/${tag}_T.freq.paired.Q30.txt
    cp /scratch/users/tsears/WholeExome/Mouse/bwa_aligned_deduped/$tag/${tag}_L.sorted.readgroup.freq.paired.Q30.txt $working_dir/$tag/${tag}_L.freq.paired.Q30.txt
    
    #Annotate with Normal Reads
    python $annotate_vars_script $working_dir/variants/$tag.variant.annotated.vcf germline $working_dir/$tag/
    python $annotate_vars_script $working_dir/variants/$tag.variant.annotated.annoGermline.txt tumor $working_dir/$tag/ 
    
    cp $working_dir/variants/$tag.variant.annotated.annoGermline.annoTumor.txt /scratch/users/tsears/WholeExome/VariantResults_11_13/x$tag.variant.results
done
}
    


##########################Script starts here########################################
start=$(date)
start_time=$(date +%s)
runid=$(date | awk '{print $2$3$NF$4}')
samplefile=samples.$runid.txt
logfile=log.$runid.txt

create_dir_structure #creates directory structure
echo -e "Start time $start" > $working_dir/$logfile

make_bam_links #sym link to input bams

# SECONDS=0   #mapping functions not used in this script
# cd $bam_dir
# map_samples $tumor_samples $normal_samples #maps with BWA mem to sam
# mapping_duration=$(displaytime $SECONDS)
# echo "Mapping duration =$mapping_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# sort_bams #sorts sam and converts to bam
# sort_duration=$(displaytime $SECONDS)
# echo "Sorting duration =$sort_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# remove_dups #removes dups with picard
# rm_dup_duration=$(displaytime $SECONDS)
# echo "Deduping duration=$rm_dup_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# echo "realigning bam files" | tee -a $working_dir/$logfile
# realign_bams #realigns around indels with picard
# realign_duration=$(displaytime $SECONDS)
# echo "realigning duration=$realign_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# if [ $recalibrate_switch == "on" ]	
# 	then
# 	recalibrate_bases #recalibrates bases
# 	recalibrate_bases_duration=$(displaytime $SECONDS)
# 	echo "recalibrating bases duration=$recalibrate_bases_duration" | tee -a $working_dir/$logfile
# fi

# SECONDS=0
# echo "Indexing bam and creating symlinks" | tee -a $working_dir/$logfile
# index_and_symlink #creates bam indexes and symlinks to sample directories
# index_duration=$(displaytime $SECONDS)
# echo "Indexing and symlink duration=$index_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# echo "Starting Varscan variant calling" | tee -a $working_dir/$logfile
# call_varscan_variants  #calls variants using varscan plus realign and fpf filter
# varscan_duration=$(displaytime $SECONDS)
# echo "Varscan duration=$varscan_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# echo "Starting Mutect variant calling" | tee -a $working_dir/$logfile
# call_mutect_variants #calls variants using mutect  
# mutect_duration=$(displaytime $SECONDS)
# echo "Mutect duration=$mutect_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# echo "Starting Strelka variant calling" | tee -a $working_dir/$logfile
# call_strelka_variants #calls variants using strelka
# strelka_duration=$(displaytime $SECONDS)
# echo "Strelka duration=$strelka_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# echo "Calling consensus_variants"
# consensus_variants #add VAFs, reformat variant call outputs.  R script to determine consensus calls, create the annovar input files
# echo "finished consensus_variants"
# consensus_duration=$(displaytime $SECONDS) 
# echo "Consensus calling duration=$consensus_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# annotate_annovar #annotate all calls and consensus calls
# annovar_duration=$(displaytime $SECONDS)
# echo "Annovar duration=$annovar_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# run_VEP #run VEP and annotate variant output files
# echo "finished  VEP annotation"
# VEP_duration=$(displaytime $SECONDS)
# echo "VEP duration=$VEP_duration" | tee -a $working_dir/$logfile

# SECONDS=0
# annotate_depth
# echo "Done annotating tumor and germline depths to variant output"
# annotate_depth_dur=$(displaytime $SECONDS)
# echo "annotate depth duration=$annotate_depth_dur" | tee -a $working_dir/$logfile   

# SECONDS=0
# if [ $QC_stats_switch == "on" ]
# 	then
# 	echo "Performing QC calculations" | tee -a $working_dir/$logfile
# 	perform_runstats  #calculates coverage and ontarget
# 	runstats_duration=$(displaytime $SECONDS)
# 	echo "Calculating QC stats duration=$runstats_duration" | tee -a $working_dir/$logfile
# fi

# SECONDS=0
# if [ $hla_switch == "on" ]
# 	then
# 	echo "Calling HLA genotypes and mutations" | tee -a $working_dir/$logfile
# 	call_HLA  #calls HLA type, HLA mutations, and annotates
# 	HLA_duration=$(displaytime $SECONDS)
# 	echo "Calculating HLA duration=$HLA_duration" | tee -a $working_dir/$logfile
# fi

# SECONDS=0
# process_annovar #annotation processing - reformat annotations, output annotations to file, predict coding changes, create fasta, create iedb input, create realignment file for post iedb annotation
# process_annovar_duration=$(displaytime $SECONDS)
# echo "Annovar processing duration=$process_annovar_duration" | tee -a $working_dir/$logfile

#Copy number
SECONDS=0
if [ $CNV_kit_switch == "on" ]
	then
	echo "Detecting copy number variation" | tee -a $working_dir/$logfile
	cnv_analysis #run cnvkit of each tumor against pooled germlines background
	cnv_duration=$(displaytime $SECONDS)
	echo "Copy number calling duration=$cnv_duration" | tee -a $working_dir/$logfile
fi

###Run factera
# SECONDS=0
# if [ $factera_switch == "on" ]
# 	then
# 	run_factera
# 	factera_duration=$(displaytime $SECONDS)
# 	echo "Factera duration=$factera_duration" | tee -a $working_dir/$logfile
# fi

#####End script
stop=$(date)
echo -e "Completed time $stop" >> $working_dir/$logfile
stop_time=$(date +%s)
duration=$(($stop_time-$start_time))
total_duration=$(displaytime $duration)
echo "Total pipeline duration=$total_duration" | tee -a $working_dir/$logfile
exit
