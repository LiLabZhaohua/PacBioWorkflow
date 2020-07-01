#!/usr/bin/env bash
#PacbioManualFilterByCAGE.sh
#This pipeline performs CAGE correction after getting assembly
#This script only reports isoforms that can be supported by CAGE peaks
#Modified from Anqi Wang's workflow
#Version: Yu Sun, 2020-05-04
#Version: Yu Sun, 2020-05-15, keep GPD output without * in column 7-8

if [ ! -n "$2" ]
then
  echo "*************************************************************************************************"
  echo "*                    Welcome to use PacbioManualFilterByCAGE pipeline                           *"
  echo "*This pipeline performs CAGE correction after getting assembly                                  *"
  echo "*This script only reports isoforms that can be supported by CAGE peaks                          *"
  echo "*This is the 2020 version of this pipeline, and please run after running PacbioGetAssembly.sh   *"
  echo "*This workflow requires several critical files to run:                                          *"
  echo "*1, Genome annotation: \${genome}_refFlat_modified.gpd                                           *"
  echo "*2, CAGE peaks: \${CAGE_Peaks}                                                                   *"
  echo "*3, Update all names in InputParameters.txt file, and use that as input to this pipeline        *"
  echo "*Critical information in the InputParameters.txt file:                                          *"
  echo "*        genome=\"hg38\"                                                                          *"
  echo "*        CAGE_Peaks=\"merged.CAGE.sorted.bed6.5end.sorted.peaks\"                                 *"
  echo "*In addition, this workflow also reads the previous output from PacbioGetAssembly.sh:           *"
  echo "*   \${OutputPrefix}.LRs.ccs.declip.sam                                                          *"
  echo "*   \${OutputPrefix}.LRs.ccs.declip.gpd                                                          *"
  echo "*   \${OutputPrefix}.annotation_intact.txt                                                       *"
  echo "*   \${OutputPrefix}.StringTie.filtered.assgn.chtnm.gpd                                          *"
  echo "*   \${OutputPrefix}.StringTie.filtered.assgn.chtnm.LR.annotate_intact.txt                       *"
  echo "*Usage: `basename $0` [InputParameters] [OutputPrefix]                            *"
  echo "*Output: 1, *.Assembly.CAGECorrected.Sorted.gpd                                                 *"
  echo "*        2, *.Assembly.CAGECorrected.IntactReadsSubisoforms.txt                                 *"
  echo "*************************************************************************************************"
else

source $1
OutputPrefix=$2
echo "*************************************************************************************************"
echo "*             Welcome to use PacbioManualFilterByCAGE single PacBio data pipeline               *"
echo "Starting the pipeline with following parameters:"
echo "Genome annotation: "${genome}_refFlat.txt" or "${genome}_refFlat_modified.gpd
echo "LRs bam file: "${LR_Bam}
echo "StringTie gpd file: "${StringTie_Assembly}
echo "OutputPrefix: "${OutputPrefix}
echo "In addition, this workflow also reads the previous output from PacbioGetAssembly.sh:"
echo "   \${OutputPrefix}.LRs.ccs.declip.sam"
echo "   \${OutputPrefix}.LRs.ccs.declip.gpd"
echo "   \${OutputPrefix}.annotation_intact.txt"
echo "   \${OutputPrefix}.StringTie.filtered.assgn.chtnm.gpd"
echo "   \${OutputPrefix}.StringTie.filtered.assgn.chtnm.LR.annotate_intact.txt"

module load samtools
module load boost

echo "Step 01, Modify CAGE peaks"
awk '{if($2=="+") print $1"\t"$2"\t"$8; else print $1"\t"$2"\t"$9}' ${CAGE_Peaks} > ${OutputPrefix}.CAGE.list.txt
sort -k 1,1 -k 3,3n ${OutputPrefix}.CAGE.list.txt > ${OutputPrefix}.CAGE.list.sorted.txt

echo "Step 02, Find CAGE supported LRs, and identify novel (not in exon-usage intact LRs defined by RefSeq and SR assembly) LRs"
echo "   Reading in previous output: "${OutputPrefix}.LRs.ccs.declip.sam, ${OutputPrefix}.LRs.ccs.declip.gpd
# Find the long reads whose 5' ends are supported by CAGE peaks, record their names.
SupportByCAGE ${OutputPrefix}.CAGE.list.sorted.txt ${OutputPrefix}.LRs.ccs.declip.sam > ${OutputPrefix}.LR.SupportedByCAGE.txt

echo "   Reading in previous output: "${OutputPrefix}.annotation_intact.txt, ${OutputPrefix}.StringTie.filtered.assgn.chtnm.LR.annotate_intact.txt
# Record the names of long reads that are already determined as intact at exon usage level: both RefSeq and SR assembly
cat ${OutputPrefix}.annotation_intact.txt ${OutputPrefix}.StringTie.filtered.assgn.chtnm.LR.annotate_intact.txt | cut -f1 > ${OutputPrefix}.LR.IntactExonLevel.txt

# Remove the exon-usage-intact long reads from those whose 5'ends are supported by CAGE peaks. The retained ones are to be clustered to identify novel isoforms.
MinusSet ${OutputPrefix}.LR.SupportedByCAGE.txt ${OutputPrefix}.LR.IntactExonLevel.txt ${OutputPrefix}.LR.Novel.txt
# Pull out the GPD records of the long reads to be clustered.
SearchLine ${OutputPrefix}.LR.Novel.txt 1 ${OutputPrefix}.LRs.ccs.declip.gpd 1 > ${OutputPrefix}.LR.Novel.gpd

echo "Step 03, Clustering novel LRs, then get augmented annotation by combining RefSeq and SR assembly"
#Perform cluster on the long reads that correspond to novel isoforms. Each cluster represents a novel isoform.
sort -k 3,3 -k 5,5n ${OutputPrefix}.LR.Novel.gpd > ${OutputPrefix}.LR.Novel_sorted.gpd
ClusterUpdate ${OutputPrefix}.LR.Novel_sorted.gpd ${OutputPrefix}.LR.Novel.pre_isoform ${OutputPrefix}.LR.Novel.pre_intact

# Combine the hg38 annotation and the annotation obtained through stringtie (with gene and transcripts name changed). Make an augmented annotation.
echo "   Reading previous file: "${genome}_refFlat_modified.gpd, ${OutputPrefix}.StringTie.filtered.assgn.chtnm.gpd
cat ${genome}_refFlat_modified.gpd ${OutputPrefix}.StringTie.filtered.assgn.chtnm.gpd | sort -u | sort -k 1,1 > ${OutputPrefix}.AugmentedAnnotation.gpd   #very similar to ${OutputPrefix}.MergedIsoform.gpd
sort -k 3,3 -k 5,5n ${OutputPrefix}.AugmentedAnnotation.gpd > ${OutputPrefix}.AugmentedAnnotation_sorted.gpd

# Assign the gene names of the novel isoforms obtained by clustering.
echo "   Outputting: "${OutputPrefix}.CAGE.NovelIsoform.gpd
AssignGene ${OutputPrefix}.AugmentedAnnotation_sorted.gpd ${OutputPrefix}.LR.Novel.pre_isoform > ${OutputPrefix}.CAGE.NovelIsoform.gpd


echo "Step 04, Combine the exon-usage-level intact isoforms and the CAGE novel isoforms, obtain the assembled isoforms by GDP_annotate again, just to get clean results"
cat ${OutputPrefix}.annotation_intact.txt ${OutputPrefix}.StringTie.filtered.assgn.chtnm.LR.annotate_intact.txt | cut -f3 | sort -u > ${OutputPrefix}.Selected.AnnoIsoform.txt
SearchLine ${OutputPrefix}.Selected.AnnoIsoform.txt 1 ${OutputPrefix}.AugmentedAnnotation.gpd 2 > ${OutputPrefix}.Selected.AnnoIsoform.gpd
cat ${OutputPrefix}.Selected.AnnoIsoform.gpd ${OutputPrefix}.CAGE.NovelIsoform.gpd > ${OutputPrefix}.Assembly.gpd
sort -k 3,3 -k 5,5n ${OutputPrefix}.Assembly.gpd > ${OutputPrefix}.Assembly_sorted.gpd

# Align the long read isoforms to the obtained assembly, find out the exon-usage-level intact isoforms again. This step seems to be repetitive, just for clean results.
GPD_annotate ${OutputPrefix}.Assembly_sorted.gpd ${OutputPrefix}.LRs.ccs.declip.gpd ${OutputPrefix}.Assembly.annotate
rm -rf ${OutputPrefix}.Selected.AnnoIsoform.gpd ${OutputPrefix}.Selected.AnnoIsoform.txt

echo "Step 05, Pull out the sam records for the above exon-usage-level intact transcripts (or say candidate intact transcripts)."
sort -k 1,1 ${OutputPrefix}.Assembly.annotate_intact.txt > ${OutputPrefix}.Assembly.annotate_intact_sorted.txt
SearchLine ${OutputPrefix}.Assembly.annotate_intact_sorted.txt 1 ${OutputPrefix}.LRs.ccs.declip.sam 1 > ${OutputPrefix}.LRs.ccs.declip.intact.sam

# Create the augmented CAGE list by incorporating the 5' ends in hg38 annotation.
echo "   Outputting a combined 5ends: "${OutputPrefix}.CAGE.WithRefSeq.5ends.txt
awk '{if($4=="+") print $3"\t"$4"\t"$5; else $3"\t"$4"\t"$6}' ${genome}_refFlat_modified.gpd > ${genome}_5ends.txt
cat ${OutputPrefix}.CAGE.list.sorted.txt ${genome}_5ends.txt | sort -u | sort -k 1,1 -k 3,3n > ${OutputPrefix}.CAGE.WithRefSeq.5ends_sorted.txt


echo "Step 06, Run SupportByCAGEOutputPeak, and record 5' and 3' ends"
# Align the 5' ends of the candidate intact transcripts to CAGE peaks, output the coordinates of 5' and 3'ends of each candidate intact transcript. If the 5' end is not supported by CAGE peak, then "*" marks are output as place holder. Hanging ends within reads are considered in this step.
SupportByCAGEOutputPeak ${OutputPrefix}.CAGE.WithRefSeq.5ends_sorted.txt ${OutputPrefix}.LRs.ccs.declip.intact.sam | sort -k 1,1 > ${OutputPrefix}.Intact.TSS_polyA_ends.txt

# Record the gene name, isoform name, orientation, and end coordinates of each candidate intact transcript.
paste ${OutputPrefix}.Assembly.annotate_intact_sorted.txt ${OutputPrefix}.Intact.TSS_polyA_ends.txt | cut -f1,2,3,5,10,11,12,13 > ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read.txt

# Sort the records for 3' ends polish. Sort according to gene names first. Within the records for each gene, sort the 3' coordinates.
sort -k 2,2 -k 8,8n ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read.txt > ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_sorted_3end.txt

echo "Step 07, Polish 3ends, and then 5ends, then record all possible ends"
echo "   Output: "${OutputPrefix}.AllPossibleEnds_isoform.txt
# Polish the 3' ends.
PolishEndsPolyAUpdate ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_sorted_3end.txt > ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_polished_3end.txt

# Sort the records for 5' ends polish. Sort according to gene names first. Within the records for each gene, sort the 5' coordinates.
sort -k 2,2 -k 6,6n ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_polished_3end.txt > ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_sorted_5end.txt

# Polish the 5' ends.
PolishEndsTSSUpdate ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_sorted_5end.txt > ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_polished.txt

# For each isoform, record all the possibilities of 5' end and 3' end coordinates. TEMP_ENDS.txt records all the known subisoforms.
awk '{if($5!="*") print}' ${OutputPrefix}.Isoform.TSS_polyA_ends_per-read_polished.txt | cut -f3,4,6,8 | sort -u > ${OutputPrefix}.AllPossibleEnds_isoform.txt

echo "Step 08, Modify assembly GPD file, and get polished ends assembly"
ModifyGPDFile ${OutputPrefix}.Assembly_sorted.gpd ${OutputPrefix}.AllPossibleEnds_isoform.txt > ${OutputPrefix}.Assembly_modified_isoform.gpd
cut -f1 ${OutputPrefix}.AllPossibleEnds_isoform.txt | sort -u > ${OutputPrefix}.AllPossibleEnds_need_modification.txt
CounterSearchLine ${OutputPrefix}.AllPossibleEnds_need_modification.txt 1 ${OutputPrefix}.Assembly_sorted.gpd 2 > ${OutputPrefix}.Assembly_unmodified_isoform.gpd
cat ${OutputPrefix}.Assembly_modified_isoform.gpd ${OutputPrefix}.Assembly_unmodified_isoform.gpd > ${OutputPrefix}.Assembly.EndsPolish.gpd
sort -k 3,3 -k 5,5n ${OutputPrefix}.Assembly.EndsPolish.gpd > ${OutputPrefix}.Assembly.EndsPolish_sorted.gpd

echo "Step 09, Further adjustments on assembly"
# Further adjustments on assembly, including the correction of 5-6th colums and redundancy deletion.
CorrectGPD56Column ${OutputPrefix}.Assembly.EndsPolish_sorted.gpd > ${OutputPrefix}.Assembly.EndsPolish_corrected.gpd
DetectRedundant ${OutputPrefix}.Assembly.EndsPolish_corrected.gpd > ${OutputPrefix}.Assembly.EndsPolish.redundant_list.txt
DeleteRedundancyInGPD ${OutputPrefix}.Assembly.EndsPolish_corrected.gpd ${OutputPrefix}.Assembly.EndsPolish.redundant_list.txt > ${OutputPrefix}.Assembly.CAGECorrected.gpd
sort -k 3,3 -k 5,5n ${OutputPrefix}.Assembly.CAGECorrected.gpd |awk '{OFS="\t";print $1,$2,$3,$4,$5,$6,$5,$6,$9,$10,$11}' > ${OutputPrefix}.Assembly.CAGECorrected.Sorted.gpd

echo "Step 10, Finishing"
# Modify the GPD file for long reads to add the clip lengths at ends.
AddClipToGPD ${OutputPrefix}.LRs.ccs.declip.gpd ${OutputPrefix}.LRs.ccs.declip.sam > ${OutputPrefix}.LRs.ccs.correct.clip.gpd

# Align the reads to assembly, find out the final intact transcript list.
GPD_annotate_stringent ${OutputPrefix}.Assembly.CAGECorrected.Sorted.gpd ${OutputPrefix}.LRs.ccs.correct.clip.gpd > ${OutputPrefix}.Assembly.CAGECorrected.IntactReadsSubisoforms.txt

echo "*************************************************************************************************"

fi