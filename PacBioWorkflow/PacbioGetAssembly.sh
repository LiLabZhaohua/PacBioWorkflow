#!/usr/bin/env bash
#PacbioGetAssembly.sh
#This pipeline gets assembly from the PacBio LRs and SRs assembly, and use RefSeq as annotation
#This script only reports LR supported reads, and won't do CAGE or PAS correction
#Modified from Anqi Wang's workflow
#Version: Yu Sun, 2020-05-01
#Version: Yu Sun, 2020-05-20, remove EndsCompare, AntisenseCandidate
#Version: Yu Sun, 2020-07-01, update MapReadToIsoform

if [ ! -n "$2" ]
then
  echo "*************************************************************************************************"
  echo "*                          Welcome to use PacbioGetAssembly pipeline                            *"
  echo "*This pipeline gets assembly from the PacBio LRs and SRs assembly, and use RefSeq as annotation *"
  echo "*This script only reports LR supported reads and isoforms, but won't do CAGE or PAS correction  *"
  echo "*This is the 2020 version of this pipeline, modified from Anqi's workflow                       *"
  echo "*This workflow requires several critical files to run:                                          *"
  echo "*1, Genome annotation: \${genome}_refFlat.txt or \${genome}_refFlat_modified.gpd                  *"
  echo "*2, LRs bam/bai files: \${LR_Bam}                                                                *"
  echo "*3, StingTie SR assembly: \${StringTie_Assembly}                                                 *"
  echo "*4, Update all names in InputParameters.txt file, and use that as input to this pipeline        *"
  echo "*Critical information in the InputParameters.txt file:                                          *"
  echo "*        genome=\"hg38\"     \${genome}_refFlat.txt or \${genome}_refFlat_modified.gpd              *"
  echo "*        LR_Bam=\"sort_flnc_ccs.bam\"                                                             *"
  echo "*        StringTie_Assembly=\"stringtie_hm.gpd\"                                                  *"
  echo "*Usage: `basename $0` [InputParameters] [OutputPrefix]                                   *"
  echo "*Outputs that will be used further:                                                             *"
  echo "*        1, *annotation_intact.txt                                                              *"
  echo "*        2, *LRs.ccs.declip.gpd                                                                 *"
  echo "*        3, *LRs.ccs.declip.sam                                                                 *"
  echo "*        4, *StringTie.filtered.assgn.chtnm.LR.annotate_intact.txt                              *"
  echo "*        5, *StringTie.filtered.assgn.chtnm.gpd                                                 *"
  echo "*Other useful outputs if you are not going to perform CAGE and PAS correction                   *"
  echo "*        6, *IntactIsoform_sorted.gpd (Final assembly at exon usage level)                      *"
  echo "*        7, *IntactIsoform_sorted.ReadsMapped.txt (Mapping of intact+decay reads to assembly)   *"
  echo "*************************************************************************************************"
else

source $1
OutputPrefix=$2
echo "*************************************************************************************************"
echo "*               Welcome to use PacbioGetAssembly single PacBio data pipeline                    *"
echo "Starting the pipeline with following parameters:"
echo "Genome annotation: "${genome}_refFlat.txt" or "${genome}_refFlat_modified.gpd
echo "LRs bam file: "${LR_Bam}
echo "StringTie gpd file: "${StringTie_Assembly}
echo "OutputPrefix: "${OutputPrefix}

module load samtools
module load boost

echo "Step 01, Modify annotation refFlat file"
if [[ ! -s ${genome}_refFlat_modified.gpd ]];then
  echo "  No modified refFlat, generating"
  if [ -s ${genome}_refFlat.txt ];then
    ModifyAnnotationGene ${genome}_refFlat.txt > ${genome}_refFlat_temp
    ModifyAnnotationTranscript ${genome}_refFlat_temp > ${genome}_refFlat_modified.gpd
    echo "  Modified refFlat generated"
  fi
else
  echo "  Modified refFlat found"
fi
sort -k3,3 -k5,5n ${genome}_refFlat_modified.gpd > ${genome}_refFlat_modified_sortpos.gpd
sort -k 1,1 ${genome}_refFlat_modified.gpd > ${genome}_refFlat_modified_sortname.gpd

echo "Step 02, Filter clipping and correct LRs"
echo "  This step outputs: "${OutputPrefix}.annotation_intact.txt", "${OutputPrefix}.annotation_nonintact.gpd
samtools view ${LR_Bam} > ${OutputPrefix}.LRs.ccs.sam
ClipFilter ${OutputPrefix}.LRs.ccs.sam | awk '{print $1}' > ${OutputPrefix}.LRs.ccs.short_clip
SAMtoGPD ${OutputPrefix}.LRs.ccs.sam > ${OutputPrefix}.LRs.ccs.gpd
SearchLine ${OutputPrefix}.LRs.ccs.short_clip 1 ${OutputPrefix}.LRs.ccs.gpd 1 > ${OutputPrefix}.LRs.ccs.declip.gpd
#GPD_annotate requires sorted GPD file, and unsorted file won't generate output
GPD_annotate ${genome}_refFlat_modified_sortpos.gpd ${OutputPrefix}.LRs.ccs.declip.gpd ${OutputPrefix}.annotation

echo "Step 03, Filter StringTie results"
#stringtie ../alignment/sort_hsa_sperm_cutadapt.bam -o stringtie_hm.gtf -p 10 -A stringtie_hm_abundance.txt
#gtf_to_genepred.py stringtie_hm.gtf > ${StringTie_Assembly}
#FilterStringtieResult uses qsort, which is not working on bluehive. So use sorted refFlat
FilterStringtieResult ${genome}_refFlat_modified_sortpos.gpd ${StringTie_Assembly} ${OutputPrefix}.StringTie.filtered.gpd

echo "Step 04, Assign novel gene and transcript names to StringTie assembly"
#Here we use annotation sorted based on 
AssignGene ${genome}_refFlat_modified_sortname.gpd ${OutputPrefix}.StringTie.filtered.gpd > ${OutputPrefix}.StringTie.filtered.assgn.gpd
AssignTranscript ${genome}_refFlat_modified_sortname.gpd ${OutputPrefix}.StringTie.filtered.assgn.gpd ${OutputPrefix}.StringTie.filtered.assgn.chtnm.temp.gpd
#This step is modified from Anqi, we have to sort by location, rather than names
sort -k 3,3 -k5,5n ${OutputPrefix}.StringTie.filtered.assgn.chtnm.temp.gpd | uniq > ${OutputPrefix}.StringTie.filtered.assgn.chtnm.gpd
rm -rf ${OutputPrefix}.StringTie.filtered.assgn.gpd
rm -rf ${OutputPrefix}.StringTie.filtered.assgn.chtnm.temp.gpd

echo "Step 05, Rescue non-intact transcripts using StingTie SR assembly, and get Intact isoforms, Intact reads"
GPD_annotate ${OutputPrefix}.StringTie.filtered.assgn.chtnm.gpd ${OutputPrefix}.annotation_nonintact.gpd ${OutputPrefix}.StringTie.filtered.assgn.chtnm.LR.annotate
cat ${OutputPrefix}.annotation_intact.txt ${OutputPrefix}.StringTie.filtered.assgn.chtnm.LR.annotate_intact.txt > ${OutputPrefix}.IntactReadsSplicingLevel.txt
cat ${genome}_refFlat_modified.gpd ${OutputPrefix}.StringTie.filtered.assgn.chtnm.gpd | sort -u > ${OutputPrefix}.MergedIsoform.gpd
awk '{print $3}' ${OutputPrefix}.IntactReadsSplicingLevel.txt | sort -u > ${OutputPrefix}.IntactIsoform.txt
SearchLine ${OutputPrefix}.IntactIsoform.txt 1 ${OutputPrefix}.MergedIsoform.gpd 2 > ${OutputPrefix}.IntactIsoform.gpd
sort -k 3,3 -k 5,5n ${OutputPrefix}.IntactIsoform.gpd > ${OutputPrefix}.IntactIsoform_sorted.gpd
SearchLine ${OutputPrefix}.IntactReadsSplicingLevel.txt 1 ${OutputPrefix}.LRs.ccs.declip.gpd 1 > ${OutputPrefix}.IntactReadsSplicingLevel.gpd
echo "  Intact isoforms assembly, splicing level: "${OutputPrefix}.IntactIsoform_sorted.gpd
echo "  Isoform counts: "
wc -l ${OutputPrefix}.IntactIsoform_sorted.gpd

#echo "Step 06, Get the ends from Intact isoforms and Intact reads: EndsIntact files"
#awk '{if($4=="+") print $1"\t"$5"\t"$6; else print $1"\t"$6"\t"$5}' ${OutputPrefix}.IntactReadsSplicingLevel.gpd > ${OutputPrefix}.EndsIntactRead.txt
#awk '{if($4=="+") print $2"\t"$5"\t"$6; else print $2"\t"$6"\t"$5}' ${OutputPrefix}.IntactIsoform.gpd > ${OutputPrefix}.EndsIntactIsoform.txt

#sort -k 1,1 ${OutputPrefix}.EndsIntactRead.txt > aaa
#sort -k 1,1 ${OutputPrefix}.IntactReadsSplicingLevel.txt > bbb
#paste aaa bbb | head
#paste aaa bbb | awk '{print $1"\t"$6"\t"$2"\t"$3}' > ${OutputPrefix}.EndsIntactRead_new.txt
#rm -rf ${OutputPrefix}.EndsIntactRead.txt
#mv ${OutputPrefix}.EndsIntactRead_new.txt ${OutputPrefix}.EndsIntactRead.txt

#EndsCompare script is not important, script missing
#EndsCompare ${OutputPrefix}.EndsIntactIsoform.txt ${OutputPrefix}.EndsIntactRead.txt | awk '{if($3<100 && $4<30) print}' > ${OutputPrefix}.EndsDiff_filter.txt
#awk '{print $2}' ${OutputPrefix}.EndsDiff_filter.txt | sort -u > aaa
#SearchLine aaa 1 ${OutputPrefix}.IntactIsoform.gpd 2 > ${OutputPrefix}.IntactIsoform_filter.gpd
#rm -rf aaa

#This StatQuality is missing
#StatQuality $LR_Bam 15 > ${OutputPrefix}.LR.bam.stat_quality.txt

echo "Step 06, Map reads to assembly, get Compatible reads (both intact and decay), and antisense reads"
echo "  Mapping of reads to intact isoforms (splicing level): "${OutputPrefix}.IntactIsoform_sorted.ReadsMapped.txt
#MapReadToIsoform: This is a critical script to map LRs to the assembly GPD, and get intact and decay reads
MapReadToIsoform ${OutputPrefix}.IntactIsoform_sorted.gpd ${OutputPrefix}.LRs.ccs.declip.gpd > ${OutputPrefix}.IntactIsoform_sorted.ReadsMapped.txt
cat ${OutputPrefix}.IntactIsoform_sorted.ReadsMapped.txt | grep Compatible | awk '{print $1}' | sort -u | wc -l

#AntisenseCandidate ${genome}_refFlat_modified.gpd ${OutputPrefix}.LRs.ccs.declip.gpd > ${OutputPrefix}.antisense.txt
#awk '{if($11>200 || $12>0.1) print}' ${OutputPrefix}.antisense.txt | wc -l

echo "Step 07, Reads, Species calculation"
SearchLine ${OutputPrefix}.LRs.ccs.declip.gpd 1 ${OutputPrefix}.LRs.ccs.sam 1 > ${OutputPrefix}.LRs.ccs.declip.sam
ClusterUniqueReads ${OutputPrefix}.LRs.ccs.declip.gpd ${OutputPrefix}.LRs.ccs.declip.sam > ${OutputPrefix}.read_species.txt
MapReadToIsoform ${OutputPrefix}.IntactIsoform_sorted.gpd ${OutputPrefix}.LRs.ccs.declip.gpd | grep Compatible | awk '{print $1}' | sort -u > aaa
echo "  Count lines, reads:"
SearchLine aaa 1 ${OutputPrefix}.read_species.txt 1 | wc -l
rm -rf aaa

echo "  Count lines, species:"
SearchLine ${OutputPrefix}.IntactReadsSplicingLevel.gpd 1 ${OutputPrefix}.read_species.txt 1 | wc -l


echo "*************************************************************************************************"

fi