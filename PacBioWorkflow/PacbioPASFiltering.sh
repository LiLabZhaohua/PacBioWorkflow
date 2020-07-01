#!/usr/bin/env bash
#PacbioPASFiltering.sh
#This pipeline performs PAS correction after CAGE correction
#This script only reports isoforms that can be supported by PAS peaks or RefSeq
#Modified from Anqi Wang's workflow
#Version: Yu Sun, 2020-05-11, updated on 2020-06-30

if [ ! -n "$2" ]
then
  echo "*************************************************************************************************"
  echo "*                        Welcome to use PacbioPASFiltering pipeline                             *"
  echo "*This pipeline performs PAS correction after CAGE correction                                    *"
  echo "*This script only reports isoforms that can be supported by PAS peaks or RefSeq                 *"
  echo "*This is the 2020 version, and please run after running PacbioManualFilterByCAGE.sh             *"
  echo "*This workflow requires several critical files to run:                                          *"
  echo "*1, Genome annotation: \${genome}_refFlat_modified.gpd                                           *"
  echo "*2, PAS peaks: \${PAS_Peaks}, or you can also use PAS sam file: \${PAS_SAM}                       *"
  echo "*3, Update all names in InputParameters.txt file, and use that as input to this pipeline        *"
  echo "*Critical information in the InputParameters.txt file:                                          *"
  echo "*        genome=\"hg38\"                                                                          *"
  echo "*        PAS_SAM=\"file_selected_sgl.sam\"                                                        *"
  echo "*        #PAS_Peaks=\"PAS.peaks\"                                                                 *"
  echo "*In addition, this workflow also reads the previous output from PacbioGetAssembly.sh:           *"
  echo "*   \${OutputPrefix}.Assembly.CAGECorrected.Sorted.gpd                                           *"
  echo "*   \${OutputPrefix}.Assembly.CAGECorrected.IntactReadsSubisoforms.txt                           *"
  echo "*Usage: `basename $0` [InputParameters] [OutputPrefix]                                  *"
  echo "*Output: 1, *.Assembly.CAGECorrected.PASCorrected.Sorted.gpd                                    *"
  echo "*        2, *.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.txt                    *"
  echo "*************************************************************************************************"
else

PAS_SAM="NA"
PAS_Peaks="NA"
module load boost
module load kentutils
source $1
OutputPrefix=$2

echo "*************************************************************************************************"
echo "*               Welcome to use PacbioPASFiltering single PacBio data pipeline                   *"
echo "Starting the pipeline with following parameters:"
echo "Genome annotation: "${genome}_refFlat.txt" or "${genome}_refFlat_modified.gpd
echo "PAS peaks, one file exists will be enough:"
echo "  PAS peaks sam: "${PAS_SAM}
echo "  PAS peaks file: "${PAS_Peaks}
echo "OutputPrefix: "${OutputPrefix}

echo "Step 01, Getting PAS peaks"
if [[ $PAS_Peaks == "NA" ]];then
  echo "  No PAS_Peaks, generating"
  if [[ ! $PAS_SAM == "NA" ]];then
    PASReadEnd $PAS_SAM  | sed '/random/d' | sort -k 1,1 -k 3,3n | uniq > ${OutputPrefix}.PAS_Peaks.txt
    PAS_Peaks=${OutputPrefix}.PAS_Peaks.txt
    echo "  PAS_Peaks generated"
  else
  	echo "  PAS information missing... exit"
  	exit 1
  fi
else
  echo "  PAS_Peaks found"
fi

echo "Step 02, Create a combination of PAC_Peaks and RefSeq 3ends"
Annotated3PrimeEnd ${genome}_refFlat_modified.gpd | sort -u > ${OutputPrefix}.PAS.annotated_3end.txt
cat ${PAS_Peaks} ${OutputPrefix}.PAS.annotated_3end.txt | sort -u > ${OutputPrefix}.PAS.RefSeq.3end.txt

echo "Step 03, Getting filtered transcripts by PAS peaks"
PASFilterGPDByPeaks ${OutputPrefix}.PAS.RefSeq.3end.txt ${OutputPrefix}.Assembly.CAGECorrected.Sorted.gpd > ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Temp.gpd
sort -k 3,3 -k 5,5n ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Temp.gpd > ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Sorted.gpd

SearchLine ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Sorted.gpd 2 ${OutputPrefix}.Assembly.CAGECorrected.IntactReadsSubisoforms.txt 3 > ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.txt
awk '{print $3}' ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.txt|sort|uniq > ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.list
SearchLine ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.list 1 ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Sorted.gpd 2 | sort -k 3,3 -k 5,5n > ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Final.Sorted.gpd
awk '{print $3"-"$2}' ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.txt|sort|uniq -c|awk '{OFS="\t";print $2,$1}' > ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.count.txt
echo "Final intact reads count:"
wc -l ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.count.txt
echo "Final isoform number:"
wc -l ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Final.Sorted.gpd

echo "Step 04, Final"
awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Final.Sorted.gpd > ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Final.Sorted.gpd.cut
genePredToGtf file ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Sorted.gpd.cut ${OutputPrefix}.Assembly.CAGECorrected.PASCorrected.Final.Sorted.gtf

fi
