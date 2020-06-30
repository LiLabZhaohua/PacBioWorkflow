# PacBioWorkflow
A collection of scripts for our PacBio paper:

Sun, Y.H., Wang, A., Song, C., Srivastava, R.K., Au, K.F., Li, X.Z. (2020) Single-molecule long-read sequencing reveals a conserved selection mechanism retaining intact long RNAs and miRNAs in sperm. bioRxiv, [doi:10.1101/2020.05.28.122382](https://www.biorxiv.org/content/10.1101/2020.05.28.122382v1)

This workflow contains C++ codes and pipeline scripts, for mouse and human PacBio assembly.

## Before running the workflow

1, Download codes

```
git clone https://github.com/LiLabZhaohua/PacBioWorkflow.git
```

1, Compile C++ codes

The compilation of C++ codes requires g++ and [boost library](https://www.boost.org/).

Add compiled file directory to $PATH

2, Test the pipelines

```
chmod +x PacBioWorkflow_Mouse/*sh
```

Then type pipeline names to get manual page (these three files need to be run in the following order):

```
PacbioGetAssembly.sh
PacbioManualFilterByCAGE.sh
PacbioPASFiltering.sh
```

## Prepare input files and run the workflow

All required files for mouse workflow can be found in ./PacBioData_Mouse, except the long reads (LRs) [bam file]()

The follwoing three workflows need to be run in the same directory (input files should also be in the same directory).

### 1, PacbioGetAssembly.sh

This pipeline requires RefSeq annotation in refFlat format, long reads (LRs) bam file (with bam index file), and short reads (SRs) assembly by StringTie.

```
#Write the input file names in InputParameters1.txt file
cat InputParameters1.txt
$ #Required parameters for PacbioGetAssembly.sh
$ genome="mm10"                           #Files: ${genome}_refFlat.txt, if not, then ${genome}_refFlat_modified.gpd is also acceptable
$ LR_Bam="PacBioSpermTestisMerged_sgl.bam"
$ StringTie_Assembly="stringtie_merge.gpd"

#Run the pipeline:
PacbioGetAssembly.sh InputParameters1.txt MouseCombined
```

### 2, PacbioManualFilterByCAGE.sh

This pipeline requires RefSeq annotation in refFlat format, CAGE peak data, and outputs from the previous pipeline.

```
#Write the input file names in InputParameters2.txt file
cat InputParameters2.txt
$ #Required parameters for PacbioManualFilterByCAGE.sh
$ genome="mm10"
$ CAGE_Peaks="CAGE_peaks.txt"

#Run the pipeline:
PacbioManualFilterByCAGE.sh InputParameters2.txt MouseCombined
```

### 3, PacbioPASFiltering.sh

This pipeline requires RefSeq annotation in refFlat format, CAGE peak data, and outputs from the previous pipelines.

```
#Write the input file names in InputParameters3.txt file
cat InputParameters3.txt
$ #Required parameters for PacbioPASFiltering.sh
$ genome="mm10"
$ PAS_Peaks="PAS_peaks.txt"
  
#Run the pipeline:
PacbioPASFiltering.sh InputParameters3.txt MouseCombined
```

### Pipeline output files:

1, *.Assembly.CAGECorrected.PASCorrected.Sorted.gpd

This file contains final assembly in GenePred (gpd) format, supported by LRs, CAGE and PAS peaks or RefSeq ends.

2, *.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.txt

This file contains the mapping of intact LRs to assembly.
