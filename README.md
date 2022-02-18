# PacBioWorkflow
A collection of scripts for our PacBio paper:

Sun, Y.H., Wang, A., Song, C., Shankar, G., Srivastava, R.K., Au, K.F. and Li, X.Z. Single-molecule long-read sequencing reveals a conserved intact long RNA profile in sperm. *Nature Communications* 12, 1361 (2021)., [https://doi.org/10.1038/s41467-021-21524-6](https://www.nature.com/articles/s41467-021-21524-6)

This workflow contains C++ codes and pipeline scripts, for mouse and human PacBio assembly.

If CAGE and PAS data supported, run all the three workflows for intact isoform assembly with end correction. If not, only need to run the first workflow to define intact isoforms at exon-usage level.

![](/images/Workflow.png)

## Before running the workflow

1, Download codes

```
git clone https://github.com/LiLabZhaohua/PacBioWorkflow.git
```

1, Compile C++ scripts

The compilation of C++ codes requires g++ and [boost library](https://www.boost.org/).

```
chmod +x Makefile.sh
Makefile.sh
```

Add compiled file directory to $PATH

If you are using Linux system, all compiled binary files are now in the release, and you can download them here:

https://github.com/LiLabZhaohua/PacBioWorkflow/archive/refs/tags/v1.0.0.tar.gz

2, Test the pipelines

```
chmod +x PacBioWorkflow/*sh
```

Add PacBioWorkflow directory to $PATH

Then type pipeline names to get manual page (these three files need to be run in the following order):

```
PacbioGetAssembly.sh
PacbioManualFilterByCAGE.sh
PacbioPASFiltering.sh
```

## Prepare input files and run the workflow

All required files for mouse workflow can be found in ./PacBioData_Mouse except the long reads (LRs) [bam file](https://github.com/LiLabZhaohua/PacBioWorkflow/releases/download/v1.0/PacBioSpermTestisMerged_sgl.bam) was uploaded in the release separately.

All required files for mouse workflow can be found in ./PacBioData_Human, except the long reads (LRs) [bam file](https://github.com/LiLabZhaohua/PacBioWorkflow/releases/download/v1.0/sort_flnc_ccs.bam) was uploaded in the release separately.

The follwoing three workflows need to be run in the same directory (input files should also be in the same directory).

Here, we use mouse mm10 as an example. For human, change the InputParameter files accordingly.

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

## Pipeline output files:

### Exon usage level intact isoform identification

If no CAGE and PAS data provided, only need to run the first workflow: PacbioGetAssembly.sh

1, *.IntactIsoform_sorted.gpd

This is the final assembly of exon usage level isoforms in GenePred (gpd) format.

2, *.IntactIsoform_sorted.ReadsMapped.txt

This file contains the mapping of LRs (intact + decay + unmapped) to intact assembly (*.IntactIsoform_sorted.gpd)

### Intact isoform identification with CAGE and PAS data

1, *.Assembly.CAGECorrected.PASCorrected.Final.Sorted.gpd

This file contains final assembly in GenePred (gpd) format, supported by LRs, CAGE and PAS peaks or RefSeq ends.

2, *.Assembly.CAGECorrected.PASCorrected.IntactReadsSubisoforms.txt

This file contains the mapping of intact LRs to assembly.

