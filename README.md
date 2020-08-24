# PacBioWorkflow
A collection of scripts for our PacBio paper:

Sun, Y.H., Wang, A., Song, C., Srivastava, R.K., Au, K.F., Li, X.Z. (2020) Single-molecule long-read sequencing reveals a conserved selection mechanism retaining intact long RNAs and miRNAs in sperm. bioRxiv, [doi:10.1101/2020.05.28.122382](https://www.biorxiv.org/content/10.1101/2020.05.28.122382v1)

This workflow contains C++ codes and pipeline scripts, for mouse and human PacBio assembly.

If CAGE and PAS data supported, run all the three workflows for intact isoform assembly with ends correction. If not, only need to run the first workflow to define intact isoforms at exon-usage level.

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

All required files for mouse workflow can be found in ./PacBioData_Mouse, except the long reads (LRs) [bam file]()

All required files for mouse workflow can be found in ./PacBioData_Human, except the long reads (LRs) [bam file]()

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

