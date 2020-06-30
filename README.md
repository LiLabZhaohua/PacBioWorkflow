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

## Input files for the pipelines

### 1, PacbioGetAssembly.sh

This pipeline requires RefSeq annotation in refFlat format, long reads (LRs) bam file (with bam index file), and short reads (SRs) assembly by StringTie.

```
#Write the input file names in InputParameters1.txt file
cat InputParameters1.txt
  #Required parameters for PacbioGetAssembly.sh
  genome="mm10"                           #Files: ${genome}_refFlat.txt, if not, then ${genome}_refFlat_modified.gpd is also acceptable
  LR_Bam="PacBioSpermTestisMerged_sgl.bam"
  StringTie_Assembly="stringtie_merge.gpd"

#Run the pipeline:
PacbioGetAssembly.sh InputParameters1.txt MouseCombined
```




