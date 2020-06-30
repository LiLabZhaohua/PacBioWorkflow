# PacBioWorkflow
A collection of scripts for our PacBio paper:

Sun, Y.H., Wang, A., Song, C., Srivastava, R.K., Au, K.F., Li, X.Z. (2020) Single-molecule long-read sequencing reveals a conserved selection mechanism retaining intact long RNAs and miRNAs in sperm. bioRxiv, [doi:10.1101/2020.05.28.122382](https://www.biorxiv.org/content/10.1101/2020.05.28.122382v1)

This workflow contains C++ codes and pipeline scripts, for mouse and human PacBio assembly.

## Before running the workflow

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
