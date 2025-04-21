# Ervrmancer

ERVmancer is a bioinformatics conda package that quantifies Human Endogenous Retrovirus (HERV) short read RNA sequencing expression data by aligning short reads to a curated subset of HERVs and then resolving ambiguity in alignment using a pre-computed HERV phylogenetic tree.

## Installation and Usage

```bash
conda install -c bioconda ervmancer
```

## Maintainer Directions
### Setting up Your Local Environment (Linux)

* Install the following:
    * Anaconda
    * conda-build
    * git

From the Conda Docs, the most straightforward way to do this is to first install Anaconda, then use conda through CLI to install conda-build and Git.

```
conda create -n <environment_name> python=3.8
conda install -n <environment_name> conda-build git conda-verify
```

The bowtie2 indices used 

### Dependencies

The following are dependencies for ERVmancer (and are also in bioconda)
* bowtie2 >= 2.42
* samtools >= 1.2
* bedtools >= 2.29.2

## Troubleshooting

See internal cluster docs or Notion for more information. For default testing of ervmancer on your local machine, see below:

To set up a fresh environment to build/test ervmancer locally:

```bash
conda create -n ervmancer_local
conda install -n ervmancer_local
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

To iteratively test ervmancer as a local package and clear your cache:

```bash
# Purge previous local copies
conda build purge
# Build locally and install
conda build . -c bioconda
conda install -y <path_to_local_ervmancer_build>
```

Usually, if a versioning conflict happens when installing your local test build please look for a similar excerpt to below, near the ending of the build standard output:
```
# To have conda build upload to anaconda.org automatically, use
# conda config --set anaconda_upload yes
anaconda upload \
    /home/<user>/conda-bld/noarch/ervmancer-0.0.x-py_x.tar.bz2
anaconda_upload is not set.  Not uploading wheels: []
```

The path after anaconda upload, should be the manual install path to the local ervmancer package build.

## Authors
* Andrew Patterson
* Bryant Duong