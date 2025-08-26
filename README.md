# ERVmancer

[![Anaconda-Server Badge](https://anaconda.org/bioconda/ervmancer/badges/version.svg)](https://anaconda.org/bioconda/ervmancer) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ervmancer/badges/platforms.svg)](https://anaconda.org/bioconda/ervmancer)

ERVmancer is a bioconda package that quantifies Human Endogenous Retrovirus (HERV) short read RNA sequencing expression data by aligning short reads to a curated subset of HERVs and then resolving ambiguity in alignment using a pre-computed HERV phylogenetic tree.

## Paper

Andrew Patterson, Bryant Duong, Leena Yoon, Maya Foster, Lauren MacMullen, Samantha Soldan, Jayamanna Wickramasinghe, Anastasia Lucas, Avi Srivastava, Maureen Murphy, Paul Lieberman, Noam Auslander. ERVmancer: A phylogenetic approach to quantification of human endogenous retrovirus expression. In Preparation.

## User Installation and Usage

In your desired conda virtual environment, install using the below commands.
```bash
# necessary for ervmancer and its dependencies, python >= 3.8
conda create --name ervmancer_env python=3.8
conda activate ervmancer_env
conda config --add channels bioconda
conda config --add channels conda-forge
conda install ervmancer
```

Download ```clean_kmer_31_60_percent_cutoff.pkl``` and ```GRCh38_noalt_as.tar.gz``` from [Zenodo](https://zenodo.org/records/15231904). Move these files from your download directory to your desired directory and extract the GRCh38 index.

```tar -xzf <path to tar.gz> -C <desired extracted folder output path>```

### Usage Options

See below documentation for examples on how to run the possible parameters/entrypoints:

* [Entrypoint 1 (from fastq file) and Entrypoint 2 (from aligned SAM file)](https://github.com/AuslanderLab/ERVmancer/wiki/Primary-Entrypoint-Guide)

* [Advanced: Entrypoint 3 (mapping counts to phylogenetic tree)](https://github.com/AuslanderLab/ERVmancer/wiki/Advanced-Entrypoint-Guide)


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

### Dependencies

The following are dependencies for ERVmancer (and are also in the bioconda/conda-forge channel)
* bowtie2 >= 2.42
* samtools >= 1.2
* bedtools >= 2.29.2

The metadata used in this package is hosted on [Zenodo](https://zenodo.org/records/15231904) and can be downloaded here.

### Releasing Versions

This project uses [SemVer](https://semver.org/) for versioning. In order to create a release package, click [create a new release](https://github.com/AuslanderLab/ERVmancer/releases/new) and create a new tag with the corresponding bumped version ```(vX.y.z)```. Once a release has been published, update ervmancer on [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/ervmancer) using your own remote fork with the proper versioning and sha256 hash. Then, create a PR in bioconda-recipes to send your changes for review with the bioconda team.

To grab the sha256 hash of the new release, download the tar file from the ervmancer releases page and use the following command:

```wget -O- https://github.com/AuslanderLab/ERVmancer/archive/refs/tags/v<version>.tar.gz | shasum -a 256```

## Conda Troubleshooting

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
conda uninstall -y ervmancer
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

## Package Authors
* Andrew Patterson
* Bryant Duong
