# BloomGenome

This program is aiming to apply bloom filter to genomic sequence processing and kinship identification, currently under development.

## Dependencies
The program requires following dependencies. If you use anaconda/miniconda, you can run the following to initialize the environment.

```bash
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment
conda create --name bloom-env

# activate conda environment
conda activate bloom-env

conda install -c bioconda -c conda-forge -c anaconda python=3 bcftools numpy bitarray mmh3
```
If you face any issues on installing any dependencies via conda, try to install them directly and input their execute path into the program manually.

Otherwise, please manually install the following dependencies.
* bcftools
* numpy
* bitarray
* mmh3