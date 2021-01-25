# **_Bloom_: A computational framework to reveal occult patterns in 3C data**

### Bloom is a computational framework which takes a 3C-like contact map as input and will produce:
* A new contact matrix showing potential occult patterns
* A list of putative loops with an intrinsic scoring method (IFS; Interaction Frequency Score); shown to highgly correlate with function

### Bloom is a very powerful framework, especially - but not limited to - the following types of analyses:
* Sparse Data: Bloom is able to reveal patterns in sparse data (e.g. from single-cell Hi-C, Dip-seq, etc.)
* Enhancer Function: Bloom has been shown to be able to correlate its results with elements' impact on gene expression. Allowing the creation of a cell's full enhancer atlas

### Bloom allows the following formats as input and output:
* Juicer's ".hic" (single or multiple resolutions)
* Cooler's ".(m)cool" (single or multiple resolutions)
* Bedgraph ".bg2" plain text (single or multiple resolutions)

---

## Requirements

Bloom has very few requirements, most of which are usually used by any standard scientific python environment:
* Numpy (v. 1.17.x or higher): https://numpy.org/
* Scipy (v. 1.4.x or higher): https://www.scipy.org/
* Setuptools (python installation): https://pypi.org/project/setuptools/

The installation of all the requirements can be easily done with python package manager PyPi:

```
pip install -U setuptools
python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
```

---

## Installation

After making sure all the dependencies are installed. Simply execute the following steps:

```
git clone https://github.com/CostaLab/reg-gen.git
cd reg-gen
python setup.py install --user
```

This installation will automatically enable all C/C++ modules given your PC requirements and create a folder under your $HOME directory called _bloom_data_ containing extra data and scripts needed to obtain larger data, required to execute Bloom (check step 2 below). A more thorough tutorial will be made available soon.

:warning: Bloom has been tested and is supported only for python 3. Bloom is **not** supported on Windows systems.

---

## Simple Usage

Here we show a simple step-by-step example to execute Bloom. A complete tutorial and the full documentation of Bloom's package will be available soon.

### 1. Downloading barcodes

TBD

### 2. Downloading a sample dataset (from the original _Bloom_ manuscript)

TBD

### 3. Execution of Bloom

TBD

---

## Help and Citation

For usage help, comments, suggestions and error reports, please send an email to:
eduardo.gadegusmao at med.uni-goettingen.de

Bloom's original manuscript is still under editorial process. If you wish to cite this tool, or the findings on the paper, please refer to the pre-print on BioRxiv: https://www.biorxiv.org/content/10.1101/2020.11.10.376533v1.full

Textual citation:

Gade Gusmao E, Mizi A, Brant L, Papantonis A. Retrieving high-resolution chromatin interactions and decoding enhancer regulatory potential in silico doi: https://doi.org/10.1101/2020.11.10.376533

BibTex:

```
@article {Gusmao2021,
	author = {Gusmao, Eduardo Gade and Mizi, Athanasia and Brant, Lilija and Papantonis, Argyris},
	title = {Retrieving high-resolution chromatin interactions and decoding enhancer regulatory potential in silico},
	elocation-id = {2020.11.10.376533},
	year = {2021},
	doi = {10.1101/2020.11.10.376533},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2020/11/10/2020.11.10.376533},
	eprint = {https://www.biorxiv.org/content/early/2020/11/10/2020.11.10.376533.full.pdf},
	journal = {bioRxiv}
}
```

---
