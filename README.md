# **_Livros_:**

:warning: Eu vou colocar os links para os livros aqui temporariamente, caso você venha no futuro e os links não mais estejam aqui, basta falar comigo.

---

### Os links para cada 'pasta', contêm os seguintes livros: 

```bash
Livros
├── livros_algoritmos_e_Programacao
│             ├── Cormen_Introduction_to_Algorithms.pdf
│             ├── Culhane_Introduction_to_Programming_in_R.pdf
│             ├── Halsvorsen_Python_Programming.pdf
│             ├── Kernighan_Ritchie_EEE_C_Programming_Language.pdf
│             ├── Kochan_Programming_in_C.pdf
│             ├── Paradis_R_for_Beginners.pdf
│             └── Udayan_Das_Introduction_to_Python_Programming.pdf
├── livros_aprendizagem_de_maquina
│             ├── An_Introduction_to_Statistical_Learning_with_Applications_in_Python_James.pdf
│             ├── Artificial_Intelligence_A_Modern_Approach_Russel.pdf
│             ├── Deep_Learning_Goodfellow.pdf
│             ├── Deep_Learning_for_Computer_Vision_Rosebrock.pdf
│             ├── Deep_Learning_with_Python_Chollet.pdf
│             ├── Durr_Probabilistic_Deep_Learning.pdf
│             ├── Livros_Aprendizado_de_Maquina.txt
│             ├── Pattern_Classification_Duda_Hart_Stork.pdf
│             ├── Pattern_Recognition_and_Machine_Learning_Bishop.pdf
│             ├── The_Hundred_Page_Machine_learning_Book_Burkov.pdf
│             └── livro_eiben.pdf
├── livros_bioinformatica
│             ├── 2012_Stuart_Brown_Next_Generation_DNA_Sequencing_Informatics.pdf
│             ├── 2023_Hamid_D_ Ismail_Bioinformatics_ A_Practical_Guide to_Next_Generation.pdf
│             ├── Bioinformatics - Sequence and Genome Analysis, Second Edition by David Mount.pdf
│             ├── Biotech_Consortium_India_Bioentrepreneurship_Development.pdf
│             ├── DW_Mount_BioinformaticsSequenceandGenomeAnalysis.pdf
│             ├── Deep_Learning_with_Python.pdf
│             ├── Janeway's Immunobiology, 8th Edition - Kenneth M. Murphy.pdf
│             ├── Koche_Fundamentos_de_Metodologia_Científica.pdf
│             ├── Low_Tammi_Bioinformatics_A_Practical_Handbook_of_Next_Generation_Sequencing_and_Its_Applications.pdf
│             ├── Marconi_Fundamentos_de_metodologia_científica.pdf
│             ├── Pevsner2011.pdf
│             ├── Practical Guide to ChIP-seq Data Analysis.pdf
│             ├── RNAseq data analysis a practical approach Eija Korpenlainen.pdf
│             ├── [Vince_Buffalo]_Bioinformatics_Data_Skills_Reprod(b-ok.org).pdf
│             ├── bioinformatics-the-machine-learning-approach-second-edition-pierre-baldi-soren-brunak-1.pdf
│             ├── durbin_book.pdf
│             ├── pevsner.pdf
│             ├── ramezani_arxiv.pdf
│             └── roy2016.pdf
├── livros_biologia
│             ├── Kandel2021.pdf ("Principles of Neural Science" - 1693 págs.)
│             ├── Lodish2021.pdf ("Molecular Cell Biology [9a edição]" - 4610 págs.)
│             └── Murphy2012_Janeways_Immunobiology_8ed.pdf
├── livros_matematica
│             ├── 1_Precalculus_Stewart.pdf
│             ├── 2_Calculus_Thomas.pdf
│             ├── 3_Elementary_Linear_Algebra_Larson.pdf
│             ├── 4_Mathematical Statistics_Wackerly.pdf
│             ├── 5_Elementary_Statistics_Weiss.pdf
│             ├── 6.Vector_Calculus_Colley.pdf
│             ├── Applied_Linear_Regression_Models_Neter.pdf
│             ├── Applied_Linear_Statistical_Models_Kutner.pdf
│             ├── Applied_Multivariate_Statistical_Analysis_Johnson.pdf
│             ├── Applied_Regression_Analysis_and_Other_Multivariable_Methods_Kleinbaum.pdf
│             ├── Livros_Matematica_TheMathSorcerer.txt
│             ├── Methods_of_Multivariate_Analysis_Rencher.pdf
│             ├── Nonparametric_Statistical_Methods_Hollander.pdf
│             └── Nonparametric_Statistics_Theory_and_Methods_Deshpande.pdf
├────────────────────────────────────────────────────────────────────────────────────────────────────
```

### Links:

* [livros_algoritmos_e_Programacao](https://drive.google.com/file/d/12bCRQSRuNuw_GY1dRo38wQej9ZQvgE7t/view?usp=sharing)
* [livros_aprendizagem_de_maquina](https://drive.google.com/file/d/1kMXcdREf6Inpp3SqjYgelH3freVal1dj/view?usp=sharing)
* [livros_bioinformatica](https://drive.google.com/file/d/13s9iiPfqVLTiSOgTZ97kf6XN4qtnWOg5/view?usp=sharing)
* [livros_biologia](https://drive.google.com/file/d/1BPIcIYIiUULjkM0D5wPaT2-LtpElYJHU/view?usp=sharing)
* [livros_matematica](https://drive.google.com/file/d/1QGzU000ojQ4Co0_eE_j9DpiNoHEhnyMl/view?usp=sharing)


---
---
---
---
---
---

# **_Bloom_: A computational framework to reveal occult patterns in 3C data**

:warning: **Bloom is undergoing some changes**: For the time being, the code is not usable. We have noticed all emails that were sent to us and we will reply to all of them as soon as Bloom is in a working condition again. The current changes are being made to acommodate for some suggestions made by colleagues and reviewers.

---

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
git clone https://github.com/eggduzao/Bloom.git
cd Bloom
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
