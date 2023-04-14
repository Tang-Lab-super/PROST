# PROST: A quantitative pattern recognition framework for spatial transcriptomics 
## Overview
`PROST` is a flexible framework to quantify gene spatial expression patterns and detect spatial tissue domains using spatially resolved transcriptomics with various resolutions. `PROST` consists of two independent workflows: **PROST Index (PI)** and **PROST Neural Network (PNN)**. 


![figure1](./docs/imgs/figure/figure1.png)

Using `PROST` you can do:
* Quantitative identification of spatial patterns of gene expression changes by the proposed **PROST Index (PI)**.

* Unsupervised identification of spatial tissue domains using a **PROST Neural Network (PNN)**. 
---

## Installation
### 1. Prepare `Python` environment
To install `PROST`, we recommend using the [Anaconda](https://anaconda.org/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) Python Distribution and creating an isolated environment, so that the `PROST` and dependencies don't conflict or interfere with other packages or applications. Please install the `conda` in advance. 


### Create environment 
First, please download or clone the `PROST-master.zip` file from github `Code` and unzip it. 

    git clone https://github.com/Tang-Lab-super/PROST.git

The entire installation process takes place in the `PROST-master` directory, so first go to that directory by:
   
    cd PROST-master

We recommend using a conda environment to configure PROST. To create and activate the environment `PROST_ENV`:  
**a.** If you use `Linux` or `Windows`, Run the following command in `bash` or `Anaconda Powershell Prompt`:

    conda create -n PROST_ENV python=3.7
    conda activate PROST_ENV

**b.** If you use `MacOS`, you may run the command instead:

    conda config --env --set subdir osx-64
    conda create -n PROST_ENV python=3.7 -c conda-forge
    conda activate PROST_ENV

### 2. Prepare `R` environment
The `PROST` uses the `mclust` package in the `R` language environment, and links it in a `Python` environment via `rpy2`. You can install the `R` language environment under `PROST_ENV` environment by:

    conda install -c conda-forge r-base
    conda install -c conda-forge r-mclust=5.4.10


### 3. Install dependency packages 
**a.** If you want to install `PROST` in `Linux` environment, you can install the dependency packages using `pip` by:
   
    pip install -r requirements.txt

**b.** If you want to install `PROST` in `Windows` environment, you can install the dependency packages using `pip` by:

    pip install -r requirements_win.txt
    pip install rpy2-2.9.5-cp37-cp37m-win_amd64.whl

**c.** If you want to install `PROST` in `MacOS` environment, you can install the dependency packages using `pip` by:
    
    pip install -r requirements_mac.txt
    RPY2_CFFI_MODE=BOTH pip3 install rpy2

### 4. Install `PROST`
Install the `PROST` package under `PROST_ENV`environment by:

    pip install setuptools==58.2.0
    python setup.py build
    python setup.py install

Here, the environment configuration is completed！If you are using `MacOS`, you may get a `segmentation fault` when running `PNN`, which may be a compatibility issue. You may try it running as a `script`. 

---

## How to use `PROST`
Before you use `PROST`, you have to make sure that the following two steps are taken place:  
* `PROST_ENV` environment is activated; 
* Add `environment variables` using the following `python code` before using `PROST`：

        import os
        ENVpath = "your path of PROST_ENV"  
        os.environ['R_HOME'] = f'{ENVpath}/lib/R'
        os.environ['R_USER'] = f'{ENVpath}/lib/python3.7/site-packages/rpy2'

**Note: using `conda info -e` can get `your path of PROST_ENV`.**

---

## Turorials
### Quick Start
* After `PROST` installation, we suggest downloading the complete tutorial examples from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7827565.svg)](https://doi.org/10.5281/zenodo.7827565) (The dataset is too large to upload to github, there only 1 case (151672) of DLPFC data).   
Similarly, you can download the dataset for each turorial individually via the `google drive` below. After you have downloaded the folder, unzip it and move the files to the datasets.

* Before running the tutorial, make sure your current path is in `/PROST-master/test`. **This facilitates you to perform PROST tests quickly.** In this folder, `data` in `datasets` and `result` will be stored under the corresponding data folder `results`. 

        cd ./test

For more flexibility, you may need to modify the `rootdir` to make sure the paths are correct.

### [1.Application on 10x Visium human dorsolateral prefrontal cortex (DLPFC) dataset.](./docs/tutorials/DLPFC.md "In this vignette, we analyzed tissue section from the human dorsolateral prefrontal cortex (DLPFC) 10x Visium ST dataset, which was manually annotated as the cortical layers and white matter (WM)") 
* **1.1** We performed `PROST` on the 10x Visium human dorsolateral prefrontal cortex (DLPFC) dataset from [(Pardo B. et al. 2022)](https://doi.org/10.1186/s12864-022-08601-w).
* **1.2** The original data and manual annotation are prepared and aviliable at (https://drive.google.com/drive/folders/1a_barvf0oUXDqyivWUUmv5bu9KxUlYoo) .
* **1.3** A Jupyter Notebook of the tutorial is accessible from : ([DLPFC.ipynb](./docs/vignettes/DLPFC.ipynb))

### [2. Application on Stereo-seq mouse olfactory bulb dataset.](./docs/tutorials/Stereo-seq.md "In this vignette, we analysis an ST dataset with cellular resolution (~14 μm in diameter per spot) generated by the Stereo-seq platform from mouse olfactory bulb tissue (add citation) to evaluate the performance of PROST on ST datasets with single-cell resolution.")
* **2.1** We performed `PROST` on the Stereo-seq mouse olfactory bulb dataset from [(Huazhu Fu. et al. 2021)](https://doi.org/10.1101/2021.06.15.448542) generated by [Stereo-seq](https://doi.org/10.1016/j.cell.2022.04.003).
* **2.2** We followed [(Kangning Dong. et al. 2022)](https://doi.org/10.1038/s41467-022-29439-6) to remove the spots outside the tissue section. The processed data can be downloaded from (https://drive.google.com/drive/folders/1lcGdArvtBcXkA0vE-v2GyQgeBrlZUpto?usp=share_link)
* **2.3** A Jupyter Notebook of the tutorial is accessible from : ([Stereo-seq.ipynb](./docs/vignettes/Stereo-seq.ipynb))

### [3. Application on SeqFISH mouse embryo dataset.](./docs/tutorials/SeqFISH_mouse_embryo.md "In this vignette, we applied PROST onto a SeqFISH-profiled dataset to evaluate its general applicability.")
* **3.1** We performed `PROST` on the processed mouse embryo ST data from [(Lohoff, T. et al. 2022)](https://doi.org/10.1038/s41587-021-01006-2) generated by [SeqFISH](https://spatial.caltech.edu/seqfish/).
* **3.2** The original data can be downloaded from (https://drive.google.com/drive/folders/1r5Cuo4YqZVtFPpB-VqDoHkMW6vhl9glf?usp=share_link)
* **3.3** A Jupyter Notebook of the tutorial is accessible from : ([SeqFISH_mouse_embryo.ipynb](./docs/vignettes/SeqFISH_mouse_embryo.ipynb))

---

## Improvements
We welcome any comments about `PROST`, and if you find bugs or have any ideas, feel free to leave a comment [FAQ](https://github.com/Tang-Lab-super/PROST/labels/FAQ).

---
