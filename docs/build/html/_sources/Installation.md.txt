# Installation
## 1. Prepare `Python` environment
To install `PROST`, we recommend using the [Anaconda](https://anaconda.org/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) Python Distribution and creating an isolated environment, so that the `PROST` and dependencies don't conflict or interfere with other packages or applications. Please install the `conda` in advance. 


### Create environment 
First, please download or clone the `PROST-master.zip` file from github `Code` and unzip it. 

    git clone https://github.com/Tang-Lab-super/PROST.git

The entire installation process takes place in the `PROST-master` directory, so first go to that directory by:
   
    cd PROST-master

We recommend using a conda environment to configure PROST. To create and activate the environment `PROST_ENV`, run the following command in `bash` or `Anaconda Powershell Prompt`:  

    conda create -n PROST_ENV python=3.7
    conda activate PROST_ENV

## 2. Prepare `R` environment
The `PROST` uses the `mclust` package in the `R` language environment, and links it in a `Python` environment via `rpy2`. You can install the `R` language environment under `PROST_ENV` environment by:

    conda install -c conda-forge r-base
    conda install -c conda-forge r-mclust=5.4.10


## 3. Install dependency packages 
**a.** If you want to install `PROST` in `Linux` environment, you can install the dependency packages using `pip` by:
   
    pip install -r requirements.txt

**b.** If you want to install `PROST` in `Windows` environment, you can install the dependency packages using `pip` by:

    pip install -r requirements_win.txt
    pip install rpy2-2.9.5-cp37-cp37m-win_amd64.whl

## 4. Install `PROST`
Install the `PROST` package under `PROST_ENV`environment by:

    pip install setuptools==58.2.0
    python setup.py build
    python setup.py install

Here, the environment configuration is completed！

---