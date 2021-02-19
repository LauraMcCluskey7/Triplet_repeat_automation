# Triplet_repeats



## Download the directory:

```
git clone https://github.com/LauraMcCluskey7/Triplet_repeats.git 

```


##Create and activate the conda environment

```
conda env create -f triplets_environment.yml

conda activate triplets_environment

```



## Create the executable:


The executable file must be created on a windows platform


```
pyinstaller triplet_repeat_automation.py

```


## Run the executable:


* Click on the executable 
* Input the gene name and worksheet number into the console





## Run the tests:

```

python -m unittest test_triplet_repeats.py

```

