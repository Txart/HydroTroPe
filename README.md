# Installation
## venv + pip
This is the preferred method now.
(Instructions for Linux machines)

Create a virtual environment with venv:
`python3 -m venv [virtual_environment_name]`

Activate the virtual environment
`source [virtual_environment_name]/bin/activate`

Install the requirements
`python3 -m pip install -r ./requirements.txt`


## Conda
It was developed in Windows 10 and tested in Ubuntu with the conda environment that results from the following installation (some of the packages contained in the yaml might not be strictly necessary):
```
conda env create -f environment.yml
```

After the installation, activate the environment with 
```
conda activate hydrotrope
```




-------------
(The name for this repo was created using openAI's ChatGPT. It contains a joke, and I think that's pretty noteworthy.)
