# Notebook dedicated to making plots with GO terms / GOSLIM terms, given a list of TAIR IDs

## Create necessary environment

If micromamba is not installed (it takes 1 minute to install)                                                                                                                                                                                                                                                                                                                                                                                                                                    
````bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
````
full documentation at:
https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

````bash
micromamba create -n <env name> -f conda-linux-64.lock
micromamba activate <env name>
# install additional packages with pip
pip install -r requirements.txt
````

## Running

Replace in the notebook the path to the file containing your TAIR IDs.
````bash
cat tair_ids.csv
AT3G28510.1
AT1G34790
AT3G54150.1
````
.1, .2, etc. (splicing variant numbers) will be removed during the process

Then run the whole notebook!