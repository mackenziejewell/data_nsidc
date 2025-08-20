# data_nsidc
Scripts for opening, read data files from the National Snow and Ice Data Center (NSIDC)


### To build matching environment:

(1) Navigate to folder in terminal and run:
conda env create --file=environment.yml

(2) Optional: if running in Jupyter Notebook environment
After creating the environment, you will next need to add it to the list of available kernels. 
This does not happen automatically. In the command line, run:
python -m ipykernel install --user --name <NAME>
* Make sure to replace <NAME> with the name of the environment

If conda command not recognized by default with terminal launch:
source /opt/anaconda3/bin/activate

(3) Optional, to update environment file after modifying packages:
conda env export > environment.yml


# icedrift.py

Polar Pathfinder Daily 25 km EASE-Grid Sea Ice Motion Vectors (NSIDC-0116)
This opens nc files from the Polar Pathfinder Daily 25 km EASE-Grid Sea Ice Motion Vectors, Version 4

Data set id:
NSIDC-0116
DOI: 10.5067/INAWUWO7QH7B

https://nsidc.org/data/nsidc-0116/versions/4#anchor-data-access-tools 


# icesat2.py

Various data products for ICESat-2 data hosted on NSIDC, including:

ICESat-2 L4 Monthly Gridded Sea Ice Thickness, Version 3
Data set id: IS2SITMOGR4
DOI: 10.5067/ZCSU8Y5U1BQW
