# Hairpin-Vortex-Extraction

> This repository is a work in progress...

We have uploaded the code for the interactive visualization system so far. We will upload the code for the remaining modules soon.

### How to run the interactive visualization code
In the root folder, do the following to install all the requirements for Python.
```
pip install -r requirements.txt
```
then in the vis_sys folder, open the command prompt and run the file *pyecharts_QtWebEngine.py*.
- Use the browse button to open the file from the folder *full_data*. 
- Doing so will automatically read the tree and profiles information from the *.json* file and the clustering infromation from the *_clusters.txt* associated with each dataset. 
- Currently *.vtk*, *.json* and *_clusters.txt* are only provided for Benard and Cylinder flow datasets. These files are the output of our vortex extraction and separation strategies explained in the paper.
