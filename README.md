# Hairpin-Vortex-Extraction

> This repository is a work in progress...

### How to run the code
- Download the code or clone the repository. In the first step, go the the folder "splitting" and follow the readme instructions.
- Then go to the clustering folder and perform clustering. Read the comments in the code for additional information.
- Move the clustering results, the tree json file and the vtk file from step 1 to the same folder and follow the instructions below to visualize the results in the system.

**_Required_ Python >= 3.9**

### How to run the interactive visualization code
In the root folder, do the following to install all the requirements for Python.
```
pip install -r requirements.txt
```
then in the vis_sys folder, run the file *pyecharts_QtWebEngine.py*.
- Use the browse button to open the file from the folder *full_data*. 
- Doing so will automatically read the tree and profiles information from the *.json* file and the clustering infromation from the *_clusters.txt* associated with each dataset. 
- Currently *.vtk*, *.json* and *_clusters.txt* are only provided for **Benard** and **Cylinder** flow datasets. These files are the output of our vortex extraction and separation strategies explained in the paper.
