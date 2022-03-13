**XRFast Package for Dictionary learning of MA-XRF Images** <br>

***

This repository is a companion to the following paper:

**XRFast a New Software Package for Processing of MA-XRF Datasets using Machine Leaning** 

Marc Vermeulen1, Alicia McGeachy1, Bingjie Xu1, Henry Chopp2 Aggelos Kataggelos2, Rebecca Meyers3, Matthias Alfeld4, and Marc Walton5,a,*
1 Northwestern University / Art Institute of Chicago Center for Scientific Studies in the Arts (NU-ACCESS), 2145 Sheridan Road, Evanston, IL, United States
2 Department of ECE, Northwestern University, Evanston, IL, United States
3 National Museum of Mexican Art, 1852 W. 19th street, Chicago, IL, United States
4 Delft University of Technology, Department of Materials Science and Engineering (MSE), 2628 CN Delft, Netherlands.
5 M+, 38 Museum Drive, West Kowloon Cultural District, Hong Kong
a formerly: NU-ACCESS
* Corresponding author: marc.walton@mplus.org.hk

**Abstract:** *To be put here*

***

The aim of this repository is to give full access to the Jupyter notebook, associated macros, test dataset, and instructions allowing the reduction and visualization of visible spectral images of works of art as well as pigments identification and their spatial distribution using UMAP, 2D histogram peak finding algorithm and non-negative least square fitting to the original data. 

The intention is that the presented research can be fully replicated and implemented by other scientists in their host institutions. 

Efforts have been made to the best of the authors' abilities to make the data processing procedures accessible and reusable in support of the growing Open Science movement. 

1. Overview of contents

1.1. Jupyter notebook (Spectral Imaging Data Treatment - UMAP)

Most of the data treatment is undertaken in the Jupyter notebook. To successfully produce the final text figures, cells within the notebook should be run sequentially. Further pertinent instructions are provided throughout the notebook.

Required versions: Jupyter notebook >=5.5, Python >=3.4.

1.2. Install requirements

All python packages required to run the notebook can be install through the Python 3 or Anaconda navigator command prompt.

1.3. Test dataset

The Gauguin Poèmes Barbares spectral stack used in the associated article can be downloaded here 

1.4. RGB conversion scripts

LambdaStack_to_XYZ 
XYZ_to_RGB 
Refer to NU-ACCESS/Spectral-Microscope-Tools to download the required scripts

LambdaStack_to_XYZ 

converts a stack of monochromatic wavelength images into an XYZ tristimulus space. This is accomplished by multiplying each wavelength image by the CIE 1931 Standard Observer matching function, thus producing X, Y, and Z tristimulus value images. The input parameters are the spacing between wavelengths (dependant on the camera used for acquisition, default: 2) as well as the starting and ending wavelengths (default = 393 and 750 nm).
XYZ_to_RGB 

converts the XYZ stack into adobe RGB color space assuming a standard D65 Illuminant.
The scripts here were written in Python for use in ImageJ (Fiji) and are based on algorithms detailed in Oakley et. al. "Improved spectral imaging microscopy for cultural heritage through oblique illumination", Heritage Science, 8:1 

To run the scripts, download the latest version of Fiji and place the scripts into the plugins folder. Launch Fiji or "refresh menus" if already running. 
Both scripts will be available in the Plugins drop-down menu. 
The scripts should be ran in the following order: 1) LambdaStack_to_XYZ and 2) XYZ_to_RGB. 
The final RGB composite image should be transformed into a RGB color image (Image>Type>RGB color) prior to be saved as a TIF image (File>Save as>Tiff)

2. Notebook

Markdown cells provide the instructions and expected output of the subsequent cells of code. When numbers are provided in the markdown cells, they correspond to the position of the cell in the subsequent block of cells of code: "(1) xxxxx" corresponds to the instructions for the first following code cell whereas "(2) xxxx" corresponds to the instruction of the second cell of code.

3. Endmember spectra

Endmember spectra are exported as a .csv file using space as delimiter. They can be imported in Excel for further data processing.

4. Creating the endmember distribution maps from the .txt text image files saved

A stack of endmember distribution maps can be created in ImageJ (Fiji). To create them, follow the following steps:

Open the save txt file as a text image (File>Import>Text Image...)
Select the desired txt file
Transform the txt montage to a stack (Image>Stacks>Tools>Montage to Stack)
In the Stack Maker pop-up window, input: 
Columns: the Y value of the original data set dimensions (it will be given when running the Import the .tif spectral data cell of the Notebook) 
Rows: 1 
Click "OK" 
Reslice the stack from the top to create a stack of the distribution maps (Image>Stacks>Reslice).
Select outspacing (pixels): 1 
Start at: Top 
Tick "Avoid interpolation (use 1 pixel spacing)
Click OK 
The stack of images can be broken into individual images (Image>Stack>Stack to Images)
[![Build Status](https://travis-ci.com/NU-ACCESS/XRFast2.jl.svg?branch=master)](https://travis-ci.com/NU-ACCESS/XRFast2.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/NU-ACCESS/XRFast2.jl?svg=true)](https://ci.appveyor.com/project/NU-ACCESS/XRFast2-jl)
[![Codecov](https://codecov.io/gh/NU-ACCESS/XRFast2.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NU-ACCESS/XRFast2.jl)
[![Coveralls](https://coveralls.io/repos/github/NU-ACCESS/XRFast2.jl/badge.svg?branch=master)](https://coveralls.io/github/NU-ACCESS/XRFast2.jl?branch=master)
