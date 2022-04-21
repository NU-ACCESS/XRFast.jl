# The XRFast Package for Dictionary Learning of MA-XRF Images <br>

***

This repository is a companion to the following paper:

**XRFast a New Software Package for Processing of MA-XRF Datasets using Machine Leaning** 
> *Journal where the paper will be pubished*,  <br>

Marc Vermeulen (1), 
Alicia McGeachy (1), 
Bingjie Xu (1), 
Henry Chopp (2),
Aggelos Kataggelos (2), 
Rebecca Meyers (3), 
Matthias Alfeld (4), 
Marc Walton (5)<br>


> 1. Northwestern University / Art Institute of Chicago Center for Scientific Studies in the Arts (NU-ACCESS), 2145 Sheridan Road, Evanston, IL, United States
> 2. Department of ECE, Northwestern University, Evanston, IL, United States
> 3. National Museum of Mexican Art, 1852 W. 19th street, Chicago, IL, United States
> 4. Delft University of Technology, Department of Materials Science and Engineering (MSE), 2628 CN Delft, Netherlands.
> 5. M+ Museum, 38 Museum Drive, West Kowloon Cultural District, Hong Kong*

*corresponding author: marc.walton@mplus.org.hk 

**Abstract:** *XRF is a common technique in the field of heritage science. However, data processing and data interpretation remains a challenge because often requires a-priori knowledge of the pigments present in the analyzed objects. For this reason, we developed an open-source unsupervised dictionary learning algorithm reducing the complexity of large datasets containing 10’s of thousands of spectra and identifying patterns. This approach reduces the number of variables and creates correlated elemental maps, characteristic for pigments containing various elements or for pigments mixtures. This alternative approach creates an overcomplete dictionary which is learned which is based on the input data itself, therefore reducing the a-priori user input. In this study, we apply this method to a MA-XRF data set obtained on an 18th century Mexican painting, and positively identified smalt (pigment characterized by the co-occurrence of cobalt, arsenic, bismuth, nickel, and potassium), mixtures of vermilion and lead white, identified two complex conservation materials/interventions, and identified correlated elements that were not identified using the traditional elemental maps approach without image processing. This approach proved very useful as it yields the same conclusions as the traditional elemental maps approach without requiring further image processing or user manipulation to understand elemental correlation. This is therefore a gain of time and a useful resource to understand better the pigments mixtures used in historical paintings.*

***

The aim of this repository is to give full access to a Julia package called XRFast designed for dictionary learning and sparse representation of XRF signals using KSVD. Also included is a Jupyter notebook showing how to use this package, a test dataset, and instructions allowing for the identification of spectral endmembers and their spatial distribution across a 2D work of art.

The intention is that the presented research can be fully replicated and implemented by other scientists in their host institutions. 

Efforts have been made to the best of the authors' abilities to make the data processing procedures accessible and reusable in support of the growing Open Science movement. <br>

## 1. Overview of contents
#### 1.1. Jupyter notebook (XRF Dictionary learning)
Most of the data treatment is undertaken in the Jupyter notebook. To successfully produce the final text figures, cells within the notebook should be run sequentially. Further pertinent instructions are provided throughout the notebook.

Required versions: [Jupyter notebook](https://jupyter.org/)>=5.5, [Julia](https://julialang.org/downloads/)>=1.5 

#### 1.2. Install requirements
All packages required to run the notebook can be installed within the Jupyter notebook or using the Julia command prompt. <br>
To install a package, use the following lines: <br>
> using Pkg <br>
> Pkg.add("Package Name") <br>

To install the <b>XRFast package</b>, use the following import line: 
>Pkg.add(url="https://github.com/NU-ACCESS/XRFast.jl") <br>

#### 1.3. Test dataset
The XRF data cube from a selected area of an untitled 18th century Mexican painting belonging to the permanent collection of the National Museum of Mexican Art (Chicago, IL, USA) used in the article is [provided](https://.../) as a test dataset. <br>

## 2. Notebook
Markdown cells provide the instructions and expected output of the subsequent cells of code. 

## 3. Endmember distribution maps and spectra
Endmember distribution maps and spectra are exported as a .csv files using comma as delimiter. <br>

<li><b>Endmember spectra</b> can be imported in Excel or other graphing applications for further data processing. <br>
<li>A stack of <b>endmember distribution maps</b> can be created in ImageJ (Fiji) by following the steps below (also found in the Jupyter Notebook): </ol><br>

<ol>
<li><b>Import data</b>: File>>Import>>text Image...<br>
<i>Select the desired .csv file</i></li>
<li><b>Create a stack</b>: Image>>Stacks>>Tools>>Montage to Stack...<br>
<ins>Columns</ins>: last number of the filename, before the .csv <br>
<ins>Rows</ins>: 1 <br>
<ins>Border width</ins>: 0</li>
<li><b>Transform compressed stack into stack of images</b>: Image>>Stacks>>Reslice[/]... (or Shortcut Stk>>Reslice[/]...) on the <i>Stack</i> window created</li>
<ins>Output spacing (pixels)</ins>: 1.000 <br>
<ins>Start at</ins>: Top <br>
<ins>Avoid interoplation</ins> should be checked <br>
<i>Flip vertically</i> and <i>Rotate 90 degrees</i> can be checked or unchecked depending on the orientation of the original data. 
<li><b>Save as a set of individual images</b>: File>>Save As>>Image Sequence... <br>
<ins>Dir.</ins>: Select the folder where images will be saved <br>
<ins>Format</ins>: select the desired format (TIFF is the default option) <br>
<ins>Name</ins>: Indicate the root name of the images + EM (for endmember, recommended)<br>
<ins>Start at</ins>: Indicate 1 (this correspond to the first endmember, which is the first slice)<br>
<ins>Digits (1-8)</ins>: we recommend 2. This will save each image as EM01, EM02, [...] EM20).

***

[![Build Status](https://travis-ci.com/NU-ACCESS/XRFast2.jl.svg?branch=master)](https://travis-ci.com/NU-ACCESS/XRFast2.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/NU-ACCESS/XRFast2.jl?svg=true)](https://ci.appveyor.com/project/NU-ACCESS/XRFast2-jl)
[![Codecov](https://codecov.io/gh/NU-ACCESS/XRFast2.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NU-ACCESS/XRFast2.jl)
[![Coveralls](https://coveralls.io/repos/github/NU-ACCESS/XRFast2.jl/badge.svg?branch=master)](https://coveralls.io/github/NU-ACCESS/XRFast2.jl?branch=master)
