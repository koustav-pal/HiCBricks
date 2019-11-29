Travis CI: [![Build Status](https://travis-ci.com/koustav-pal/HiCBricks.svg?branch=master)](https://travis-ci.com/koustav-pal/HiCBricks)   

[<img src="https://www.bioconductor.org/images/logo/jpg/bioconductor_logo_rgb.jpg" width="300" align="right"/>](https://bioconductor.org/)

<br/><br/><br/><br/>

<a href="https://bioconductor.org/checkResults/devel/bioc-LATEST/HiCBricks/"><img src="https://bioconductor.org/shields/build/devel/bioc/HiCBricks.svg" align="right"/></a> <a href="https://bioconductor.org/packages/devel/bioc/html/HiCBricks.html#since"><img src="https://bioconductor.org/shields/years-in-bioc/HiCBricks.svg" align="right"/></a> <a href="https://bioconductor.org/packages/devel/bioc/html/HiCBricks.html#archives"><img src="https://bioconductor.org/shields/availability/3.9/HiCBricks.svg" align="right"/></a>

<br/><br/>

![My Image](https://user-images.githubusercontent.com/20904402/55158335-ef97f600-515e-11e9-8b84-f8557428da70.png)

# HiCBricks

HiCBricks is a **R/Bioconductor** package for handling high-resolution Hi-C datasets through HDF (Hierarchical Data Format) files. Read more about HDF [here](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) 

- HiCBricks greatly simplifies user handling of Hi-C contact matrices. 
- Forces users to adhere to a set of Hi-C analysis good-practices.
- HiCBricks simplifies how users interaction with HDF files containing Hi-C contact matrices.

## Features

- Import Hi-C data in multiple data formats. Currently, NxN dimensional matrices and _mcool_ files are supported, with more to come.
- Fetch different subset of the Hi-C data by their features with easy to use functions. Feature examples: by distance, matrix squares, rows or columns.
- Keep user-defined annotations associated to the HDF files.
- Use HiCBricks accessors to build more complex analysis such as TAD calling and visualizations.

## Installation

To install the most stable development version from Bioconductor, run this from a R console. **Note:** `R version >= 3.5` is required. This command will first installs `BiocManager` from CRAN. `BiocManager` is a convenient utility to install `Bioconductor` packages. Then, we install `HiCBricks` through `BiocManager`.
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("HiCBricks", version = "3.10")
```


To install the stable release version from Bioconductor, run this from a R console. 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("HiCBricks", version = "3.10")
```

To install the most cutting-edge stable version of HiCBricks, do this from a R console to download it directly from GitHub. 
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_git("https://github.com/koustav-pal/HiCBricks")
```

## Getting Started

To start working with **HiCBricks**, please checkout the vignette (tutorial) [here](http://bioconductor.org/packages/devel/bioc/vignettes/HiCBricks/inst/doc/IntroductionToHiCBricks.html), at the main Bioconductor website. It contains an in-depth walkthrough of almost all functions in **HiCBricks** and will guide users through the process of 

- Loading data from text 2D files.
- Loading data from mcool files.
- Making TAD calls and spohisticated heatmaps with example functions built using HiCBricks accessor functions.

## Development Notes

- **HiCBricks** API is now stable. While we may move to sparse or feather representations later, this API will not change.
- With Bioconductor release 3.10, a formal S4 class has been implemented for a better user experience.

## Future Roadmap

There are many new developments which are planned for future releases of HiCBricks. Broadly speaking, 

- In v1, I will try to implement `read and export functions` for as many new Hi-C data formats as possible. On top priority is the `sparse matrix`, HiCExplorer and diffHiC, in that exact order. Support for `.hic` is planned and requires some time.


## Contributing 

If you would like to help out, let me know via email.

