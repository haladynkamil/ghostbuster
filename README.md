# ![](https://i.imgur.com/oDOcrQA.png)
![Python](https://img.shields.io/badge/python-v3.6+-blue.svg)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
## Introduction
Ghostbuster is a python software that allows you to plate solve images taken by savart plates. It takes your 'doubled' image, finds parity stars and solves astrometry for one set of stars.
## Table of contents
* [Introduction](#introduction)
* [Technologies](#technologies)
* [Setup](#setup)
## Technologies

- Python 3.6.8
- [ds9](http://ds9.si.edu/site/Home.html) (or equivalent)
- [astropy](https://www.astropy.org/)
- [numpy](https://numpy.org/)
- [photutils](https://photutils.readthedocs.io/en/stable/)
- [matplotlib](https://matplotlib.org/)
- [docker](https://www.docker.com/)
- [scipy](https://www.scipy.org/)

You can install all packages using requirements.txt:

```
pip install --user --requirement requirements.txt
```

## Setup
First of all you need python 3.6.8. To install python on your PC use [this guide](https://realpython.com/installing-python/)
After that you will have to download all libraries pointed in "Technologies" chapter. Each of them has its own documentary that has an install guide. All of them can be installed using PIP.

After that, you need to download indexes for astrometry. Indexes are slices of the sky that are compared to our image by astrometry software. You can find them [here](http://broiler.astrometry.net/~dstn/4200/) or here if you have [wide-field images](http://broiler.astrometry.net/~dstn/4100/)

Third important thing that is required to use this software is Docker. Docker allows you to run a container with wrap-up software, and here we are going to use astronova container. [Here is the guide](https://www.docker.com/get-started) how to get started with docker.

After we set up docker we have to run nova container and mount the volumens.
To do this use this code: 

```
sudo docker run -d --name nova -v <PATH_TO_INDEX>:/usr/local/astrometry/data -v <PATH_TO_DATA>:/data_market -p 8000:8000 michalzg/nova
```

**PATH_TO_INDEX** --> Provide a system path to index files downloaded earlier,

**PATH_TO_DATA** --> This is a folder where Ghostbuster will output coordination tables, and from where astronova will take them. I highly reccomend using output_tables/ folder just like in this git repository.


## Example
Start with setting up your config.ini file. 
To get x_distance and y_distance we reccomend to use ds9 and simply to create the rectangular box between star and her ghost.

![ds9 setup](https://i.imgur.com/XL2l18v.png)

Try to be as precise as possible in finding those pixel distances. It's one of the most important thing.

Here are some default values:
```
[DEFAULT]

fwhm = 3  
threshold = 200  <- Lower this value if you want to find more stars on the image.
sigma = 3  
angle = 30      
savart_distance = 39.2  
min_star_limit = 10  <- minimal limit of stars found on image to process it  
approx_distance = 10   
norm_percent = 99.9  

[SAVART_PARAMETERS_P1]  
x_distance = 27  <- Distance between real star and her ghost in X dimension (in pixels)  
y_distance = 32  <- Same for Y dimension  
  
[SAVART_PARAMETERS_P3]           <- Parameters for 2nd savart plate  
x_distance = 1  
y_distance = 38  

[EXECUTABLE_PARAMETERS]  
plot_no_mask = False    <- if you want to plot your picture without mask  
plot_mask = False   <- if you want to plot your picture with masked ghosts  
mask_square_size = 7  <- width of the masking square  
search_distance = 3  <- distance in pixels, this is the error measure for ghost position (we can diverge from our position by 3px                 in this example)  
save_wcs = False  
  
[PLOT_PARAMETERS]    
figsize_rcparams_width = 10  
figsize_rcparams_height = 20  
alpha = 0.45  
draw_aperture = True  
aperture_size = 4  
 ```
 
After creating your config.ini we can process our picture.  
To do so, simply run the script with path and type of savart plate (p1 or p3)  
Example:  

```python3 savart.py ../images/observations/savartp1/image0001.fits p1```

Results will depend on your config file. 
Assuming that all plots will be created (you can skip unwanted plots using config file), ghostbuster will start by plotting your picture.

After that, masks will be added to every star, and then masks will be moved by direction provided by config file.

![Moved masks](https://i.imgur.com/qgvB33o.png)

All unmasked stars are going to be found again. Those are our primary stars. Software wil provide those stars as X, Y points for astrometry solution. After the solution, WCS file will be loaded and all ghosts and real stars will have their RA and DEC. Results will be outputed in a table format. 

![Example of table content](https://i.imgur.com/55WR8B9.png)

**If astrometry fails, try changing the 'treshold' value** 


