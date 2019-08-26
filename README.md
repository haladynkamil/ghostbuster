# ghostbuster
## Introduction
Ghostbuster is a python software that allows you to plate solve images taken by savart plates. It takes your 'doubled' image, finds parity stars and solves astrometry for one set of stars.
## Table of contents
* [Introduction](#introduction)
* [Technologies](#technologies)
* [Setup](#setup)
## Technologies

- Python 3.6.8
- ds9 (or equivalent)
- astropy
- numpy
- photutils
- matplotlib

## Setup

## Example
Start with setting up your config.ini file. 
Here are some default values:

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
  
After creating your config.ini we can process our picture.  
To do so, simply run the script with path and type of savart plate (p1 or p3)  
####Example:  
python3 savart.py ../images/observations/savartp1/image0001.fits p1


