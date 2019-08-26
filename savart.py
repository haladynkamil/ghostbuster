import random
import astropy
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from os import listdir
from os.path import isfile, join
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import simple_norm
from astropy.table import Table, setdiff
import math
import configparser
import os
import sys
import logging
def open_image(my_path):
    """ Returns fits data from path provided by the user.
    """

    image_data = fits.getdata(my_path)
    return image_data

def create_mask(image_data):
    mask = np.zeros_like(image_data, dtype=bool)
    mask[...] = False
    return mask


def star_data_gatherer(image_data, mask, **kwargs):
    """ finds stars based on config parameters and returns list sorted by flux
    """

    fwhm = kwargs['fwhm']
    threshold = kwargs['threshold']
    sigma = kwargs['sigma']

    mean, median, std = sigma_clipped_stats(image_data, sigma=sigma) #calculation of mean, median and standard deviation from the data
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(image_data - median, mask=mask)
    norm = simple_norm(image_data, 'sqrt', percent=99.9)
    sources.sort('flux')

    return sources, norm, mask


def plot_photo(image_data, mask, norm, **kwargs):
    """ Plots a photo from fits image data.
    """
    plt.rcParams['figure.figsize'] = (20,10)
    plt.imshow(image_data, cmap='Greys',  norm=norm)
    plt.imshow(mask, interpolation='none',  alpha=0.45)
    plt.show()


def starmasking(coordinates_chart, mask):
    for row in coordinates_chart:
        print(row)
        mask[int(row['ghost_ycentroid'])-7:int(row['ghost_ycentroid'])+7, int(row['ghost_xcentroid'])-7:int(row['ghost_xcentroid'])+7] = False
    return mask

    
def get_ghost_distance(**kwargs):

    savart_radian_angle = np.deg2rad(kwargs['angle'])  #changes angle value of star line to radians
    x_distance = kwargs['savart_distance'] * np.cos(savart_radian_angle) #calculates the distance in X axis and Y axis from possible star ghosts
    y_distance = kwargs['savart_distance'] * np.sin(savart_radian_angle)

    return x_distance, y_distance

def ghostbuster_2(sources, mask, **kwargs):
    right_x, right_y = kwargs['right_x'],kwargs['right_y']
    for row in sources:
        mask[int(row['ycentroid']+right_y)-8:int(row['ycentroid']+right_y)+8, int(row['xcentroid']+right_x)-8:int(row['xcentroid']+right_x)+8] = True
    return mask

    
def sources_to_list(sources, ghost_xy, file_name):
    file_name = 'coordinates_output_table/' + os.path.basename(file_name)
    
    
    sources.write(file_name, format='ascii', include_names=['xcentroid','ycentroid'], overwrite=True)
    ghost_xy.write(file_name + '_ghosts', format='ascii', overwrite=True)


#-------------------DIFFERENT-APPROACH------------------------->


def ghostbuster(all_sources, clean_sources, **kwargs):
    """ Ghostbuster finds stars duplicated by savart plate
    """
    
    chosen_stars = clean_sources.copy() #creates empty astropy.table object, where we can store our chosen stars 

    coordinates_chart = Table(names=('real_xcentroid', 'real_ycentroid', 'ghost_xcentroid', 'ghost_ycentroid'))
    ghost_sources = setdiff(all_sources, clean_sources, keys=['xcentroid', 'ycentroid'])
    
    for star in chosen_stars:
        """ using ghostfinder to find a star that has te closest values of x and y distances
        """
        
        real, ghost = ghostfinder(star, ghost_sources, **kwargs)
        if ghost == 'no_match':
            pass
        else:
            coordinates_chart.add_row([real['xcentroid'], real['ycentroid'], ghost['xcentroid'], ghost['ycentroid']])
    

    return(coordinates_chart)


def ghostfinder(star, ghost_sources, **kwargs):
    """ Used by ghostbuster to find the star and define which one is a ghost
    """

    distance = 3  #variable for cleaner coding
    right_x, right_y, = star['xcentroid'] + kwargs['right_x'], star['ycentroid'] + kwargs['right_y']
    i = 0
    for possible_ghost in ghost_sources:
        i += 1
        if abs(possible_ghost['xcentroid'] - right_x) < distance and abs(possible_ghost['ycentroid'] - right_y) < distance:
            return star, possible_ghost
        
    return star, 'no_match'
    



def get_hideout_coordinates(star, x_distance, y_distance, **kwargs):
    """ Finds new coordinates in appropiate quadrants around the star
    """
    
    angle = kwargs['angle']
    
    if angle > 0 and angle < 90:

        right_x = star['xcentroid'] + x_distance
        left_x = star['xcentroid'] - x_distance
        right_y = star['ycentroid'] + y_distance
        left_y = star['ycentroid'] - y_distance

    elif angle > 90:
        
        right_x = star['xcentroid'] - x_distance
        left_x = star['xcentroid'] - x_distance
        right_y = star['ycentroid'] + y_distance
        left_y = star['ycentroid'] + y_distance

    return right_x, right_y, left_x, left_y
    

def get_apertures(image_data, sources):
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.)
    apertures.plot()
    return apertures

def read_config():
    path = 'config.ini'
    config = configparser.ConfigParser()
    config.readfp(open('config.ini'))
    config_dict = {
            'fwhm': int(config['DEFAULT']['fwhm']),
            'threshold': int(config['DEFAULT']['threshold']),
            'sigma': int(config['DEFAULT']['sigma']),
            'angle': config['DEFAULT']['angle'],
            'savart_distance': config['DEFAULT']['savart_distance'],
            'min_star_limit': config['DEFAULT']['min_star_limit'],
            'right_x': float(config['SAVART_PARAMETERS']['x_distance']),
            'right_y': float(config['SAVART_PARAMETERS']['y_distance']),
            'plot_no_mask': config['EXECUTABLE_PARAMETERS']['plot_no_mask'],
            'plot_mask': config['EXECUTABLE_PARAMETERS']['plot_mask'],
            }

    return config_dict

def exec():

    logging.basicConfig(filename='logger.log', level=logging.DEBUG)
    file_name = sys.argv[1]
    if os.path.isfile(file_name):
        logging.info('Opening file: ', os.path.basename(file_name))
        print(os.path.basename(file_name))
    else:
        logging.ERROR('File does not exist')
        print('File does not exist')
        sys.exit()

    image_data = open_image(file_name)
    logging.info('Opened image_file')
    config = read_config()
    all_sources, norm, mask = star_data_gatherer(image_data, create_mask(image_data), **config)
    logging.info('Created all sources list')
    if config['plot_no_mask'] == 'True':
        plot_photo(image_data, mask, norm, **config)
        logging.info('Ploted image without masks')
    mask = ghostbuster_2(all_sources, mask, **config)
    logging.info('Masked ghosts')
    clean_sources, norm, mask = star_data_gatherer(image_data, mask, **config)
    logging.info('Found real stars')
    get_apertures(image_data, clean_sources)
    if config['plot_mask'] == 'True':
        plot_photo(image_data, mask, norm, **config)
        logging.info('Plotted image with masked ghosts')
    ghost_xy = ghostbuster(all_sources, clean_sources, **config)
    logging.info('Found ghosts')
    sources_to_list(clean_sources, ghost_xy, file_name)
    logging.info('Saved normal and ghosts list')
if __name__ == '__main__':
    exec()








