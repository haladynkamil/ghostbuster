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
from datetime import date
from astropy.coordinates import SkyCoord
import subprocess as sub
from astropy import wcs

def get_hdr(hdr):
    ra = hdr['RA']
    dec = hdr['DEC']
    coo = SkyCoord(ra=ra, dec=dec, frame='icrs', unit=('hour','deg'))
    return coo

def run_astrometry(image_path, coo, **kwargs):
    logging.info('ASTROMETRY: start')
    image_base_name = os.path.basename(image_path)
    solve_field_command = [
        'docker', 'exec', 'nova', 'solve-field',
        '--ra', '%s' % coo.ra.deg,
        '--dec', '%s' % coo.dec.deg,
        '--radius', '1',
        '--depth', '700',
        '--cpulimit', '3000',
        '--scale-units', 'arcsecperpix',
        '--scale-low', '0.8',
        '--scale-high', '1.3',
        '--x-column', 'x',
        '--y-column', 'y',
        '--width', '1024',
        '--height', '1024',
        '--overwrite', '--no-verify', '--no-plots',
        '/data_market/' + image_base_name + 'new.fits']
    text_fits_command = [
        'docker', 'exec', 'nova', 'text2fits',
        '-f', 'ff',
        '-H', 'x y',
        '/data_market/' + image_base_name,
        '/data_market/' + image_base_name + 'new.fits']
    print(image_base_name)

    print(text_fits_command) 
    p = sub.Popen(text_fits_command, stdout=sub.PIPE,
            stderr=sub.PIPE)
    output, error = p.communicate()
    print(output)
    print(error)
    
    sub.Popen(solve_field_command, stdout=sub.PIPE,
            stderr=sub.PIPE).communicate()
    wcs_file_path = './coordinates_output_table/' + image_base_name + 'new.wcs'
    if os.path.exists(wcs_file_path):
        logging.info('ASTROMETRY: success')
        return wcs_file_path
    
    logging.info('ASTROMETRY: FAILED')
    return None, False

def get_ra_dec_table(ghost_table, wcs_path, file_name):
    print(wcs_path)
    hdr = fits.open(wcs_path)[0].header
    w = wcs.WCS(hdr)
    x = ghost_table['real_xcentroid'].data
    y = ghost_table['real_ycentroid'].data
    xyli = np.c_[x,y]
    world = w.wcs_pix2world(xyli, 1)
    ghost_table['RA'] = world[:,0]
    ghost_table['DEC'] = world[:,1]
    print(ghost_table)
    file_name = 'coordinates_output_table/' + os.path.basename(file_name).split('.')[0]
    ghost_table.write(file_name + '_ghosts', format='ascii', overwrite=True)

def open_image(my_path):
    """ Returns fits data from path provided by the user.
    """
    logging.info('Opening image data')
    image_data = fits.getdata(my_path)
    hdul = fits.open(my_path)
    hdr = hdul[0].header
    logging.info('Finished')
    return image_data, hdr


def create_mask(image_data):
    mask = np.zeros_like(image_data, dtype=bool)
    mask[...] = False
    return mask


def star_data_gatherer(image_data, mask, **kwargs):
    """ finds stars based on config parameters and returns list sorted by flux
    """
    logging.info('Finding stars on image...')
    fwhm = kwargs['fwhm']
    threshold = kwargs['threshold']
    sigma = kwargs['sigma']

    mean, median, std = sigma_clipped_stats(image_data, sigma=sigma) #calculation of mean, median and standard deviation from the data
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(image_data - median, mask=mask)
    norm = simple_norm(image_data, 'sqrt', percent=kwargs['norm_percent'])
    sources.sort('flux')
    logging.info('Found %s stars', str(len(sources)))
    logging.info('Finished')
    return sources, norm, mask


def plot_photo(image_data, mask, norm, **kwargs):
    """ Plots a photo from fits image data.
    """
    logging.info('Ploting...')
    plt.rcParams['figure.figsize'] = (20,10)
    plt.imshow(image_data, cmap='Greys',  norm=norm)
    plt.imshow(mask, interpolation='none',  alpha=kwargs['alpha'])
    plt.show()

    
def get_ghost_distance(**kwargs):
    savart_radian_angle = np.deg2rad(kwargs['angle'])  #changes angle value of star line to radians
    x_distance = kwargs['savart_distance'] * np.cos(savart_radian_angle) #calculates the distance in X axis
    y_distance = kwargs['savart_distance'] * np.sin(savart_radian_angle) #and Y axis from possible star ghosts
    return x_distance, y_distance


def ghostbuster_2(sources, mask, **kwargs):
    logging.info('Masking ghost with ghostbuster_2')
    right_x, right_y = kwargs['right_x'],kwargs['right_y']
    logging.info('X, Y shifts = %s %s', str(right_x), str(right_y))
    square_size = kwargs['square_size']
    for row in sources:
        mask[int(row['ycentroid']+right_y)-10:int(row['ycentroid']+right_y)+10,
                int(row['xcentroid']+right_x)-10:int(row['xcentroid']+right_x)+10] = True
    
    logging.info('Ghosts masked')
    return mask

    
def sources_to_list(sources, ghost_xy, file_name):
    logging.info('Saving normal and ghosts list')
    file_name = 'coordinates_output_table/' + os.path.basename(file_name).split('.')[0]
    sources.write(file_name, format='ascii', include_names = ['xcentroid','ycentroid'],overwrite=True)

    with open(file_name, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(file_name, 'w') as fout:
        fout.writelines(data[1:])
    ghost_xy.write(file_name + '_ghosts', format='ascii', overwrite=True)
    logging.info('Saved normal and ghosts list')
    return file_name


def ghostbuster(all_sources, clean_sources, **kwargs):
    """ Ghostbuster finds stars duplicated by savart plate
    """
    logging.info('Looking for ghosts started..')    
    chosen_stars = clean_sources.copy() #creates empty astropy.table object, where we can store our chosen stars 
    coordinates_chart = Table(names=('ID','real_xcentroid', 'real_ycentroid',
                                    'ghost_xcentroid', 'ghost_ycentroid'), dtype=('i4', 'f8', 'f8', 'f8', 'f8'))
    ghost_sources = setdiff(all_sources, clean_sources, keys=['xcentroid', 'ycentroid'])
    i = 0 
    for star in chosen_stars:
        """ using ghostfinder to find a star that has te closest values of x and y distances
        """
        real, ghost = ghostfinder(star, ghost_sources, **kwargs)
        if ghost == 'no_match':
            pass
        else:
            coordinates_chart.add_row([int(i), real['xcentroid'], real['ycentroid'],
                ghost['xcentroid'], ghost['ycentroid']])
            i += 1
    logging.info('Found %s ghosts', str(len(coordinates_chart)))
    logging.info('Finished')
    return(coordinates_chart)


def ghostfinder(star, ghost_sources, **kwargs):
    """ Used by ghostbuster to find the star and define which one is a ghost
    """
    distance = kwargs['distance']  #variable for cleaner coding
    right_x, right_y, = star['xcentroid'] + kwargs['right_x'], star['ycentroid'] + kwargs['right_y']
    for possible_ghost in ghost_sources:
        if abs(possible_ghost['xcentroid'] - right_x) < distance and abs(possible_ghost['ycentroid'] - right_y) < distance:
            return star, possible_ghost 
    return star, 'no_match'
    

def get_apertures(sources):
    logging.info('Creating apertures')
    positions_real = np.transpose((sources['real_xcentroid'], sources['real_ycentroid']))
    positions_ghost = np.transpose((sources['ghost_xcentroid'], sources['ghost_ycentroid']))
    for number, (pos_real, pos_ghost) in enumerate(zip(positions_real, positions_ghost)):
        apertures = CircularAperture([pos_real, pos_ghost], r=4)
        plt.text(pos_ghost[0]+5, pos_ghost[1]+5, number)
        plt.text(pos_real[0]+5, pos_real[1]+5, number)
        apertures.plot(color=np.random.rand(3,))
    
    logging.info('Plotting apertures for %s stars', str(len(sources)))
    logging.info('Finished')
    return apertures

def read_config():
    logging.info('Reading config.ini file')
    path = 'config.ini'
    config = configparser.ConfigParser()
    config.readfp(open('config.ini'))
    config_dict = {
            'fwhm': int(config['DEFAULT']['fwhm']),
            'threshold': int(config['DEFAULT']['threshold']),
            'sigma': int(config['DEFAULT']['sigma']),
            'min_star_limit': int(config['DEFAULT']['min_star_limit']),
            'norm_percent': float(config['DEFAULT']['norm_percent']),
            'plot_no_mask': config['EXECUTABLE_PARAMETERS']['plot_no_mask'],
            'plot_mask': config['EXECUTABLE_PARAMETERS']['plot_mask'],
            'distance': float(config['EXECUTABLE_PARAMETERS']['search_distance']),
            'square_size': int(config['EXECUTABLE_PARAMETERS']['mask_square_size']),
            'alpha': float(config['PLOT_PARAMETERS']['alpha']),
            'figsize': (config['PLOT_PARAMETERS']['figsize_rcparams_width'], config['PLOT_PARAMETERS']['figsize_rcparams_height']),
            'draw_aperture': config['PLOT_PARAMETERS']['draw_aperture'],
            'save_no_mask': config['PLOT_PARAMETERS']['save_plot_no_mask'],
            'save_mask': config['PLOT_PARAMETERS']['save_plot_mask'],
            'save_aperture': config['PLOT_PARAMETERS']['save_plot_aperture'],
            }
    if sys.argv[2] == 'p1':
        config_dict['savart_distance'] = config['SAVART_PARAMETERS_P1']['savart_distance'],
        config_dict['angle'] = config['SAVART_PARAMETERS_P1']['angle'],
        config_dict['right_x'] = float(config['SAVART_PARAMETERS_P1']['x_distance']),
        config_dict['right_y'] = float(config['SAVART_PARAMETERS_P1']['y_distance']),
    elif sys.argv[2] == 'p3':
        config_dict['right_x'] = float(config['SAVART_PARAMETERS_P3']['x_distance']),
        config_dict['angle'] = config['SAVART_PARAMETERS_P1']['angle'],
        config_dict['savart_distance'] = config['SAVART_PARAMETERS_P3']['savart_distance'],
        config_dict['right_y'] = float(config['SAVART_PARAMETERS_P3']['y_distance']),
    
    logging.info('Finished')
    return config_dict

def exec():
    logging.basicConfig(filename='logger.log', level=logging.DEBUG)
    logging.info('----------- %s -----------', date.today())
    file_path = sys.argv[1]
    if os.path.isfile(file_path):
        logging.info('Opening file: %s', str(os.path.basename(file_path)))
        print(os.path.basename(file_path))
    else:
        logging.ERROR('File does not exist')
        print('File does not exist')
        sys.exit()

    image_data, hdr = open_image(file_path)
    coo = get_hdr(hdr)
    file_name = os.path.basename(file_path)
    file_name = file_name.split('.')[0]
    config = read_config()
    all_sources, norm, mask = star_data_gatherer(image_data, create_mask(image_data), **config)
    if config['plot_no_mask'] == 'True':
        plot_photo(image_data, mask, norm, **config)
        logging.info('Ploted image without masks')
        if config['save_no_mask'] == True:
            plt.savefig('plots_no_mask/' + str(file_name.strip) + '.png')
    mask = ghostbuster_2(all_sources, mask, **config)
    clean_sources, norm, mask = star_data_gatherer(image_data, mask, **config)
    logging.info('Found real stars')
    if config['plot_mask'] == 'True':
        plot_photo(image_data, mask, norm, **config)
        logging.info('Plotted image with masked ghosts')
        if config['save_mask'] == 'True':
            plt.savefig('plots_mask/' + str(file_name) + '_mask.png') 
    ghost_xy = ghostbuster(all_sources, clean_sources, **config)
    if config['draw_aperture'] == 'True':
        get_apertures(ghost_xy)
        plot_photo(image_data,create_mask(image_data), norm, **config)
        if config['save_aperture'] == 'True':
            plt.savefig('plots_aperture/' + str(file_name) + '_ap.png')
    nova_path = sources_to_list(clean_sources, ghost_xy, file_path)
    wcs = run_astrometry(nova_path, coo)
    get_ra_dec_table(ghost_xy, wcs, file_path)
if __name__ == '__main__':
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)
    
    exec()








