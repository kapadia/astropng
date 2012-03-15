import os
import numpy
import pyfits
from AstroPNGWriter import AstroPNGWriter

def fits2png(filename, bit_depth = 16, clip_on_percentiles = False):
    """
    :param filename:            Path to FITS image
    :param bit_depth:           Either 8 or 16 bit
    :param clip_on_quantiles:   Clip data based on computed quantiles
    
    Converts a FITS image into a PNG, while testing if the image may be represented by integers.
    """
    
    hdu = pyfits.open(filename)
    
    header = hdu[0].header
    data = hdu[0].data
    
    y_dim, x_dim = data.shape
    
    # Prep data for packing
    data = numpy.flipud(data)
    data = data.flatten()
    
    # Determine the minimum pixel and maximum pixel value
    min_pix = data.min()
    max_pix = data.max()
    
    # Calculate vmin and vmax
    vmin, vmax = compute_percentile(data)
    
    # Clip data
    if clip_on_percentiles:
        data = trim_data(data, vmin, vmax)
        min_pix = vmin
        max_pix = vmax
        
    # Scale image to 8 bit integer space
    if bit_depth == 8:
        range_of_pixels = max_pix - min_pix
        data = 255 * ( (data - min_pix) / range_of_pixels )
        min_pix, max_pix = 0, 255
    
    # Reshape the data to its original dimensions
    data = data.reshape(y_dim, x_dim)
    
    # Create a PNG writer object with the appropriate settings
    png_writer = AstroPNGWriter(
        width = x_dim,
        height = y_dim,
        greyscale = True,
        alpha = False,
        bitdepth = bit_depth,
    )
    
    # Set metadata in PNG
    png_writer.set_header(header)
    
    # Create the file name
    out_file = os.path.splitext(filename)[0] + '.png'
    
    f = open(out_file, 'wb')
    png_writer.write(f, data)
    f.close()
    
    
