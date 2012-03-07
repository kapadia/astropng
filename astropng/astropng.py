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
    
    
def compute_percentile(data, lower=0.0025, upper=0.9975):
    """
    :param data:    Array of pixel intensities
    :param lower:   Lower percentile
    :param upper:   Upper percentile

    Computes the pixel value of data at the given percentiles
    """

    data = numpy.sort(data)
    num_elements = data.size
    vmin_index = numpy.floor( lower*(num_elements-1) + 1 )
    vmax_index = numpy.floor( upper*(num_elements-1) + 1 )
    lower_index = numpy.floor( (lower - 0.0015)*(num_elements - 1) + 1 );
    upper_index = numpy.floor( (upper + 0.0015)*(num_elements - 1) + 1 );

    vmin = data[vmin_index]
    vmax = data[vmax_index]

    return vmin, vmax


def trim_data(data, vmin, vmax):
    """
    :param data: 1D array
    :param vmin: Minimum value for pixels
    :param vmax: Maximum value for pixels

    Trim data based on vmin and vmax.
    """
    min_indexes = numpy.where( data < vmin )[0]
    max_indexes = numpy.where( data > vmax )[0]
    data[min_indexes] = vmin
    data[max_indexes] = vmax

    return data
