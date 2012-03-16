import numpy

def to_scaled_integers(data):
    """
    :param data:    2D array of pixel intensities
    
    Scales a 2D float array of fluxes to integers,
    returning the array, an array of scale factors,
    and an array with locations of nans.
    """
    # Get the dimensions
    height, width = data.shape
    
    # Generate 10000 random numbers
    random_numbers = random_number_generator(N_RANDOM = width * height).reshape( (height, width) )
    
    # Get the zeros and scales of each row
    z_zeros     = numpy.nanmin(data, axis = 1)
    z_scales    = zscales(data)
    
    # Quantize the data
    quantized_data = numpy.round( ( data - numpy.vstack(z_zeros) ) / numpy.vstack(z_scales) + random_numbers - 0.5 )
    
    return z_zeros, z_scales, quantized_data

def nan_locations(data):
    """
    Determines the location of nans in a 2D array, and represent them with a 1D array of x, y coordinates.
    
    e.g. data[y, x] = nan
    """
    y, x = numpy.where(numpy.isnan(data))
    locations = numpy.empty(x.size + y.size, dtype = numpy.int16)
    locations[0::2], locations[1::2] = x, y
    
    return locations
      
def median_absolute_deviations(fluxes):
    """
    Determines the noise level of each row in a 2D array of fluxes.
    """
    n = len(fluxes[0])
    noise = 0.6052697 * numpy.median(numpy.abs(2.0 * fluxes[:, 2:n-2] - fluxes[:, 0:n-4] - fluxes[:, 4:n]), axis = 1)
    return noise
    
def zscales(data, D = 100):
    """
    Computes the scale parameter that is used to represent 32 bit float as integers.
    """
    sigmas = median_absolute_deviations(data)
    return sigmas / float(D)

def random_number_generator(N_RANDOM = 10000):
    """
    Random number generator as described in:
    
    http://fits.gsfc.nasa.gov/registry/tilecompression/tilecompression2.2.pdf
    """
    a, m, seed = 16807.0, 2147483647.0, 1
    random_numbers = numpy.empty(N_RANDOM)
    
    for i in xrange(N_RANDOM):
        temp = a * seed
        seed = temp - m * numpy.floor(temp / m)
        random_numbers[i] = seed / m
    return random_numbers
    
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
