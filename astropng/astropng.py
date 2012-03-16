import os
import numpy
import pyfits
from AstroPNGWriter import AstroPNGWriter

class AstroPNG(object):
    """
    Yes, this is a class for representing astronomical data in the PNG format.
    We utilize the 16 bit integer format of PNGs, and custom chunks to store
    various metadata such as the FITS header.
    """
    
    def __init__(self, f):
        """        
        :param f:   The path of the image file
        """
        filename, extension = os.path.splitext(f)
        if extension.lower() == '.png':
            self.png = f
            self.fits = None
        else:
            self.fits = f
            self.png = None
        
        self.quantized = False
        self.clipped = False
    
    def to_png(self, out_file, bit_depth = 16, clip_on_percentiles = False):
        """
        Converts the FITS file to a PNG.
        
        :param out_file:            User specified filename for the PNG
        :param clip_on_percentiles: Clips flux values to a lower and upper percentile
        
        .. warning:: Setting bit_depth to 8 may reduce the dynamic range of the image.
        """
        if not self.fits:
            raise ValueError, "AstroPNG was not initialized with a FITS image"
        
        hdu = pyfits.open(self.fits)
        header, fluxes = hdu[0].header, hdu[0].data
        
        height, width = fluxes.shape
        
        # Prepare data
        fluxes = numpy.flipud(fluxes)
        fluxes = fluxes.flatten()
        
        # Determine the minimum pixel and maximum pixel value
        min_pix, max_pix = fluxes.min(), fluxes.max()
        
        # Clip data
        if clip_on_percentiles:
            min_pix, max_pix = self.__compute_percentile(fluxes)
            fluxes = self.__clip(fluxes, min_pix, max_pix)
    
        # Scale image to 8 bit integer space
        if bit_depth == 8:
            range_of_pixels = max_pix - min_pix
            fluxes = 255 * ( (data - min_pix) / range_of_pixels )
            min_pix, max_pix = 0, 255
        
        # Reshape the data to its original dimensions
        fluxes = fluxes.reshape(height, width)
        
        # If non-integer data type then quantize pixels
        if header['BITPIX'] in (-32, -64):
            nans = self.__find_nans(fluxes)
            z_zeros, z_scales, fluxes = self.__quantize(fluxes)
            
            # Replace the nans with zeros
            
            nan_indices = numpy.where(numpy.isnan(fluxes))
            fluxes[nan_indices] = 0
        
        # Create a PNG writer object with the appropriate settings
        png_writer = AstroPNGWriter(
            width = width,
            height = height,
            greyscale = True,
            alpha = False,
            bitdepth = bit_depth,
        )

        # Set various metadata in PNG
        png_writer.set_header(header)
        if self.quantized:
            png_writer.set_quantization_parameters(z_zeros, z_scales)
            png_writer.set_nans(nans)

        f = open(out_file, 'wb')
        png_writer.write(f, fluxes)
        f.close()
    
    def to_fits(self, out_file):
        """
        Converts a PNG to file a FITS file.
        
        :param out_file:    User specified filename for the FITS
        """
        raise NotImplementedError
    
    def __quantize(self, fluxes):
        """
        Quantizes the fluxes in the FITS image when it is represented by floats.
        
        .. warning:: Calling this method results in lossy compression
        """
        # Get the dimensions
        height, width = fluxes.shape
        
        # Generate 10000 random numbers
        random_numbers = self.__random_number_generator(N = width * height).reshape( (height, width) )
        
        # Get the zeros and scales of each row then quantize with dithering
        z_zeros     = numpy.nanmin(fluxes, axis = 1)
        z_scales    = self.__zscales(fluxes)
        quantized_data = numpy.round( ( fluxes - numpy.vstack(z_zeros) ) / numpy.vstack(z_scales) + random_numbers - 0.5 )
        
        self.quantized = True
        return z_zeros, z_scales, quantized_data
    
    def __find_nans(self, fluxes):
        """
        Locates the NANs in the image.  NANs are present in data represented by floats.
        Returns a 1D array of x, y coordinates.
        
        e.g. data[y, x] = nan
        """
        y, x = numpy.where(numpy.isnan(fluxes))
        locations = numpy.empty(x.size + y.size, dtype = numpy.int16)
        locations[0::2], locations[1::2] = x, y
        
        return locations
    
    def __median_absolute_deviations(self, fluxes):
        """
        Determines the noise level of each row in a 2D array of fluxes.  This is used when
        quantizing data.
        
        Stoehr, F. et al. 2007, ST-ECF Newsletter. 42, 4
        """
        n = len(fluxes[0])
        noise = 0.6052697 * numpy.median(numpy.abs(2.0 * fluxes[:, 2:n-2] - fluxes[:, 0:n-4] - fluxes[:, 4:n]), axis = 1)
        return noise

    def __zscales(self, fluxes, D = 100):
        """
        :param fluxes:  2D array of fluxes
        :param D:       Relates to noise bits per pixel.  Suggested range is
                        between (10, 100) which corresponds to 5 to 8 noise
                        bits per pixel.
                        
        Press, Chicago, and Astronomical Society. 2012.  Lossless Astronomical
        Image Compression and the Effects of Noise.  Astronomical Society of
        the Pacific 121 (878): 414-427.
        
        Computes the scale parameter that is used to represent 32 bit float as integers.
        """
        sigmas = self.__median_absolute_deviations(fluxes)
        return sigmas / float(D)

    def __random_number_generator(self, N = 10000):
        """
        Random number generator as described in:

        http://fits.gsfc.nasa.gov/registry/tilecompression/tilecompression2.2.pdf
        """
        a, m, seed = 16807.0, 2147483647.0, 1
        random_numbers = numpy.empty(N)

        for i in xrange(N):
            temp = a * seed
            seed = temp - m * numpy.floor(temp / m)
            random_numbers[i] = seed / m
        return random_numbers

    def __compute_percentile(self, fluxes, lower=0.0025, upper=0.9975):
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

        vmin, vmax = data[vmin_index], data[vmax_index]
        return vmin, vmax
    
    def __clip(self, fluxes, vmin, vmax):
        """
        :param data: 1D array
        :param vmin: Minimum value for pixels
        :param vmax: Maximum value for pixels

        Clip data on vmin and vmax.
        """
        min_indexes = numpy.where( data < vmin )[0]
        max_indexes = numpy.where( data > vmax )[0]
        data[min_indexes] = vmin
        data[max_indexes] = vmax

        self.clipped = True
        return data
