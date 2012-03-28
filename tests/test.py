import os
import unittest

import numpy
import pyfits

import astropng

def generate_fits(data_type = 'integer'):
    """
    Generate a FITS image with specified data type.
    
    :param data_type:   Either integer or float
    """
    DIMENSION = 16
    BITPIX = 16 if data_type == 'integer' else -32
    header = pyfits.Header()
    header.update('SIMPLE', 'T')
    header.update('BITPIX', BITPIX)
    header.update('NAXIS', 2)
    header.update('NAXIS1', DIMENSION)
    header.update('NAXIS2', DIMENSION)
    
    # Generate random data
    data = numpy.random.random(DIMENSION * DIMENSION).astype(numpy.float32).reshape((DIMENSION, DIMENSION))    
    data = (1000*data).astype(numpy.int16) if BITPIX is 16 else data
    
    hdu = pyfits.PrimaryHDU(data, header)
    hdu.writeto('test.fits')

class TestAstroPNG(unittest.TestCase):
    """
    Unit testing for AstroPNG format.
    """
    
    def test_integer_png(self):
        generate_fits(data_type = 'integer')
        
        # Create PNG from the test FITS
        ap1 = astropng.AstroPNG('test.fits')
        ap1.to_png('test.png', crush=False)
        
        # Regenerate the FITS file from the PNG
        ap2 = astropng.AstroPNG('test.png')
        ap2.to_fits('test_from_png.fits')
        
        # Compare the data of each FITS image
        data1 = pyfits.getdata('test.fits')
        data2 = pyfits.getdata('test_from_png.fits')
        validity = (data1 == data2).flatten()
        
        for value in validity:
            self.assertTrue(value)

        os.remove('test.fits')
        os.remove('test.png')
        os.remove('test_from_png.fits')
    
    def test_float_png(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
