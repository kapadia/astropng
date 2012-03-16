import numpy
import pyfits

def generate_fits():
    """
    Generate a FITS image with BITPIX = -32
    """
    DIMENSION = 4
    
    header = pyfits.Header()
    header.update('SIMPLE', 'T')
    header.update('BITPIX', -32)
    header.update('NAXIS', 2)
    header.update('NAXIS1', DIMENSION)
    header.update('NAXIS2', DIMENSION)
    
    # Generate random data
    data = numpy.random.random(DIMENSION * DIMENSION).astype(numpy.float32).reshape((DIMENSION, DIMENSION))
    
    hdu = pyfits.PrimaryHDU(data, header)
    hdu.writeto('test.fits')
    
generate_fits()