astropng
========

Yes, this is a Python library for storing astronomical data in the PNG format.  We take advantage of the supported bit depths of PNGs, which can store either 8 bit or 16 bit integers.  For the cases when a FITS image is stored as a floating point number, we quantize the data and store a zero point and scale factor in a custom PNG chunk.  We also preserve the FITS header in another custom chunk.

Dependencies
------------
Numpy (http://numpy.scipy.org/)
PyFITS (http://www.stsci.edu/institute/software_hardware/pyfits)
PyPNG (http://code.google.com/p/pypng/)

Installation
------------
1. Download the source code from the repository (https://github.com/kapadia/astropng)

    git clone git@github.com:kapadia/astropng.git
    
2. Navigate to the root directory of the library
3. Run setup.py

    python setup.py install

Instructions
------------

    import astropng
    
    # To convert from a FITS image to a PNG
    ap = astropng.AstroPNG('blah.fits')
    ap.to_png('blah.png')
    
    # To convert from a PNG back to a FITS image
    ap = astropng.AstroPNG('blah.png')
    ap.to_fits('blah.fits')
    
References
----------
Stoehr, F. et al. 2007, ST-ECF Newsletter. 42, 4.

White, Richard L, Perry Greenfield, William Pence, Nasa Gsfc, Doug Tody, and Rob Seaman. 2011.  Tiled Image Convention for Storing Compressed Images in FITS Binary Tables: 1-17.