from distutils.core import setup

setup(
    name='astropng',
    version='0.0.1',
    description='Python package to convert FITS images to 8 or 16 bit PNGs',
    author='Amit Kapadia',
    author_email='amit@zooniverse.org',
    url='https://github.com/kapadia/astropng',
    requires=['png', 'pyfits'],
    packages=['astropng']
)