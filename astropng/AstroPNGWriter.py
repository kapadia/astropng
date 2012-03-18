import zlib
import struct
from array import array
from png import Writer
from png import write_chunk

import numpy

try:  # see :pyver:old
    array.tostring
except:
    def tostring(row):
        l = len(row)
        return struct.pack('%dB' % l, *row)
else:
    def tostring(row):
        """Convert row of bytes to string.  Expects `row` to be an
        ``array``.
        """
        return row.tostring()

_signature = struct.pack('8B', 137, 80, 78, 71, 13, 10, 26, 10)

class AstroPNGWriter(Writer):
    """
    Custom PNG writer to incorporate non-standard chunks for storing additional
    metadata from a FITS header metadata.  We implement 3 custom chunks:
    
    * fITS
        Stores the entire FITS header as an ASCII string.
    * qANT
        Stores 2 arrays representing ZZERO and ZSCALE when fluxes have been 
        quantized to integers.
    * nANS
        Stores the pixel coordinates of nans
    """
    
    def __init__(self, *args, **kwargs):
        Writer.__init__(self, *args, **kwargs)
        self.fITS = False
        self.qANT = False
        self.nANS = False

    def set_header(self, header):
        """
        Store metadata from the FITS header into the PNG chunk zOON
        """
        self.fits_header = str(header)
        self.fITS = True
    
    def set_quantization_parameters(self, zzero, zscale):
        """
        Store quantization parameters, ZZERO and ZSCALE.
        """
        self.zzero = zzero
        self.zscale = zscale
        self.nan_representation = None
        
        values = numpy.concatenate((zzero, zscale))
        if numpy.isnan(values).any():
            self.nan_representation = 0
            while True:
                if numpy.where( values == self.nan_representation)[0].size > 0:
                    self.nan_representation += 1
                else:
                    break
        nan_indexes = numpy.where(numpy.isnan(self.zzero))[0]
        self.zzero[nan_indexes] = self.nan_representation
        nan_indexes = numpy.where(numpy.isnan(self.zscale))[0]
        self.zscale[nan_indexes] = self.nan_representation
        
        self.qANT = True
    
    def set_nans(self, nans):
        """
        Store the location of NANs when quantizing float data.
        """
        self.nans = nans
        self.nANS = True
    
    def write_passes(self, outfile, rows, packed=False):
        """
        Write a PNG image to the output file.

        Most users are expected to find the :meth:`write` or
        :meth:`write_array` method more convenient.

        The rows should be given to this method in the order that
        they appear in the output file.  For straightlaced images,
        this is the usual top to bottom ordering, but for interlaced
        images the rows should have already been interlaced before
        passing them to this function.

        `rows` should be an iterable that yields each row.  When
        `packed` is ``False`` the rows should be in boxed row flat pixel
        format; when `packed` is ``True`` each row should be a packed
        sequence of bytes.

        """

        # http://www.w3.org/TR/PNG/#5PNG-file-signature
        outfile.write(_signature)

        # http://www.w3.org/TR/PNG/#11IHDR
        write_chunk(outfile, 'IHDR',
                    struct.pack("!2I5B", self.width, self.height,
                                self.bitdepth, self.color_type,
                                0, 0, self.interlace))

        # See :chunk:order
        # http://www.w3.org/TR/PNG/#11gAMA
        if self.gamma is not None:
            write_chunk(outfile, 'gAMA',
                        struct.pack("!L", int(round(self.gamma*1e5))))

        # See :chunk:order
        # http://www.w3.org/TR/PNG/#11sBIT
        if self.rescale:
            write_chunk(outfile, 'sBIT',
                struct.pack('%dB' % self.planes,
                            *[self.rescale[0]]*self.planes))

        # :chunk:order: Without a palette (PLTE chunk), ordering is
        # relatively relaxed.  With one, gAMA chunk must precede PLTE
        # chunk which must precede tRNS and bKGD.
        # See http://www.w3.org/TR/PNG/#5ChunkOrdering
        if self.palette:
            p,t = self.make_palette()
            write_chunk(outfile, 'PLTE', p)
            if t:
                # tRNS chunk is optional.  Only needed if palette entries
                # have alpha.
                write_chunk(outfile, 'tRNS', t)

        # http://www.w3.org/TR/PNG/#11tRNS
        if self.transparent is not None:
            if self.greyscale:
                write_chunk(outfile, 'tRNS',
                            struct.pack("!1H", *self.transparent))
            else:
                write_chunk(outfile, 'tRNS',
                            struct.pack("!3H", *self.transparent))

        # http://www.w3.org/TR/PNG/#11bKGD
        if self.background is not None:
            if self.greyscale:
                write_chunk(outfile, 'bKGD',
                            struct.pack("!1H", *self.background))
            else:
                write_chunk(outfile, 'bKGD',
                            struct.pack("!3H", *self.background))

        # Write custom chunk with information from the FITS header
        if self.fITS:
            length = len(self.fits_header)
            chunk_data = struct.pack("!%dc" % length, *self.fits_header)
            write_chunk(outfile, 'fITS', chunk_data)
            
        # Write custom chunk with quantization parameters
        if self.qANT:
            parameters = numpy.empty(self.zzero.size + self.zscale.size, dtype = self.zzero.dtype)
            parameters[0::2] = self.zzero
            parameters[1::2] = self.zscale
            length = parameters.size
            
            chunk_data = struct.pack("!I", self.nan_representation)
            chunk_data += struct.pack("!%df" % length, *parameters)
            write_chunk(outfile, 'qANT', chunk_data)
        
        # Write custom chunk with nan locations
        if self.nANS:
            nans = numpy.concatenate(self.nans)
            length = nans.size
            chunk_data = struct.pack("!%dI" % length, *nans)
            write_chunk(outfile, 'nANS', chunk_data)

        # http://www.w3.org/TR/PNG/#11IDAT
        if self.compression is not None:
            compressor = zlib.compressobj(self.compression)
        else:
            compressor = zlib.compressobj()

        # Choose an extend function based on the bitdepth.  The extend
        # function packs/decomposes the pixel values into bytes and
        # stuffs them onto the data array.
        data = array('B')
        if self.bitdepth == 8 or packed:
            extend = data.extend
        elif self.bitdepth == 16:
            # Decompose into bytes
            def extend(sl):
                fmt = '!%dH' % len(sl)
                data.extend(array('B', struct.pack(fmt, *sl)))
        else:
            # Pack into bytes
            assert self.bitdepth < 8
            # samples per byte
            spb = int(8/self.bitdepth)
            def extend(sl):
                a = array('B', sl)
                # Adding padding bytes so we can group into a whole
                # number of spb-tuples.
                l = float(len(a))
                extra = math.ceil(l / float(spb))*spb - l
                a.extend([0]*int(extra))
                # Pack into bytes
                l = group(a, spb)
                l = map(lambda e: reduce(lambda x,y:
                                           (x << self.bitdepth) + y, e), l)
                data.extend(l)
        if self.rescale:
            oldextend = extend
            factor = \
              float(2**self.rescale[1]-1) / float(2**self.rescale[0]-1)
            def extend(sl):
                oldextend(map(lambda x: int(round(factor*x)), sl))

        # Build the first row, testing mostly to see if we need to
        # changed the extend function to cope with NumPy integer types
        # (they cause our ordinary definition of extend to fail, so we
        # wrap it).  See
        # http://code.google.com/p/pypng/issues/detail?id=44
        enumrows = enumerate(rows)
        del rows

        # First row's filter type.
        data.append(0)
        # :todo: Certain exceptions in the call to ``.next()`` or the
        # following try would indicate no row data supplied.
        # Should catch.
        i,row = enumrows.next()
        try:
            # If this fails...
            extend(row)
        except:
            # ... try a version that converts the values to int first.
            # Not only does this work for the (slightly broken) NumPy
            # types, there are probably lots of other, unknown, "nearly"
            # int types it works for.
            def wrapmapint(f):
                return lambda sl: f(map(int, sl))
            extend = wrapmapint(extend)
            del wrapmapint
            extend(row)

        for i,row in enumrows:
            # Add "None" filter type.  Currently, it's essential that
            # this filter type be used for every scanline as we do not
            # mark the first row of a reduced pass image; that means we
            # could accidentally compute the wrong filtered scanline if
            # we used "up", "average", or "paeth" on such a line.
            data.append(0)
            extend(row)
            if len(data) > self.chunk_limit:
                compressed = compressor.compress(tostring(data))
                if len(compressed):
                    # print >> sys.stderr, len(data), len(compressed)
                    write_chunk(outfile, 'IDAT', compressed)
                # Because of our very witty definition of ``extend``,
                # above, we must re-use the same ``data`` object.  Hence
                # we use ``del`` to empty this one, rather than create a
                # fresh one (which would be my natural FP instinct).
                del data[:]
        if len(data):
            compressed = compressor.compress(tostring(data))
        else:
            compressed = ''
        flushed = compressor.flush()
        if len(compressed) or len(flushed):
            # print >> sys.stderr, len(data), len(compressed), len(flushed)
            write_chunk(outfile, 'IDAT', compressed + flushed)
        # http://www.w3.org/TR/PNG/#11IEND
        write_chunk(outfile, 'IEND')
        return i+1

