bioformats2n5 converter
========================
_Important note_: This repository has been forked from https://github.com/glencoesoftware/bioformats2raw by Glencoe Software. 
We mainly restored the N5 support and to signal this we renamed the project. 


Java application to convert image file formats, including .mrxs,
to an intermediate N5 or Zarr structure.
The [raw2ometiff](https://github.com/glencoesoftware/raw2ometiff)
application can then be used to produce a
Bio-Formats 5.9.x ("Faas") or Bio-Formats 6.x (true OME-TIFF) pyramid.

Requirements
============

libblosc (https://github.com/Blosc/c-blosc) version 1.9.0 or later must be installed separately.
The native libraries are not packaged with any relevant jars.  See also note in jzarr readme (https://github.com/bcdev/jzarr/blob/master/README.md)

 * Mac OSX: `brew install c-blosc`
 * Ubuntu 18.04+: `apt-get install libblosc1`

Installation
============

1. Download and unpack a release artifact:

    https://github.com/glencoesoftware/bioformats2raw/releases

Development Installation
========================

1. Clone the repository:

    git clone git@github.com:keeneyetech/bioformats2n5.git

2. Run the Gradle build as required, a list of available tasks can be found by running:

    ./gradlew tasks

Usage
=====

Run the conversion:

    bioformats2n5 /path/to/file.mrxs /path/to/n5-pyramid --resolutions 6
    bioformats2n5 /path/to/file.svs /path/to/n5-pyramid --resolutions 6

Maximum tile dimensions are can be configured with the `--tile_width` and `--tile_height` options.  Defaults can be viewed with
`bioformats2n5 --help`.  `--resolutions` is optional; if omitted, the number of resolutions is set so that the smallest
resolution is no greater than 256x256.

If the input file has multiple series, a subset of the series can be converted by specifying a comma-separated list of indexes:

    bioformats2n5 /path/to/file.scn /path/to/n5-pyramid --series 0,2,3,4

By default, two additional readers (MiraxReader and PyramidTiffReader) are added to the beginning of Bio-Formats' list of reader classes.
Either or both of these readers can be excluded with the `--extra-readers` option:

    # only include the reader for .mrxs, exclude the reader for Faas pyramids
    bioformats2n5 /path/to/file.tiff /path/to/n5-pyramid --extra-readers MiraxReader
    # don't add any additional readers, just use the ones provided by Bio-Formats
    bioformats2n5 /path/to/file.mrxs /path/to/n5-pyramid --extra-readers

Reader-specific options can be specified using `--options`:

    bioformats2n5 /path/to/file.mrxs /path/to/n5-pyramid --options mirax.use_metadata_dimensions=false

Be aware when experimenting with different values for `--options` that the corresponding memo (cache) file may need to be
removed in order for new options to take effect.  This file will be e.g. `/path/to/.file.mrxs.bfmemo`.

The output in `/path/to/n5-pyramid` can be passed to `raw2ometiff` to produce
an OME-TIFF that can be opened in ImageJ, imported into OMERO, etc. See
https://github.com/glencoesoftware/raw2ometiff for more information.

Usage Changes
=============

Versions 0.2.6 and prior used the input file's dimension order to determine the output
dimension order, unless `--dimension-order` was specified.
Version 0.3.0 uses the `TCZYX` order by default, for compatibility with https://ngff.openmicroscopy.org/0.2/#image-layout.
The `--dimension-order` option can still be used to set a specific output dimension order, e.g.:

    bioformats2n5 /path/to/file.mrxs /path/to/n5-pyramid --dimension-order XYCZT

or can be set to use the input file's ordering, preserving the behavior of 0.2.6:

    bioformats2n5 /path/to/file.mrxs /path/to/n5-pyramid --dimension-order original

If a specific dimension order is passed to `--dimension-order`, it must be a valid dimension order as defined in
the [OME 2016-06 schema](https://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome_xsd.html#Pixels_DimensionOrder).
The specified dimension order is then reversed when creating Zarr arrays, e.g. `XYCZT` would become `TZCYX` in Zarr.

Performance
===========

This package is __highly__ sensitive to underlying hardware as well as
the following configuration options:

 * `--max_workers`
 * `--tile_width`
 * `--tile_height`

On systems with significant I/O bandwidth, particularly SATA or
NVMe based storage, you may find sharply diminishing returns with high
worker counts.  There are significant performance gains to be had utilizing
larger tile sizes but be mindful of the consequences on the downstream
workflow.

The worker count defaults to the number of detected CPUs.  This may or may not be appropriate for the chosen input data.
If reading a single tile from the input data requires a lot of memory, decreasing the worker count will be necessary
to prevent memory exhaustion.  JPEG, PNG, and certain TIFFs are especially susceptible to this problem.

The worker count should be set to 1 if the input data requires a Bio-Formats reader that is not thread-safe.
This is not a common case, but is a known issue with Imaris HDF data in particular.

In general, expect to need to tune the above settings and measure
relative performance.

License
=======

The converter is distributed under the terms of the GPL license.
Please see `LICENSE.txt` for further details.
