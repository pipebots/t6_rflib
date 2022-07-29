```
 __________.__             ___.     |__|  __
 \______   \__|_____   ____\_ |__   _||__/  |_  ______ (C) George Jackson-Mills 2020
  |     ___/  \____ \_/ __ \| __ \ /  _ \   __\/  ___/
  |    |   |  |  |_> >  ___/| \_\ (  O_O )  |  \___ \
  |____|   |__|   __/ \___  >___  /\____/|__| /____  >
              |__|        \/    \/                 \/
```

# rflib - RF and Microwave engineering related utility functions

## Overview

This module contains various utility functions that are used in Theme 6's work on the electromagnetic modelling of sewer pipes. A lot of them are general purpose and can be applied to other situations. The functions themselves are not particularly optimised for performance at the moment.

There are currently four submodules:

- `antennas` - Currently has a few basic functions for calculating Fresnel zones, far-field distances, and current through a Hertzian dipole. Potential to merge work on patch and horn antenna design into here.
- `conversions` - Helper functions for moving between different units popular in electromagnetics, e.g. linear magnitude to dB, dB to Np, etc.
- `dielectrics` - Functions to deal with the different representations of a material's complex relative permittivity. Used extensively in setting up and running gprMax simulations and theoretical analysis of lossy waveguides.
- `propagation` - Functions related to electromagnetic wave propagation in various media.

## Requirements

Requires some standard Python packages for scientific computing. As my current development environment is a bit polluted, I'll list the required packages here. Apologies.

- `Python>=3.6`
- `numpy`
- `scipy`

## Tests

To be added. Bad practice, I know.

## Installation

Use `pip install -e .` in the folder to which you clone or download this. This will install `rflib` as an "editable" package in your current environment, meaning you should just do a `git pull` in the future to get any updates.

## Contributing

Contributions are more than welcome and are in fact actively sought! Please contact Viktor at [eenvdo@leeds.ac.uk](mailto:eenvdo@leeds.ac.uk).

## Acknowledgements

This work is supported by the UK's Engineering and Physical Sciences Research Council (EPSRC) Programme Grant EP/S016813/1
