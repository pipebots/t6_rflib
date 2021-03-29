"""Conversions submodule

Short and sweet functions to go between dB, Np, and magnitude. Might be more
added in the future, when and if the need for those arises.
"""

import numpy as np


np.seterr(divide='raise', invalid='raise')


def nepers_to_db(nepers: float) -> float:
    return (nepers * 8.685889638)


def db_to_nepers(db: float) -> float:
    return (db * 0.115129255)


def db_to_mag(value: float, mode: str = 'power') -> float:
    """Convert from dB to magnitude

    Quick conversion between dB, i.e. logarithmic, and magnitude, i.e. linear
    values. Supports both power and amplitude conversions.

    Args:
        value: A `float` with the value to convert.
        mode: A `str` which determines whether to use 10 * log10 or 20 * log10.
              Must be either `power` or `amplitude`.

    Returns:
        The magnitude as a `float` number.

    Raises:
        ValueError: In case a mode different than `power` or `amplitude` is
                    specified.
    """

    if 'power' == mode.lower():
        mag = np.float_power(10, value / 10)
    elif 'amplitude' == mode.lower():
        mag = np.float_power(10, value / 20)
    else:
        raise ValueError('Mode must be power or amplitude')

    return mag


def mag_to_db(value: float, mode: str = 'power') -> float:
    """Convert from magnitude to dB

    Quick conversion between dB, i.e. logarithmic, and magnitude, i.e. linear
    values. Supports both power and amplitude conversions.

    Args:
        value: A `float` with the value to convert.
        mode: A `str` which determines whether to use 10 * log10 or 20 * log10.
              Must be either `power` or `amplitude`.

    Returns:
        The dB value as a `float` number.

    Raises:
        ValueError: In case a mode different than `power` or `amplitude` is
                    specified. Additionally raised if `value` is given as a
                    number that is <= 0.
    """

    if 'power' == mode.lower():
        mag = 10
    elif 'amplitude' == mode.lower():
        mag = 20
    else:
        raise ValueError('Mode must be power or amplitude')

    try:
        mag *= np.log10(value)
    except FloatingPointError as error:
        raise ValueError('Magnitude must be > 0'). \
              with_traceback(error.__traceback__)

    return mag
