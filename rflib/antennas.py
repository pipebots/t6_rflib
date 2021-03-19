"""[summary]
"""

import numpy as np
import warnings
from scipy.constants import speed_of_light


np.seterr(divide='raise', invalid='raise')


def hertzian_dipole_current(freq: float, power: float, length: float,
                            units: str = 'dBm') -> float:
    """Convert Hertzian dipole radiated power to current

    This is a short helper function for gprMax that converts desired radiated
    power by a Hertzian dipole into the current amplitude that is required
    to be passed through the dipole.

    Notes:
        1. Strictly speaking, this formula is valid for a time-harmonic current
        excitation. Other, stranger excitations require more maths.

    Args:
        freq: A `float` with the frequency at which the dipole is radiating.
              Units are GHz.
        power: A `float` with the desired radiated power. Can be in either
               dBm or W units.
        length: A `float` with the physical length of the Hertzian dipole.
                Units are metres.
        units: A `str` with the units for `power`. Supported ones are 'W' and
               'dBm', with the latter converted to the former internally.

    Returns:
        A `float` value with the required current, in A, to be passed through
        the dipole to result in the desired radiated power.

    Raises:
        RuntimeWarning: If the physical length of the dipole is not much
                        smaller than the wavelength.
        ZeroDivisionError: If `power` or `freq` evaluate to zero.
        RuntimeError: If the `dipole_current` evaluates to a negative number.
    """

    if 'dbm' == units.lower():
        power = np.float_power(10, power / 10)
        power /= 1e3
    elif 'w' == units.lower():
        if 0 == power:
            raise ZeroDivisionError('Power in absolute units must be > 0')
    else:
        raise RuntimeError('Unsupported power units')

    freq *= 1e9

    try:
        wavelength = speed_of_light / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    if length > (wavelength / 10):
        warnings.warn('Dipole is not electrically small',
                      category=RuntimeWarning)

    try:
        dipole_current = 40 * np.float_power(np.pi, 2) * \
                         np.float_power(length / wavelength, 2)
        dipole_current = power / dipole_current
    except (ZeroDivisionError, FloatingPointError) as error:
        raise ZeroDivisionError('Dipole length must be > 0'). \
              with_traceback(error.__traceback__)

    dipole_current = np.sqrt(dipole_current)

    if np.isnan(dipole_current):
        raise RuntimeError('Dipole current somehow ended negative')

    return dipole_current


def max_antenna_separation_full(freq: float, radius: float,
                                mode: str = 'normal') -> float:
    """Maximum antenna separation for a given 1st Fresnel Zone

    This function calculates the maximum allowable distance between two
    antennas such that the radius of the first Fresnel zone is less than or
    equal to a specified value.

    Notes:
        1. This is derived from the full equation for the first Fresnel zone
        radius, as opposed to using the approximate formula for large antenna
        separation.
        2. Depending on the combination of input values you might get negative
        separation. While mathematically correct, this does not have any
        physical meaning.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        radius: A `float` with the maximum desired radius of the first
                Fresnel zone. Units are metres.
        mode: A `str` specifying the mode, either 'normal', i.e. we want the
              zone to be 100% free, or 'cheeky', for 60% of the zone.

    Returns:
        A single `float` number with the maximum separation between the two
        antennas, with units in metres.

    Raises:
        ZeroDivisionError: If the frequency has been given as zero.
    """

    freq *= 1e9

    try:
        wavelength = speed_of_light / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    # * Assume the pipe radius is 60% of the first Fresnel zone
    # * as opposed to 100%
    if 'cheeky' == mode.lower():
        radius /= 0.6

    separation = 16 * np.float_power(radius, 2) - np.float_power(wavelength, 2)
    separation /= (4 * wavelength)

    return separation


def max_antenna_separation_approx(freq: float, radius: float,
                                  mode: str = 'normal') -> float:
    """Maximum antenna separation for a given 1st Fresnel Zone

    This function calculates the maximum allowable distance between two
    antennas such that the radius of the first Fresnel zone is less than or
    equal to a specified value.

    Notes:
        1. This is derived from the approximate equation for the first Fresnel
        zone radius, which is the better-known and used one. However it assumes
        that the distance between transmitter and receiver is much, much larger
        than the wavelength, e.g. kilometres vs centimetres.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        radius: A `float` with the maximum desired radius of the first
                Fresnel zone. Units are metres.
        mode: A `str` specifying the mode, either 'normal', i.e. we want the
              zone to be 100% free, or 'cheeky', for 60% of the zone.

    Returns:
        A single `float` number with the maximum separation between the two
        antennas, with units in metres.

    Raises:
        ZeroDivisionError: If the frequency has been given as zero.
    """

    freq *= 1e9

    try:
        wavelength = speed_of_light / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    # * Assume the pipe radius is 60% of the first Fresnel zone
    # * as opposed to 100%
    if 'cheeky' == mode.lower():
        radius /= 0.6

    separation = 4 * np.float_power(radius, 2)
    separation /= wavelength

    return separation


def fresnel_zone_radius(freq: float, distance_1: float,
                        distance_2: float) -> float:
    """Calculates the radius of the 1st Fresnel zone

    Uses the well-known and used formula to find the radius of the
    first Fresnel zone at a specified point between two antennas.

    Notes:
        1. The assumption is that the distances between the antennas
        and the point of interest are much, much larger than the
        wavelength.
        2. The distances do not have to be the same, i.e. this is valid
        for any point along the length of the wireless link.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        distance_1: A `float` with the distance from the first antenna
                    to the point of interest. Units are in metres.
        distance_2: A `float` with the distance from the second antenna
                    to the point of interest. Units are in metres.

    Returns:
        A single `float` with the radius of the 1st Fresnel zone.
        Units are metres.

    Raises:
        ZeroDivisionError: In case the frequency or both distances are
                           given as zero.
    """

    freq *= 1e9

    try:
        wavelength = speed_of_light / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    try:
        radius = wavelength * distance_1 * distance_2
        radius /= (distance_1 + distance_2)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Distances must be > 0'). \
              with_traceback(error.__traceback__)

    radius = np.sqrt(radius)

    return radius


def far_field_distance(freq: float, antenna_dimension: float,
                       antenna_type: str = 'array') -> float:
    """Calculates far field boundary for an antenna

    Uses the well-established formula for far field region based on the
    largest physical size of an antenna.

    Notes:
        1. There is support for three common antenna types - monopole, dipole,
        and antenna array with half-wavelength spacing.
        2. The default behaviour assumes an array in which case the antenna
        dimension is the number of elements along x and/or y.
        3. In case of a monopole or dipole the antenna dimension is ignored as
        the antenna size is calculated from the free-space wavelength.
        4. Otherwise the units for the antenna dimension should be metres.

    Args:
        freq: A `float` with the frequency at which the far field distance
              is being calculated. Units are GHz.
        antenna_dimension : A `float` with the largest physical size of the
                            antenna. See Notes for further information
                            on units.
        antenna_type : A `str` with the type of antenna. Can be `monopole`,
                       `dipole`, or `array`. Any of these overrides the
                       `antenna_dimension` variable.

    Returns:
        A single `float` with the minimum far field distance, with units in m.

    Raises:
        ZeroDivisionError: If a frequency of zero is given.
    """

    freq *= 1e9

    try:
        wavelength = speed_of_light / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    # * A monopole is assumed to be a quarter-wavelength resonator, and
    # * a dipole to be a half-wavelength resonator
    if 'monopole' == antenna_type.lower():
        dimension = wavelength / 4.0
    elif 'dipole' == antenna_type.lower():
        dimension = wavelength / 2.0
    elif 'array' == antenna_type.lower():
        dimension = antenna_dimension * (wavelength / 2)
    else:
        dimension = antenna_dimension

    distance = 2 * np.float_power(dimension, 2)
    distance /= wavelength

    return distance
