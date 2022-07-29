"""Propagation submodule

A few functions vaguely related to electromagnetic propagation in conductors
and dielectrics.
"""

from typing import Tuple
import numpy as np
from scipy.constants import epsilon_0, mu_0, speed_of_light


np.seterr(divide='raise', invalid='raise')


def skin_depth(freq: float, conductivity: float,
               real_permeability: float) -> float:
    """Calculates skin depth for a particular metal at a particular frequency

    Uses the well-known formula for metal skin depth, with some additional
    error-checking to prevent runtime errors.

    Args:
        freq: A `float` with the frequency at which we want to know the skin
              depth. Units are GHz.
        conductivity: A `float` with the value for the metal's conductivity.
                      Units are S/m.
        real_permeability: A `float` with the relative permeability of the
                           metal. Unitless.

    Returns:
        A single `float` with the skin depth in metres, due to using base
        units in the function.

    Raises:
        RuntimeError: If for whatever reason one or more of the input
                      variables are negative.
        ZeroDivisionError: If for whatever reason one or more of the input
                           variables are zero.
    """

    # ? Add database of metals and their properties

    freq *= 1e9

    delta = np.pi * freq * conductivity * mu_0 * real_permeability
    delta = np.sqrt(delta)

    if np.isnan(delta):
        raise RuntimeError('All variables must be > 0')

    try:
        delta = 1 / delta
    except (ZeroDivisionError, FloatingPointError) as error:
        raise ZeroDivisionError('Variable values must be > 0'). \
              with_traceback(error.__traceback__)

    return delta


def metal_resistance(freq: float, conductivity: float,
                     real_permeability: float) -> float:
    """Calculates frequency-dependent metal resistance

    Uses the formula in the Waveguide Handbook by Marcuvitz to calculate the
    frequency-dependent resistance of a given metal.

    Args:
        freq: A `float` with the frequency at which we want to know the skin
              depth. Units are GHz.
        conductivity: A `float` with the value for the metal's conductivity.
                      Units are S/m.
        real_permeability: A `float` with the relative permeability of the
                           metal. Unitless.

    Returns:
        A single `float` with the frequency-dependent metal resistance in Ohms.

    Raises:
        ZeroDivisionError: If the `freq` is specified as zero.

    """

    metal_skin_depth = skin_depth(freq, conductivity, real_permeability)

    freq *= 1e9

    try:
        wavelength = speed_of_light / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    resistance = np.pi * np.sqrt(mu_0 / epsilon_0)
    resistance *= (metal_skin_depth / wavelength)

    return resistance


def plane_wave_prop_const(freq: float, real_permittivity: float,
                          imag_permittivity: float,
                          real_permeability: float) -> Tuple[float, float]:
    """Calculate the complex propagation constant in homogeneous medium

    The general-case formula is used to find the attenuation constant in Np/m
    and the phase constant in rad/m of a planar EM wave in homogeneous medium.

    Args:
        freq: A `float` with the frequency of interest. Units are GHz.
        real_permittivity: A `float` with the relative permittivity of the
                           medium. Unitless.
        imag_permittivity: A `float` with the imaginary part of the complex
                           relative pertmittivity of the medium. Unitless.
        real_permeability: A `float` with the relative permeability of the
                           medium. Unitless.

    Returns:
        A `tuple` consisting of `float` values for the attenuation constant
        `alpha` and the phase constant `beta`. Units are Np/m and rad/m.

    Raises:
        ZeroDivisionError: If the real part of the relative permittivity is
                           given as zero.
        RuntimeError: If any of the square root arguments turn out to be
                      negative
    """

    imag_permittivity = np.abs(imag_permittivity)

    real_permittivity *= epsilon_0
    imag_permittivity *= epsilon_0
    real_permeability *= mu_0

    freq *= 1e9
    ang_freq = 2 * np.pi * freq

    try:
        common_root = (
            1 + np.float_power(imag_permittivity / real_permittivity, 2)
        )
    except (ZeroDivisionError, FloatingPointError) as error:
        raise ZeroDivisionError('Real relative permittivity must be >= 1'). \
              with_traceback(error.__traceback__)

    common_root = np.sqrt(common_root)

    if np.isnan(common_root):
        raise RuntimeError('All variables must be > 0')

    common_multiplier = (real_permeability * real_permittivity) / 2.0

    alpha = ang_freq * np.sqrt(common_multiplier * (common_root - 1))

    if np.isnan(alpha):
        raise RuntimeError('All variables must be > 0')

    beta = ang_freq * np.sqrt(common_multiplier * (common_root + 1))

    if np.isnan(beta):
        raise RuntimeError('All variables must be > 0')

    return (alpha, beta)
