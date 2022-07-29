"""Dielectrics submodule

This submodule has functions that deal with moving between different ways of
representing the electromagnetic properties of dielectrics, i.e. conversions
between conductivity, loss tangent, and the imaginary part of the complex
relative permittivity.

There is also one function on calculating the equivalent real permittivity
of a multilayered material.

Finally, there are implementations of the Cole-Cole single pole relaxation
model and the Debye multipole one.
"""

from typing import List
import numpy as np
from scipy.constants import epsilon_0


np.seterr(divide='raise', invalid='raise')


def conductivity_to_tan_delta(freq: float, conductivity: float,
                              real_permittivity: float) -> float:
    """Converts between conductivity and loss tangent at a specific frequency

    This is a simple and straightforward conversion between the value of
    conductivity, in S/m, at a particular frequency, and a loss tangent.

    Args:
        freq: A `float` with the frequency, in GHz, at which to do the
              conversion
        conductivity: A `float` value for the conductivity in S/m
        real_permittivity: A `float` value for the real part of the
                           complex relative permittivity.

    Returns:
        The value for the loss tangent, as a `float` number.

    Raises:
        ValueError: If a negative value is provided for the permittivity or
                    the conductivity
        ZeroDivisionError: If you specify 0 Hz, i.e. DC, for the frequency
    """

    if (real_permittivity < 0) or (conductivity < 0):
        raise ValueError('The real part of the permittivity and the'
                         ' conductivity must be positive')

    try:
        tan_delta = 17.97591 * conductivity / (real_permittivity * freq)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Real part and frequency must be > 0'). \
              with_traceback(error.__traceback__)

    return tan_delta


def tan_delta_to_conductivity(freq: float, real_permittivity: float,
                              tan_delta: float) -> float:
    """Converts between loss tangent and conductivity

    This is a simple and straightforward conversion between the loss tangent
    and the value of conductivity, in S/m, at a particular frequency.

    Args:
        freq: A `float` with the frequency, in GHz, at which to do the
              conversion
        real_permittivity: A `float` value for the real part of the complex
                           relative permittivity.
        tan_delta: A `float` value for the loss tangent.

    Returns:
        The value for conductivity in S/m, which is equivalent to mho/m.

    Raises:
        ValueError: If you specify 0 Hz, i.e. DC, for the frequency
    """

    if real_permittivity < 0:
        raise ValueError('The real part of the permittivity must be positive')

    if np.isclose(freq, 0):
        raise ValueError('Frequency must be > 0')

    conductivity = 0.05563 * real_permittivity * tan_delta * freq

    return conductivity


def complex_permittivity_to_tan_delta(real_permittivity: float,
                                      imag_permittivity: float) -> float:
    """Computes loss tangent from complex relative permittivity

    This is a simple and straightforward calculation of a material's loss
    tangent from the real and imaginary parts of its complex relative
    permittivity.

    Args:
        real_permittivity: A `float` value for te real part of the
                           complex relative permittivity.
        imag_permittivity: A `float` value for the imaginary part of the
                           complex relative permittivity.

    Returns:
        The value for the loss tangent.

    Raises:
        ZeroDivisionError: If you specify 0 Hz, i.e. DC, for the real part
                           of the permittivity.
    """

    try:
        tan_delta = np.abs(imag_permittivity) / np.abs(real_permittivity)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Real part must be > 0'). \
              with_traceback(error.__traceback__)

    return tan_delta


def imaginary_permittivity_to_conductivity(freq: float,
                                           imag_permittivity: float) -> float:
    """Converts between imaginary permittivity and conductivity

    This is a simple and straightforward conversion between the imaginary
    part of the complex relative permittivity and the value of conductivity,
    in S/m, at a particular frequency.

    Args:
        freq: A `float` with the frequency, in GHz, at which to do the
              conversion
        imag_permittivity: A `float` value for the imaginary part of the
                           complex relative permittivity.

    Returns:
        The value for conductivity in S/m, which is equivalent to mho/m.

    Raises:
        ZeroDivisionError: If you specify 0 Hz, i.e. DC, for the frequency
    """

    if np.isclose(freq, 0):
        raise ZeroDivisionError('Frequency must be > 0')

    conductivity = 0.05563 * freq * np.abs(imag_permittivity)

    return conductivity


def conductivity_to_imaginary_permittivity(freq: float,
                                           conductivity: float) -> float:
    """Converts between conductivity and imaginary permittivity

    This is a simple and straightforward conversion between the value
    of conductivity, in S/m, and the imaginary part of the complex
    relative permittivity, at a particular frequency.

    Args:
        freq: A `float` with the frequency, in GHz, at which to do the
              conversion
        conductivity: A `float` value for the conductivity in S/m

    Returns:
        The unitless value for the imaginary part of the complex
        relative permittivity.

    Raises:
        ZeroDivisionError: If you specify 0 Hz, i.e. DC, for the frequency
    """

    try:
        imag_permittivity = 17.97591 * conductivity / freq
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    return imag_permittivity


def equivalent_relative_permittivity(epsilon_real: List[float],
                                     thicknesses: List[float]) -> float:
    """Calculate equivalent relative permittivity of multi-layer media

    Uses the low-frequency approximation for equivalent real part of the
    relative permittivity of a multi-layer dielectric media. This approximation
    is valid up to about 100 GHz.

    Args:
        epsilon_real: A `List` of `float` values for the real part of the
                      relative permittivity of the individual layers.
        thicknesses: A `List` of `float` values for the thickness of the
                     individual layers. Units do not matter as long as they
                     are all the same.

    Returns:
        A `float` with the equivalent relative permittivity.

    Raises:
        ZeroDivisionError: In case the total thickness is zero, or a
                           particular relative permittivity is given as
                           zero.
    """

    total_thickness = np.sum(thicknesses)
    epsilon_real_eff = 0.0

    try:
        for er_eff, ind_thickness in zip(epsilon_real, thicknesses):
            epsilon_real_eff += (ind_thickness / (er_eff * total_thickness))
        epsilon_real_eff = 1 / epsilon_real_eff
    except ZeroDivisionError as error:
        raise ZeroDivisionError('One or more arguments evaluate to zero'). \
              with_traceback(error.__traceback__)

    return epsilon_real_eff


def cole_cole_single(freq: float, er_static: float, er_inf: float,
                     cond_static: float, relax_time: float,
                     alpha: float) -> complex:
    """Single-pole Cole-Cole model

    This function implements the single-pole Cole-Cole dielectric relaxation
    model. This is often used to model the complex relative permittivity of
    various dielectric materials.

    Args:
        freq: A `float` with the frequency at which to calculate the complex
              relative permittivity. Units are GHz.
        er_static: A `float` with the material's relative permittivity at 0 Hz
        er_inf: A `float` with the material's relative permittivity at infinity
        cond_static: A `float` with the material's static electrical
                     conductivity. Units are S/m.
        relax_time: A `float` with the material's relaxation time. Units are
                    seconds.
        alpha: A `float` coefficient describing the pole broadening.

    Returns:
        A single `complex` number of the form `e_real - j * e_imag`.

    Raises:
        ZeroDivisionError: In case the frequency is given as 0 Hz.
    """

    freq *= 1e9
    ang_freq = 2 * np.pi * freq

    try:
        er_complex_1 = cond_static / (1j * ang_freq * epsilon_0)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0'). \
              with_traceback(error.__traceback__)

    er_complex_2 = er_static - er_inf
    temp_base = 1j * ang_freq * relax_time
    er_complex_2 /= (1 + np.power(temp_base, 1 - alpha))

    er_complex = er_inf + er_complex_2 + er_complex_1

    return er_complex


def debye_multipole(freq: float, er_inf: float, cond_static: float,
                    relax_times: List[float],
                    er_disps: List[float]) -> complex:
    """Multipole Debye model

    This function implements the multipole Debye dielectric relaxation
    model. This is often used to model the complex relative permittivity of
    various dielectric materials.

    Args:
        freq: A `float` with the frequency at which to calculate the complex
              relative permittivity. Units are GHz.
        er_inf: A `float` with the material's relative permittivity at infinity
        cond_static: A `float` with the material's static electrical
                     conductivity. Units are S/m.
        relax_times: A `List` of `float` with the material's relaxation times.
                     Units are seconds.
        er_disps: A `List` of `float` with the material's pole amplitudes.

    Returns:
        A single `complex` number of the form `e_real - j * e_imag`.

    Raises:
        RuntimeError: In case a different number of relaxation times and
                      pole amplitudes are given.
        ZeroDivisionError: In case the frequency is given as 0 Hz.
    """

    if len(relax_times) != len(er_disps):
        raise RuntimeError(
            'Need same number of relaxation times and pole amplitudes'
        )

    freq *= 1e9
    ang_freq = 2 * np.pi * freq

    try:
        er_complex_1 = cond_static / (1j * ang_freq * epsilon_0)
    except ZeroDivisionError as error:
        raise ZeroDivisionError('Frequency must be > 0.'). \
              with_traceback(error.__traceback__)

    er_complex_2 = 0 + 0j
    for t_relax, er_disp in zip(relax_times, er_disps):
        er_complex_2 += (er_disp / (1 + (1j * ang_freq * t_relax)))

    er_complex = er_inf + er_complex_2 + er_complex_1

    return er_complex
