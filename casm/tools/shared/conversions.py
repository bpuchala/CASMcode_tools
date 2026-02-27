import copy

import numpy as np


def to_kelvin(value: np.ndarray) -> np.ndarray:
    R"""Convert 3x3 symmetric tensor to Kelvin notation.
    Parameters
    ----------
    value : np.ndarray
        The 3x3 symmetric tensor.

    Returns
    -------
    v_kelvin : np.ndarray
        The tensor in Kelvin notation,
        :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
        \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.
    """
    v_kelvin = np.array(
        [
            value[0, 0],
            value[1, 1],
            value[2, 2],
            np.sqrt(2.0) * value[1, 2],
            np.sqrt(2.0) * value[0, 2],
            np.sqrt(2.0) * value[0, 1],
        ]
    )
    return v_kelvin


def from_kelvin(v_kelvin: np.ndarray) -> np.ndarray:
    R"""Convert tensor from Kelvin notation to 3x3 symmetric tensor.

    Parameters
    ----------
    v_kelvin : np.ndarray
        The tensor in Kelvin notation,
        :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
        \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.

    Returns
    -------
    value : np.ndarray
        The 3x3 symmetric tensor.
    """
    value = np.array(
        [
            [v_kelvin[0], v_kelvin[5] / np.sqrt(2.0), v_kelvin[4] / np.sqrt(2.0)],
            [v_kelvin[5] / np.sqrt(2.0), v_kelvin[1], v_kelvin[3] / np.sqrt(2.0)],
            [v_kelvin[4] / np.sqrt(2.0), v_kelvin[3] / np.sqrt(2.0), v_kelvin[2]],
        ]
    )
    return value


def voigt_to_kelvin(v_voigt: np.ndarray) -> np.ndarray:
    R"""Convert tensor from Voigt notation to Kelvin notation.

    Parameters
    ----------
    v_voigt : np.ndarray
        The tensor in Voigt notation,
        :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz},
        \sigma_{xy}]`.

    Returns
    -------
    v_kelvin : np.ndarray
        The tensor in Kelvin notation,
        :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
        \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.
    """
    v_kelvin = copy.deepcopy(v_voigt)
    v_kelvin[3:] *= np.sqrt(2.0)
    return v_kelvin


def kelvin_to_voigt(v_kelvin: np.ndarray) -> np.ndarray:
    R"""Convert tensor from Kelvin notation to Voigt notation.

    Parameters
    ----------
    v_kelvin : np.ndarray
        The tensor in Kelvin notation,
        :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sqrt(2)\sigma_{yz},
        \sqrt(2)\sigma_{xz}, \sqrt(2)\sigma_{xy}]`.

    Returns
    -------
    v_voigt : np.ndarray
        The tensor in Voigt notation,
        :math:`[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz},
        \sigma_{xy}]`.
    """
    v_voigt = copy.deepcopy(v_kelvin)
    v_voigt[3:] /= np.sqrt(2.0)
    return v_voigt
