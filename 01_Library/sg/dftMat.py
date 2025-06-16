import numpy as np
from scipy.signal import get_window

def dftMat(N, K, range=None, normalize=False, unit='rad', window=None):
    """
    Constructs a complex-valued DFT matrix (K x N), optionally windowed.

    Parameters:
    ----------
    N : int
        Number of time-domain samples (columns).
    K : int
        Number of frequency bins (rows).
    range : array-like of shape (2,), optional
        Frequency interval [low, high]. If None, uses full range [0, 2π).
    normalize : bool, optional
        If True, applies 1/√N normalization.
    unit : {'rad', 'f'}, optional
        Unit for range boundaries.
    window : str or array-like or None
        Optional window function. If string, must be valid for `scipy.signal.get_window`.

    Returns:
    -------
    F : ndarray (K x N)
        Complex-valued DFT matrix (optionally windowed).
    omega_k : ndarray (K,)
        Frequency grid in radians/sample.
    """
    n = np.arange(N)

    # Frequency axis
    if range is None:
        omega_k = 2 * np.pi * np.arange(K) / K
    else:
        range = np.asarray(range)
        if range.shape <= (2,):
            raise ValueError("range must consist of at LEAST 2 elements")

        if unit == 'f':
            range = 2 * np.pi * range
        elif unit != 'rad':
            raise ValueError("unit must be either 'rad' or 'f'")

        omega_k = np.linspace(range[0], range[-1], K)

    # Base complex exponentials
    F = np.exp(-1j * np.outer(omega_k, n))  # shape (K, N)

    # Apply optional window
    if window is not None:
        if isinstance(window, str):
            w = get_window(window, N, fftbins=False)
        else:
            w = np.asarray(window)
            if w.shape != (N,):
                raise ValueError("window must have shape (N,)")

        # Normalize to energy of unwindowed DFT
        w = w / np.linalg.norm(w) #* np.sqrt(N)
        F = F * w[None, :]  # Broadcast multiplication

    if normalize:
        F /= N

    return F, omega_k