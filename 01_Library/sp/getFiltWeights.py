import numpy as np
from scipy.interpolate import interp1d

def getFiltWeights(vW, vK, mRange=None):
    """
    Interpolates the given vector vW to a new resolution vK over specified frequency range(s).
    
    Parameters:
    -----------
    vW : array-like
        Original vector with weights (e.g., filter weights in frequency domain).
    vK : int
        Number of interpolated values in the target range.
    mRange : array-like, optional
        Frequency range(s) [start, end] for interpolation.
        Can be 1D (single range) or 2D (multiple ranges).
    
    Returns:
    --------
    vW_D : ndarray
        Interpolated vector (or stacked vectors for multiple ranges).
    """
    
    sKOrig = len(vW)  # Original length of vW
    sRes = np.pi / (sKOrig // 2)  # Frequency resolution for [0, pi)
    
    # Generate frequency axis vOmegaVw for the upper half of the spectrum (Nyquist range)
    vOmegaVw = np.linspace(0, (np.pi - sRes), sKOrig // 2)
    
    if mRange.ndim == 1:
        # Case: single frequency range
        
        # Find nearest neighbor index for range StopStart
        sRadBegStp = np.abs(vOmegaVw - mRange[0])
        sNNidxRadBegStp = np.argmin(sRadBegStp)
        
        # Find nearest neighbor index for range start
        sRadBeg = np.abs(vOmegaVw - mRange[1])
        sNNidxRadBeg = np.argmin(sRadBeg)
        
        # Find nearest neighbor index for range end
        sRadEnd = np.abs(vOmegaVw - mRange[2])
        sNNidxRadEnd = np.argmin(sRadEnd)
        
        # Find nearest neighbor index for range StopEnd
        sRadEndStp = np.abs(vOmegaVw - mRange[3])
        sNNidxRadEndStp = np.argmin(sRadEndStp)
        
        # Select the section of vW and corresponding frequency axis
        vWsec = vW[sNNidxRadBegStp:sNNidxRadEndStp]
        vOmegaVwSec = vOmegaVw[sNNidxRadBegStp:sNNidxRadEndStp]
        
        # Create scaled frequency axis for interpolation (vK points)
        vWscaled = np.linspace(vOmegaVw[sNNidxRadBegStp], vOmegaVw[sNNidxRadEndStp-1], vK)
        
        # Linear interpolation
        f_interp = interp1d(vOmegaVwSec, vWsec, kind='linear')
        vW_D = f_interp(vWscaled)
        
        return vW_D.T
    
    elif mRange.ndim > 1:
        # Case: multiple frequency ranges
        for mm in range(mRange.ndim):
            # Find nearest neighbor indices for current range
            sRadBegStp = np.abs(vOmegaVw - mRange[mm,0])
            sNNidxRadBegStp = np.argmin(sRadBegStp)
            
            # Find nearest neighbor index for range start
            sRadBeg = np.abs(vOmegaVw - mRange[mm,1])
            sNNidxRadBeg = np.argmin(sRadBeg)
            
            # Find nearest neighbor index for range end
            sRadEnd = np.abs(vOmegaVw - mRange[mm,2])
            sNNidxRadEnd = np.argmin(sRadEnd)
            
            # Find nearest neighbor index for range StopEnd
            sRadEndStp = np.abs(vOmegaVw - mRange[mm,3])
            sNNidxRadEndStp = np.argmin(sRadEndStp)
            
            # Select the section of vW and corresponding frequency axis
            vWsec = vW[sNNidxRadBegStp:sNNidxRadEndStp]
            vOmegaVwSec = vOmegaVw[sNNidxRadBegStp:sNNidxRadEndStp]
            
            # Create scaled frequency axis for interpolation (vK points)
            vWscaled = np.linspace(vOmegaVw[sNNidxRadBegStp], vOmegaVw[sNNidxRadEndStp-1], vK)
            
            # Linear interpolation
            f_interp = interp1d(vOmegaVwSec, vWsec, kind='linear')
            vW_D = f_interp(vWscaled)
            
            # Stack interpolated sections
            if mm == 0:
                vW_D = vW_range.copy()
            else:
                vW_D = np.hstack([vW_D, vW_range])
        
        return vW_D
    
    else:
        # If mRange is not specified or has incorrect dimensions
        return None