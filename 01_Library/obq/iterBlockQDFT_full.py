import numpy as np
import obq  # deine Optimierungsfunktion
import scipy
from scipy import signal
import sg, sa, globalTools, filt  # Assumes sg.dftMat(N, K) returns (K x N) matrix
import matplotlib.pyplot as plt

def iterBlockQDFT(vx, mRange, vW_D, sK, sM, sL_max, sHop, sType='grb', verbose=True, nIter=1):
    """
    Iterative block-based one-bit DFT-domain quantization with outer iteration.
    """

    #vx = np.concatenate([np.zeros(sL), vx])  # Padding
    sxLen = len(vx)
    sNumBlocks = (sxLen - sM) // sHop + 1
    
    vb = np.zeros(sxLen)
    
    if mRange.ndim == 1:
        F_D,_ = sg.dftMat(sxLen, sK, mRange) 
    elif mRange.ndim > 1:
        for mm in range(mRange.ndim):
            FD_m,_ = sg.dftMat(sxLen, sK, mRange[mm,:]) 
            if mm == 0:
                F_D = FD_m.copy()
            elif mm > 0:
                F_D = np.vstack([F_D, FD_m])
    
    F_D = np.diag(vW_D) @ F_D
    mRIF_D          = np.vstack([F_D.real, F_D.imag])
    vXRI            = mRIF_D @ vx
    vGlobalBlockErr = np.zeros(sNumBlocks)
    mErr_k_p        = np.zeros((mRIF_D.shape[0], sNumBlocks))
    
    for iIter in range(nIter):
        if verbose:
            print(f"=== Iteration {iIter+1}/{nIter} ===")
        sM += 16*iIter
        for m in range(sNumBlocks):
            if m == 0:
                progressDFTBlock = globalTools.SimpleProgressBar(
                    sNumBlocks, width=40, prefix = "DFT-BlockOpt", fill="█", empty=" ", end=" ✓")
            
            sStIdx = m * sHop
            sEndIdx = sStIdx + sM

            if sEndIdx > sxLen:
                if verbose:
                    print(f"Skipping block {m}: exceeds signal length.")
                continue
            
            if m==0:
                vEl_hat = np.zeros(sK*2)
            
            #### Calculate blocks and projection matrices for the quantization                
            mRIF_m  = mRIF_D[:,sStIdx:sEndIdx].copy()
            vx_m    = vx[sStIdx:sEndIdx].copy()
            vx_m[0:sHop] = vx_m[0:sHop].copy()
            vb_mInit = vb[sStIdx : sEndIdx].copy()            

            vb_m, _, txtOut = obq.OptDFT(vEl_hat, mRIF_m, vx_m, vb_mInit)            
            vb[sStIdx: sEndIdx] = vb_m
            
            #### Calculating indeces for locked elements
            sErrIdx = (m+1) * sHop
            sL = (m+1) * sHop
            if sL > sL_max:
                sL = sL_max
                
            mRIF_l   = mRIF_D[:,(sErrIdx-sL):sErrIdx].copy()
            vx_l     = vx[(sErrIdx-sL):sErrIdx].copy()
            vb_l     = vb[(sErrIdx-sL):sErrIdx].copy()           
            vEl_hat  = mRIF_l @ (vx_l-vb_l)
            
            mErr_k_p[:,m] = vEl_hat.copy()
            
            #### Comparing with Global Goal
            vGlobalError       = vXRI - mRIF_D @ vb
            vGlobalBlockErr[m] = vGlobalError.T @ vGlobalError
            
            txtOut += f" BlockErr: {vGlobalBlockErr[m]:.4e}"
            progressDFTBlock.update(m+1, txtOut)

    return vb, vGlobalBlockErr, mErr_k_p