import numpy as np
import gurobipy as gp
from gurobipy import GRB
from scipy.signal import get_window

import sys
sys.path.append('../../01_Library')
# individual packages
import sg, sa, sp, obq, filt


def OptDFT(vEl_hat, mRIF_m, vx_m, vb_mInit):
    """
    Args:
        vx: Input vector.
        vW: Desired spectral 
        K: K values for DFT
        
    Returns:
        vb: Quantized one-bit vector
        ve: Error vector
    """
    
    sM = len(vx_m)

    if np.all(vb_mInit == 0):
        sMean = np.mean(vx_m)
        bInit = np.full_like(vx_m, 0)
        bInit[vx_m >= sMean] = 1
    else:
        bInit = 2*vx_m.copy()+1
        
    vb_out  = np.zeros((sM,1)).flatten()
    # GUROBI
    #Mixed-Integer Quadratically Constrained Quadratic Programming (MIQP)
    
    model = gp.Model("MIQP")
    model.setParam("OutputFlag", 0)     # 0 to Suppress Gurobi output
    model.setParam("MIPFocus", 3)               # Focus on proving optimality
    model.setParam("Method", 2)                 # Barrier method for root relaxation (robust for QCQP)
    model.setParam("Presolve", 2)               # Aggressive presolve
    model.setParam("Heuristics", 0.1)          # Minimal heuristic distraction
    model.setParam("Cuts", 3)                   # Most aggressive cutting planes
    model.setParam("VarBranch", 3)              # Strong branching
    model.setParam("TimeLimit",1)
    
    #model.setParam("TimeLimit",2)
    #model.setParam("VarBranch", 3)
    #model.setParam("MIPFocus", 3)       # Shift focus to finding good feasible solutions quickly
    #model.setParam("Heuristics", 0.9)   # Increase heuristic efforts
    #model.setParam("Presolve", 2)       # More aggressive presolve
    #model.setParam('Method', -1)
    
    # Decision variables (vb) as binary, mapped to {-1, 1} in the objective
    vb = model.addVars(sM, vtype=gp.GRB.BINARY, name="vb")
    
    for j in range(sM):
         vb[j].Start = bInit[j].copy()
    
    model.update()
    
    vbDec = {j: 2 * vb[j] - 1 for j in range(sM)}
    vdHat_m = {j: vx_m[j] - vbDec[j] for j in range(sM)}

    # --- Term 1: -2 * X̃_L^T * (R * F_wM * b)
    obj = gp.QuadExpr()

    for i in range(mRIF_m.shape[0]):  # Rows of mW
        se = vEl_hat[i].copy()
        for j in range(sM):  # Elements of vb and vx
            # Contribution of each element to the quadratic term
            se = se + mRIF_m[i, j] * vdHat_m[j] 
        obj += se * se 

    # Set the objective
    model.setObjective(obj, GRB.MINIMIZE)
    model.optimize()
    
    for i in range(sM):
            vb_out[i] = (2 * vb[i].X - 1)      
    
    sBlockErr = vEl_hat + mRIF_m @ (vx_m - vb_out)
    sBlockErr = sBlockErr.T @ sBlockErr   
  
    # Output the solution
    if model.status == GRB.OPTIMAL:
        outTxt = f"Optimal solution found. ({model.SolCount})"
    else:
        outTxt = f"No Optimal solution found. ({model.SolCount})"
        
    return vb_out, sBlockErr, outTxt