from julia import Main
from math import sqrt
from scipy.optimize import minimize
import numpy as np

from typing import Union

def _is_array(obj):
    """Function to check if an object is array-like
    """

    check_type = isinstance(obj, (list,tuple,np.ndarray))
    
    return check_type

def calculate_prony(
        gp: Union[float,list,tuple,np.ndarray] = None, 
        gpp: Union[float,list,tuple,np.ndarray] = None, 
        mu: Union[float,list,tuple,np.ndarray] = None, 
        xi: Union[float,list,tuple,np.ndarray] = None, 
        w: Union[float,list,tuple,np.ndarray] = None, 
        tol=1.00e-3
        ):
    """Calculate the 1st order prony series equivalent to the complex shear modulus representation of brain tissue viscoelasticity

    Args:
        mu (float or array): Shear stiffness [MPa]
        xi (float): Damping ratio
        w (float): vibration frequency [Hz]
        tol (float, optional): Back calculation tolerance check. Defaults to 1.00e-3.

    Raises:
        ValueError: Mu or Xi back-calculation out of range

    Returns:
        Ginf (float): long-time shear modulus
        G1 (float): short-time shear modulus
        tau (float): short-time time constant
        gp (float): Storage modulus
        gpp (float): Loss modulus
    """

    #check for valid input
    first_set_valid = all(param is not None for param in [gp, gpp, w])
    second_set_valid = all(param is not None for param in [mu, xi, w])
    
    # Raise error if neither parameter set is complete
    if not (first_set_valid or second_set_valid):
        raise ValueError("You must provide either (gp, gpp, w) or (mu, xi, w) as inputs")
    

    #check if all are array like
    first_set_type = all(_is_array(o) for o in [gp,gpp,w]) 
    second_set_type = all(_is_array(o) for o in [mu,xi,w]) 

    #convert to numpy array
    if first_set_type:
        gp = np.asarray(gp)
        gpp = np.asarray(gpp)
        w = np.asarray(w)

    if first_set_type:
        mu = np.asarray(mu)
        xi = np.asarray(xi)
        w = np.asarray(w)


    # calculate complex modulus from stiffness and damping ratio (single)
    if second_set_valid and not second_set_type:
        a = sqrt(1. + 4. * xi**2)

        gp = ((1. + a) * mu) / (2. * a**2)
        gpp = 2. * xi * gp
        gmag = sqrt(gp**2 + gpp**2)

        # back-calculate mu and xi to check calculation, raise error if they don't match
        mu_bc = 2. * gmag**2 / (gp + gmag)
        xi_bc = gpp / (2. * gp)

        if (mu - mu_bc)**2 >= tol:
            raise ValueError("Error in mu back-calculation, values do not match.")
        if (xi - xi_bc)**2 >= tol:
            raise ValueError("Error in xi back-calculation, values do not match.")

    # calculate complex modulus from stiffness and damping ratio (array)
    elif second_set_valid:

        a = np.sqrt(1. + 4. * xi**2.)

        
        gp = ((1. + a) * mu) / (2. * a**2)
        gpp = 2. * xi * gp
        gmag = np.sqrt(gp**2 + gpp**2)

        # back-calculate mu and xi to check calculation, raise error if they don't match
        mu_bc = 2. * gmag**2 / (gp + gmag)
        xi_bc = gpp / (2. * gp)

        if np.sum((mu - mu_bc)**2) >= tol:
            raise ValueError("Error in mu back-calculation, values do not match.")
        if np.sum((xi - xi_bc)**2) >= tol:
            raise ValueError("Error in xi back-calculation, values do not match.")            


    # calculate te prony series constants

    Ginf, G1, tau = prony_series(gp,gpp,w)

    return Ginf, G1, tau



def prony_series(
        gp: Union[float,np.ndarray], 
        gpp: Union[float,np.ndarray],
        w: Union[float,np.ndarray]
        ):
    """
    Calculate the 1-term prony series constants for an equivalent complex shear modulus using least squares optimization

    Arguments:
        gp (float): Storage modulus
        gpp (float): Loss modulus
        w (float): Modulus frequency

    Returns:
        Ginf (float): Long-time shear modulus
        G1 (float): short-time shear modulus
        tau (float): time constant
    """

    # Ensure all inputs are ndarrays for consistent handling
    gp = np.asarray(gp)
    gpp = np.asarray(gpp)
    w = np.asarray(w)
    
    # Check if we have scalar or vector inputs
    is_scalar_input = gp.ndim == 0 and gpp.ndim == 0 and w.ndim == 0
    
    # Initial guess for optimization variables [g1, g2, tau]
    if is_scalar_input:
        x0 = np.array([0.1, 0.2, 1.0/w])
    else:
        x0 = np.array([0.1, 0.2, 1.0/np.mean(w)])
    
    # Bounds for variables (all non-negative)
    bounds = [(0, None), (0, None), (0, None)]
    
    # Objective function to minimize
    def objective(x):
        g1, g2, tau = x
        # Calculate the expressions for gp and gpp
        gpex = g1 + (g2 * w**2 * tau**2) / (1 + w**2 * tau**2)
        gppex = (g2 * w * tau) / (1 + w**2 * tau**2)
        
        # Return the sum of squared errors
        return np.sum((gp - gpex)**2 + (gpp - gppex)**2)
    
    # Constraint: gppex/gpex = gpp/gp
    def constraint(x):
        g1, g2, tau = x
        gpex = g1 + (g2 * w**2 * tau**2) / (1 + w**2 * tau**2)
        gppex = (g2 * w * tau) / (1 + w**2 * tau**2)
        
        if gpex.any() == 0.0:
            return np.mean(gppex) - np.mean(gpp/gp)
        else:
            return np.mean(gppex/gpex) - np.mean(gpp/gp)
    
    # Define the constraint as a dictionary
    cons = {'type': 'eq', 'fun': constraint}
    
    # Solve the optimization problem
    result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=cons, options={'ftol': 1e-10})
    
    # Extract the optimized values
    g1, g2, tau = result.x
    
    return g1, g2, tau

