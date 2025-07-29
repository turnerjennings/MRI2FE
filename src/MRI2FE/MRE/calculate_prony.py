from math import sqrt
from scipy.optimize import minimize
import numpy as np

from typing import Union


def _is_array(obj):
    """Function to check if an object is array-like"""

    check_type = isinstance(obj, (list, tuple, np.ndarray))

    return check_type


def calculate_prony(
    gp: Union[float, list, tuple, np.ndarray] = None,
    gpp: Union[float, list, tuple, np.ndarray] = None,
    mu: Union[float, list, tuple, np.ndarray] = None,
    xi: Union[float, list, tuple, np.ndarray] = None,
    w: Union[float, list, tuple, np.ndarray] = None,
    tol=1.00e-3,
):
    """Calculate the 1st order prony series equivalent to the complex shear modulus representation of brain tissue viscoelasticity

    Args:
        gp (array-like): storage moduli at different frequencies. If gp is provided, gpp must also be provided.
        gpp (array-like): loss moduli at different frequencies.  If gpp is provided, gp must also be provided.
        mu (array-like): Shear stiffness at different frequencies.  If mu is provided, xi must also be provided.
        xi (array-like): Damping ratio at different frequencies.  If xi, is provided, mu must also be provided.
        w (array-like): MRE frequency for each gp/gpp or mu/xi value, must follow the same order.
        tol (float, optional): Back calculation tolerance check. Defaults to 1.00e-3.

    Raises:
        ValueError: gp/gpp or mu/xi not provided
        ValueError: number of inputs on each array do not match
        ValueError: Mu or Xi back-calculation out of range

    Returns:
        Ginf (float): long-time shear modulus
        G1 (float): short-time shear modulus
        tau (float): short-time time constant
        gp (float): Storage modulus
        gpp (float): Loss modulus
    """

    # check for valid input
    first_set_valid = all(param is not None for param in [gp, gpp, w])
    second_set_valid = all(param is not None for param in [mu, xi, w])

    # Raise error if neither parameter set is complete
    if not (first_set_valid or second_set_valid):
        raise ValueError(
            "You must provide either (gp, gpp, w) or (mu, xi, w) as inputs"
        )

    # check if all are array like
    first_set_type = all(_is_array(o) for o in [gp, gpp, w])
    second_set_type = all(_is_array(o) for o in [mu, xi, w])

    if not (first_set_type or second_set_type):
        raise ValueError("Unmatched number of values for input parameters")

    # convert to numpy array
    if first_set_type:
        gp = np.asarray(gp)
        gpp = np.asarray(gpp)
        w = np.asarray(w)

    if second_set_type:
        mu = np.asarray(mu)
        xi = np.asarray(xi)
        w = np.asarray(w)

    # calculate complex modulus from stiffness and damping ratio (single)
    if second_set_valid and not second_set_type:
        mu = abs(mu)
        xi = abs(xi)

        a = sqrt(1.0 + 4.0 * xi**2)

        gp = ((1.0 + a) * mu) / (2.0 * a**2)
        gpp = 2.0 * xi * gp
        gmag = sqrt(gp**2 + gpp**2)

        # back-calculate mu and xi to check calculation, raise error if they don't match
        mu_bc = 2.0 * gmag**2 / (gp + gmag)
        xi_bc = gpp / (2.0 * gp)

        if (mu - mu_bc) ** 2 >= tol:
            raise ValueError(
                f"Error in mu back-calculation: mu = {mu}, mu_bc = {mu_bc}."
            )
        if (xi - xi_bc) ** 2 >= tol:
            raise ValueError(
                f"Error in xi back-calculation, xi = {xi}, xi_bc = {xi_bc}.."
            )

    # calculate complex modulus from stiffness and damping ratio (array)
    elif second_set_valid:
        mu = np.abs(mu)
        xi = np.abs(xi)

        a = np.sqrt(1.0 + 4.0 * np.square(xi))

        gp = ((1.0 + a) * mu) / (2.0 * np.square(a))
        gpp = 2.0 * xi * gp
        gmag = np.sqrt(np.square(gp) + np.square(gpp))

        # back-calculate mu and xi to check calculation, raise error if they don't match
        mu_bc = 2.0 * np.square(gmag) / (gp + gmag)
        xi_bc = gpp / (2.0 * gp)

        if np.sum((mu - mu_bc) ** 2) >= tol:
            raise ValueError(
                f"Error in mu back-calculation: mu = {mu}, mu_bc = {mu_bc}."
            )
        if np.sum((xi - xi_bc) ** 2) >= tol:
            raise ValueError(
                f"Error in xi back-calculation, xi = {xi}, xi_bc = {xi_bc}.."
            )

    # calculate te prony series constants

    Ginf, G1, tau = prony_series(gp, gpp, w)

    return Ginf, G1, tau


def prony_series(
    gp: Union[float, np.ndarray],
    gpp: Union[float, np.ndarray],
    w: Union[float, np.ndarray],
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
        x0 = np.array([1.0, 1.0, 1.0 / w])
    else:
        x0 = np.array([1.0, 1.0, 1.0 / np.mean(w)])

    # Bounds for variables (all non-negative)
    bounds = [(0, None), (0, None), (0, None)]

    # Objective function to minimize
    def objective(x):
        g1, g2, tau = x
        # Calculate the expressions for gp and gpp
        gpex = g1 + (g2 * w**2 * tau**2) / (1 + w**2 * tau**2)
        gppex = (g2 * w * tau) / (1 + w**2 * tau**2)

        # Return the sum of squared errors
        return np.sum((gp - gpex) ** 2 + (gpp - gppex) ** 2)

    # Constraint: gppex/gpex = gpp/gp
    def constraint(x):
        g1, g2, tau = x
        gpex = g1 + (g2 * w**2 * tau**2) / (1 + w**2 * tau**2)
        gppex = (g2 * w * tau) / (1 + w**2 * tau**2)

        if gpex.any() == 0.0:
            return np.mean(gppex) - np.mean(gpp / gp)
        else:
            return np.mean(gppex / gpex) - np.mean(gpp / gp)

    # Define the constraint as a dictionary
    cons = {"type": "eq", "fun": constraint}

    # Solve the optimization problem
    result = minimize(
        objective,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=cons,
        options={"ftol": 1e-10},
    )

    # Extract the optimized values
    g1, g2, tau = result.x

    return g1, g2, tau
