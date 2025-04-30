from julia import Main
from math import sqrt


def calculate_prony(mu: float, xi: float, w: float, tol=1.00e-3):
    """Calculate the 1st order prony series equivalent to the complex shear modulus representation of brain tissue viscoelasticity

    Args:
        mu (float): Shear stiffness [MPa]
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

    # calculate complex modulus from stiffness and damping ratio
    a = sqrt(1 + 4 * xi**2)

    gp = ((1 + a) * mu) / (2 * a**2)
    gpp = 2 * xi * gp
    gmag = sqrt(gp**2 + gpp**2)

    # back-calculate mu and xi to check calculation, raise error if they don't match

    mu_bc = 2 * gmag**2 / (gp + gmag)
    xi_bc = gpp / (2 * gp)

    if (mu_bc - mu) ** 2 >= tol:
        raise ValueError("Mu back-calculation out of tolerance range")

    if (xi_bc - xi) ** 2 >= tol:
        raise ValueError("Xi back-calculation out of tolerance range")

    # add the Julia optimization script

    Main.include("shear_lsr.jl")

    OptModel = Main.OptModel

    # calculate te prony series constants using the julia script

    Ginf, G1, tau = OptModel.optimize_model(gp, gpp, float(w))

    return Ginf, G1, tau, gp, gpp
