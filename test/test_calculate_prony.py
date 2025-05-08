import pytest
import numpy as np

from MRI2FE.MRE import calculate_prony


def generate_complex_modulus_data():
    # define time-domain prony series terms
    g0 = 10

    g1 = 15

    tau = 0.1

    # define measurement frequencies and complex shear moduli
    w = np.array([5, 10, 15, 20, 25])

    gp = g0 + (g1 * tau**2 * w**2) / (1 + w**2 * tau**2)

    gpp = (g1 * tau * w) / (1 + w**2 * tau**2)

    return gp, gpp, w, g0, g1, tau


class Test_prony_series:
    def test_full_data(self):
        tolerance = 1e-5

        gp, gpp, w, ref_g0, ref_g1, ref_tau = generate_complex_modulus_data()

        g0, g1, tau = calculate_prony(gp=gp, gpp=gpp, w=w)

        assert (g0 - ref_g0) ** 2 <= tolerance
        assert (g1 - ref_g1) ** 2 <= tolerance
        assert (tau - ref_tau) ** 2 <= tolerance

    def test_mu_xi(self):
        tolerance = 1e-5

        gp, gpp, w, ref_g0, ref_g1, ref_tau = generate_complex_modulus_data()

        gmag = np.sqrt(gp**2 + gpp**2)

        # back-calculate mu and xi to check calculation, raise error if they don't match
        mu = 2.0 * gmag**2 / (gp + gmag)
        xi = gpp / (2.0 * gp)

        g0, g1, tau = calculate_prony(mu=mu, xi=xi, w=w)

        assert (g0 - ref_g0) ** 2 <= tolerance
        assert (g1 - ref_g1) ** 2 <= tolerance
        assert (tau - ref_tau) ** 2 <= tolerance

    def test_raises(self):
        with pytest.raises(ValueError):
            g0, g1, tau = calculate_prony(mu=1.0, gp=1.0)
            g0, g1, tau = calculate_prony(mu=1.0, gp=1.0, w=1.0)
            g0, g1, tau = calculate_prony(mu=[1.0], gp=1.0, w=1.0)
