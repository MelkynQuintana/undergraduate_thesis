'''
--------------------------------------------------------
------------ COSMOLOGICAL related functions ------------
--------------------------------------------------------
'''

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.constants import c 

'''
    COSMOLOGICAL MODELS
    - E_* -> E Function of * model
    - dist_* -> Luminosity Distance Function of * model
'''

class CosmologicalModel:
    default_h0 = 73.6  #Hubble Constant for SH0ES team
    default_om_m = 0.334 #Matter density for SH0ES team

class FlatModel(CosmologicalModel):

    # There are no additional parameters for the Flat model

    @classmethod
    def E_inv(cls, z, om_m):
        return 1 / np.sqrt(om_m * (1 + z)**3 + (1 - om_m))

    @classmethod
    def dist(cls, z, h0=None, om_m=None):

        h0 = h0 if h0 is not None else cls.default_h0
        om_m = om_m if om_m is not None else cls.default_om_m

        int_flat = np.array([
            quad(cls.E_inv, 0, zi, args=(om_m,))[0] for zi in z
        ])
        com_flat = ((c / 1000) / h0) * int_flat
        return (1 + z) * com_flat

    @classmethod
    def hubble_const(cls, z=None, mu=None, om_m=None):

        om_m = om_m if om_m is not None else cls.default_om_m

        int_flat = np.array([
            quad(cls.E_inv, 0, zi, args=(om_m,))[0] for zi in z
        ])
        hubble_const = ((c / 1000) * (1 + z) / 10**((mu - 25) / 5)) * int_flat

        return hubble_const

    @classmethod
    def E(cls, om_m=None, z=None):

        om_m = om_m if om_m is not None else cls.default_om_m

        return np.sqrt(om_m * (1 + z)**3 + (1 - om_m))


        
class LCDMModel(CosmologicalModel):

    om_A = 0.666 #Dark Energy density fot SH0ES team

    @classmethod
    def dist(cls, z, h0=None, om_m=None, om_A=None):

        h0 = h0 if h0 is not None else cls.default_h0
        om_m = om_m if om_m is not None else cls.default_om_m
        om_A = om_A if om_A is not None else cls.default_om_A

        # Luminosity distance for the ΛCDM model
        int_lcdm = cumulative_trapezoid(1 / np.sqrt(om_m * (1 + z)**3 + om_A * (1 + z)**2 + (1 - om_m - om_A)), z, initial=0)
        com_lcdm = ((3e5 / 1000) / h0) * int_lcdm
        return (1 + z) * com_lcdm
    
    @classmethod
    def _E_inv(cls, z, om_m):
        return 1 / np.sqrt(om_m * (1 + z)**3 + (1 - om_m))



