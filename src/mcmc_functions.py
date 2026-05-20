'''
--------------------------------------------------------
---------------- MCMC related functions ----------------
--------------------------------------------------------
'''

import numpy as np
import emcee
import tqdm
from multiprocessing import Pool, cpu_count
from scipy.integrate import cumulative_trapezoid, quad, quad_vec
from scipy.constants import c
import my_module as mm

def E_flat(z, om_m):
    return 1 / np.sqrt(om_m * (1 + z)**3 + (1 - om_m))


# Luminosity distance function for MCMC
def dist(theta, model, data):
    # Parameters
    om_m, h0 = theta

    # Integral 
    int_flat = np.array([quad(model.E_inv, 0, z, args=(om_m,))[0] for z in data['zHD']])

    # Comoving distance
    com_flat = ((c / 1000) / h0) * int_flat

    # Luminosity distance
    return (1 + data['zHD']) * com_flat

# Priors
def prior(theta):
    om_m, h0 = theta
    if 0.1 < om_m < 1.0 and 10 < h0 < 100:
        return 0.0
    return -np.inf

def ln_L(theta, data, cov, model):
    """Log likelihood, given a cosmological model"""
    om_m, h0 = theta
    N = len(data['zHD'])

    # Compute the model distance and theoretical mu
    model_dist = dist(theta, model, data)
    
    if np.any(model_dist <= 0):
        return -np.inf

    mu_theor = 5 * np.log10(abs(model_dist)) + 25

    data.loc[:, 'mu_theor'] = mu_theor

    residuals = np.array(data['MU_SH0ES'] - data['mu_theor'])

    # Chi-squared calculation
    chi2 = residuals.T @ cov @ residuals
    
    return -0.5 * chi2

# Log probability function
def lnprob(theta, data, cov, model):
    lp = prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + ln_L(theta, data, cov, model)



'''
    MCMC FUNCTIONS
'''

base_initials = np.array([0.3, 71]) # Default value for h0

def Markov(data, cov, model, initials=base_initials, burnin=1000):

    #Set up some walkers in a gaussian ball around the maximum likelihood result
    num_params = 2
    num_walkers = 1000
    niter = 5000

    p0 = [np.array(initials) + 1e-7 * np.random.randn(num_params) for i in range(num_walkers)]

    with Pool(processes=4) as pool:   

        #Set up the sampler
        sampler = emcee.EnsembleSampler(num_walkers, num_params, lnprob, args=(data, cov, model), pool=pool)

        print("Running burn-in...")
        p0, _, _ = sampler.run_mcmc(p0, burnin, progress=True)
        sampler.reset()

        print("Running production...")
        pos, prob, state = sampler.run_mcmc(p0, niter, progress=True)

    return sampler, pos, prob, state