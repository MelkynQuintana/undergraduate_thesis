import pandas as pd
import numpy as np
from scipy.constants import c 

'''This funcition extract a dataframe with pantheon data and his covariance matrix'''

def Pantheon_data(cut=1):

    full_data = pd.read_csv('../data/Pantheon+/Pantheon_data.dat',
                            sep=' ', header=0)

    #Adding data from covariance file
    covariance = pd.read_csv('../data/Pantheon+/covariance.cov', skiprows=1, header=None)
    cov = np.array(covariance).reshape(len(full_data), len(full_data))

    full_data.loc[:, 'MU_err'] = np.diag(cov)


    #Adding necessary columns
    full_data.loc[:, 'Distance'] = 10**((full_data['MU_SH0ES'] - 25)/5)
    full_data.loc[:, 'Velocity'] = full_data['zHD'] * c * 1e-3
    full_data.loc[:, 'Dist_err'] = np.sqrt(((np.log(10) * 10**((full_data['MU_SH0ES'] - 25)/5))/5)**2 * (full_data['MU_err'])**2)

    #Adding relative errors
    full_data['z_relerr'] = full_data['zHDERR']/full_data['zHD']
    full_data['MU_relerr'] = full_data['MU_err']/full_data['MU_SH0ES']

    if cut==1:
        temp = ['CID','IDSURVEY', 'zHD', 'zHDERR', 'z_relerr', 'MU_SH0ES', 
                'MU_err', 'MU_relerr', 'Distance', 'Velocity', 'Dist_err',
                'IS_CALIBRATOR', 'CEPH_DIST']
        
        return full_data[temp], cov
    
    elif cut==2:
        temp = ['CID', 'IDSURVEY', 'zHD', 'x1', 'x1ERR', 'c', 'cERR', 'mB', 'mBERR', 'x0', 'x0ERR', 'MU_SH0ES']

        return full_data[temp], cov

    else:
        return full_data, cov

def DES_data(cut=1):

    full_data = pd.read_csv('../data/DES/DES_data.csv', sep=',', header=0)

    full_data['MU_SH0ES'] = full_data['MU']
    full_data['DEC'] = full_data['HOST_DEC']
    full_data['RA'] = full_data['HOST_RA']

    #degeneration_factor = 5 * np.log10(c / 1000 * 70)

    #full_data['MU_SH0ES'] = full_data['MU'] - degeneration_factor 

    #Adding data from covariance file
    covariance = pd.read_csv('../data/DES/cov_matrix.txt', skiprows=1, header=None)
    cov = np.array(covariance).reshape(len(full_data), len(full_data))

    full_data.loc[:, 'MU_err'] = np.diag(cov)

    #Adding necessary columns
    full_data.loc[:, 'Distance'] = 10**((full_data['MU_SH0ES'] - 25)/5)
    full_data.loc[:, 'Velocity'] = full_data['zHD'] * c * 1e-3
    full_data.loc[:, 'Dist_err'] = np.sqrt(((np.log(10) * 10**((full_data['MU_SH0ES'] - 25)/5))/5)**2 * (full_data['MU_err'])**2)

    #Adding relative errors
    full_data['z_relerr'] = full_data['zHDERR']/full_data['zHD']
    full_data['MU_relerr'] = full_data['MU_err']/full_data['MU_SH0ES']


    if cut==1:
        temp = ['CID','IDSURVEY', 'zHD', 'zHDERR', 'z_relerr', 'MU_SH0ES', 
                'MU_err', 'MU_relerr', 'Distance', 'Velocity', 'Dist_err']
        
        return full_data[temp], cov
    
    elif cut==2:
        temp = ['CID', 'IDSURVEY', 'zHD', 'x1', 'x1ERR', 'c', 'cERR', 'mB', 'mBERR', 'x0', 'x0ERR', 'MU_SH0ES']

        return full_data[temp], cov

    else:
        return full_data, cov

def cluster_data(file, base_path='../data/Clusters/'):

    filepath = base_path + file
    
    data = pd.read_csv(filepath,
                        sep='\s+', names=['ID_halo', 'Richness', 'ignore', 'RA', 'DEC', 'zHD'])

    return data

def cluster_properties(file, base_path='../data/Clusters/', cut=True):

    filepath = base_path + file
    
    data = pd.read_csv(filepath,
                        sep='\s+', names=['ID_halo', 'Richness', 'ID', 'oID', 'zHD', 'x', 'y', 'z', 'Mvir', 'Rvir', 'Lr18', 'StellarMass', 'MhStellar', 'MhL18', 'galStellarMass', 'Rzspace'])

    if cut:

        temp = ['ID_halo', 'ID', 'oID', 'zHD', 'Mvir', 'Rvir']

        return data[temp]

    else:
        return data


def galaxy_data_mass(file, base_path='../data/Clusters/', cut=True):

    filepath = base_path + file
    
    data = pd.read_csv(filepath,
                        sep='\s+', names=['ID', 'oID', 'ID_halo', 'StellarMass', 'x', 'y', 'z'])

    if cut:

        temp = ['ID', 'oID', 'ID_halo', 'StellarMass']

        return data[temp]

    else:
        return data

def galaxy_data_radec(file, base_path='../data/Clusters/', cut=True):

    filepath = base_path + file
    
    data = pd.read_csv(filepath,
                        sep='\s+', names=['oID', 'RA', 'DEC', 'Mu', 'Mg', 'Mr', 'Mi', 'Mz', 'zHD', 'Lr', 'c', 'Mr0', 'gr0'])

    if cut:

        temp = ['oID', 'RA', 'DEC', 'zHD']

        return data[temp]

    else:
        return data