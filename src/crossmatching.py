import pandas as pd 
import numpy as np 
import importlib
from astropy.coordinates import SkyCoord
from astropy import units as u
from sklearn.neighbors import BallTree

def search_neighbors(group_KDE, group_iterator, max_sep_arcsec=60, redshift_separation=0.01, name1='group_iterator_index', name2='group_KDE_index', name3='distances'):

    def convert_to_cartesian(ra, dec):
        
        x = np.cos(np.radians(dec)) * np.cos(np.radians(ra))
        y = np.cos(np.radians(dec)) * np.sin(np.radians(ra))
        z = np.sin(np.radians(dec))
        return np.vstack((x, y, z)).T

    max_sep_rad = np.deg2rad(max_sep_arcsec / 3600)

    group_KDE_cartesian = convert_to_cartesian(group_KDE['RA'].values, group_KDE['DEC'].values)
    group_iterator_cartesian = convert_to_cartesian(group_iterator['RA'].values, group_iterator['DEC'].values)

    tree = BallTree(group_KDE_cartesian, metric='euclidean')

    matches = {name1: [], name2: [], name3: []}
    for index, pos in enumerate(group_iterator_cartesian):
        

        pos = pos.reshape(1, -1)

        indices, distances = tree.query_radius(pos, r=2 * np.sin(max_sep_rad / 2), return_distance=True)

        indices = indices[0]
        distances = distances[0]

        filtered_indices = []
        filtered_distances = []

        for i, d in zip(indices, distances):
            if not np.array_equal(group_KDE_cartesian[i], pos[0]):
                filtered_indices.append(i)
                filtered_distances.append(d)

        filtered_indices = np.array(filtered_indices)
        filtered_distances = np.array(filtered_distances)

        not_z = []

        for i in range(len(filtered_indices)):
            if group_iterator.iloc[index]['zHD'] - group_KDE.iloc[filtered_indices[i]]['zHD'] > redshift_separation:
                 not_z.append(i)

        filtered_indices = np.delete(filtered_indices, not_z)
        filtered_distances = np.delete(filtered_distances, not_z)

        index = group_iterator[name1].iloc[index]

        if len(filtered_indices) != 0:
                    matches[name1].append(index) 
                    matches[name2].append(filtered_indices[np.where(filtered_distances == min(filtered_distances))[0]])
                    matches[name3].append(min(filtered_distances))       

    return matches