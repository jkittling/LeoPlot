import h5py
import glob
import numpy as np
import pandas as pd

# Particle Types: 
# 1 - DM, 2 - Disk, 3 - Bulge, 5 - SM Black Hole
PartTypesIDs = ['PartType1/ParticleIDs','PartType2/ParticleIDs','PartType3/ParticleIDs','PartType5/ParticleIDs']
PartTypesCs  = ['PartType1/Coordinates','PartType2/Coordinates','PartType3/Coordinates','PartType5/Coordinates']

# Mass per star particle
dm_star = 1e4
# Mass per dark matter particle
dm_dark = 1e5

# Particles after this should belong to Leo I
cuts = [1070000, 344000, 43000]

def split(run, PartType, j):
    """
    Returns the particle positions of MW and Leo at the jth snapshot
    based on the IC file. PartType is an integer referring to the particle
    type, which could be 1 for Dark Matter, 2 for Disk and 3 for Bulge.
    """
    #Read files
    IC = h5py.File(run+"ICs-grid.hdf5", "r")
    snaps = sorted(glob.glob(run+"output/snap_*.hdf5"))
    CC = h5py.File(snaps[j], "r")

    i = PartType - 1

    #Get IDs of Leo I via initial condition
    ids_leo = IC[PartTypesIDs[i]][cuts[i]:]
    #Also get all the current ids
    ids_all = CC[PartTypesIDs[i]][:]

    #Get indices of Leo I IDs in the current snapshot
    mask = np.isin(ids_all, ids_leo)
    #mask_mw = np.invert(mask)

    #Apply masks to current positions
    poss_leo = CC[PartTypesCs[i]][:][mask]
    #poss_mw  = CC[PartTypesCs[i]][:][mask_mw]
    return poss_leo

def M_rad(poss, cen_pos, R, PartType):
    """
    Returns the solar mass contained in a given radius R from Leo I. Position
    should be in Cartesian units. Radius should be in same units (kpc).
    PartType can be "star" or "dark" depending if the particles are stellar
    or dark matter (2 or 3 vs. 1).
    """
    #Get separation of particles from hypothesized center of Leo I
    poss_leo = poss - cen_pos[0]
    dist_leo = np.sqrt(poss_leo[:,0]**2 + poss_leo[:,1]**2 + poss_leo[:,2]**2)

    #Get those indices of particles within radius R
    counter = dist_leo < R

    #"count" how many particles within virial radius
    number  = counter.sum()

    #Multiply number of particles by mass per particle to get mass
    if PartType != 1:
        mass = number * dm_star
    else:
        mass = number * dm_dark
    return mass

def M_array(run, R, nsnaps):
    """
    Returns the mass enclosed within a radius R for every snapshot.
    Goes up to the number of finished snapshots, nsnaps, which
    should be 80 for a complete run.
    """
    Ms = np.zeros(nsnaps)

    for j in range(nsnaps):
        cen_pos    = split(run, 5, j)[0]
        star_poss1 = split(run, 2, j)
        star_poss2 = split(run, 3, j)
        star_poss = np.append(star_poss1, star_poss2, axis = 0)
        dark_poss = split(run, 1, j)

        M = M_rad(star_poss,cen_pos,R, 2) + M_rad(dark_poss,cen_pos,R, 1)
        Ms[j] += M
    return Ms

