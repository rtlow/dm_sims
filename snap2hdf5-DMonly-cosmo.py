from readgadget import *
import h5py
import numpy as np
import os,sys
from datetime import datetime

STANDARD_FORMAT = True
NPROCS = 4
INTERP = 0

pyGRopts = {'suppress':1}

snap = sys.argv[1]
outf = sys.argv[2]

if not os.path.isfile(snap):
    print ('cannot find %s') % snap
    sys.exit()
if '.hdf5' not in outf:
    outf = '%s.hdf5' % outf

TIPSY = False
if 'bin' in snap:
    TIPSY = True
    snap  = snap[:-4]

    h   = 0.67
    Ol  = 0.69
    O0  = 0.31
    BOX = 50000.  ## in kpc/h
    #BOX = 50.0 ## in Mpc/h  (MUSIC IC)

    PROTONMASS = 1.6726e-24
    H_MASSFRAC = 0.76
    BOLTZMANN  = 1.3806e-16
    GAMMA      = 5.0/3.0

    CM_PER_MPC               = 3.085678e24
    UnitLength_in_cm         = 3.085678e21 / h
    UnitMass_in_g            = 1.989e43    / h
    UnitVelocity_in_cm_per_s = 1.e5
    L = BOX / h * 1.0e-3
    UnitTime_in_s      = UnitLength_in_cm / UnitVelocity_in_cm_per_s
    UnitDensity_in_cgs = UnitMass_in_g / UnitLength_in_cm**3
    UnitEnergy_in_cgs  = UnitMass_in_g * UnitLength_in_cm**2 / UnitTime_in_s**2
    unit_Time     = np.sqrt(8. * np.pi / 3.) * CM_PER_MPC / (100. * h * UnitVelocity_in_cm_per_s)
    unit_Density  = 1.8791e-29 * h**2
    unit_Length   = L*CM_PER_MPC
    unit_Mass     = unit_Density * unit_Length**3
    unit_Velocity = unit_Length / unit_Time
    Mconvert   = unit_Mass / UnitMass_in_g
    Lconvert   = unit_Length / UnitLength_in_cm
    Vconvert   = unit_Velocity / UnitVelocity_in_cm_per_s
    RHOconvert = (UnitLength_in_cm**3 * unit_Density) / UnitMass_in_g
    Uconvert   = 1 #UnitEnergy_in_cgs
    Ladd = 0.5

if not TIPSY:
    Mconvert = 1
    Lconvert = 1
    Ladd     = 0
    Vconvert = 1
    RHOconvert = 1
    Uconvert = 1
    
    h   = readhead(snap,'h')
    O0  = readhead(snap,'O0')
    Ol  = readhead(snap,'Ol')
    BOX = readhead(snap,'boxsize')

npart    = readhead(snap,'npart')
masses   = readhead(snap,'mass')
redshift = readhead(snap,'redshift')
time     = readhead(snap,'time')
    
OUT = h5py.File(outf,'w')
header = OUT.create_group("Header")
header.attrs['BoxSize']                = BOX
header.attrs['Flag_Cooling']           = readhead(snap,'f_cooling')
header.attrs['Flag_DoublePrecision']   = 0
header.attrs['Flag_Feedback']          = readhead(snap,'f_fb')
header.attrs['Flag_Metals']            = readhead(snap,'f_metals')
header.attrs['Flag_Sfr']               = readhead(snap,'f_sfr')
header.attrs['Flag_StellarAge']        = readhead(snap,'f_age')
header.attrs['HubbleParam']            = h
header.attrs['MassTable']              = masses
header.attrs['NumFilesPerSnapshot']    = 1
header.attrs['NumPart_ThisFile']       = npart
header.attrs['NumPart_Total']          = npart
header.attrs['NumPart_Total_HighWord'] = np.zeros(6,dtype=np.uint32)
header.attrs['Omega0']                 = O0
header.attrs['OmegaLambda']            = Ol
header.attrs['Redshift']               = redshift
header.attrs['Time']                   = time

def initGroup(i):
    if npart[i]:
        p = OUT.create_group('PartType%d' % i)
        POS = readsnap(snap,'pos',i, **pyGRopts) + Ladd        
        p.create_dataset('Coordinates', data=POS * Lconvert)
        p.create_dataset('Velocities',  data=readsnap(snap,'vel',i, **pyGRopts)*Vconvert)
        p.create_dataset('ParticleIDs', data=readsnap(snap,'pid',i, **pyGRopts))
        p.create_dataset('Masses',  data=readsnap(snap,'mass',i, **pyGRopts)*Mconvert)
        return p

## non gas & star particles
for i in [0,1,2,3,4,5]:
    p = initGroup(i)



print ('FINISHED')

