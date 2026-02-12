
import sys, os, time
from subprocess import call
from mpi4py import MPI


MPIcomm = MPI.COMM_WORLD
pId = MPIcomm.Get_rank()
nProc = MPIcomm.Get_size()
currentDirectory = os.getcwd()

INBASE = sys.argv[1]

OUTBASE = f"{INBASE}/Rockstar" #where to output files

rockstarComand ='/home/r408l055/rockstar/rockstar'

rockstarConf = {
'FILE_FORMAT': '"AREPO"',
'AREPO_LENGTH_CONVERSION' :1e-3,  #convert from kpc to Mpc
'AREPO_MASS_CONVERSION': 1e+10,
'FORCE_RES': 100000 / 2**11 / 39,                 #Mpc/h set up for 811 right now
'OUTBASE': OUTBASE,
}

parallelConf = {
'PARALLEL_IO': 1,
'INBASE':  INBASE ,               #"/directory/where/files/are/located"
'NUM_BLOCKS': 1,                              # <number of files per snapshot>
'NUM_SNAPS': 128,                               # <total number of snapshots>
'STARTING_SNAP': 0,
'FILENAME': '"snap_<snap>.hdf5"',              #"my_sim.<snap>.<block>" need to include file extension
'NUM_WRITERS': 16,                             #<number of CPUs>
'FORK_READERS_FROM_WRITERS': 1,
'FORK_PROCESSORS_PER_MACHINE': 16,             #<number of processors per node>
}


if pId == 0:
  if not os.path.exists( rockstarConf['OUTBASE']): os.makedirs(rockstarConf['OUTBASE'])
  rockstarconfigFile = rockstarConf['OUTBASE'] + '/rockstar_param.cfg'
  rckFile = open( rockstarconfigFile, "w" )
  for key in list(rockstarConf.keys()):
    rckFile.write( key + " = " + str(rockstarConf[key]) + "\n" )
  for key in list(parallelConf.keys()):
    rckFile.write( key + " = " + str(parallelConf[key]) + "\n" )
  rckFile.close()
  #Run ROCKSTAR finder
  print("\nFinding halos...")
  print(" Parallel configuration")
  print("Output: ", rockstarConf['OUTBASE'] + '\n')


MPIcomm.Barrier()

if pId == 0: call([rockstarComand, "-c", rockstarconfigFile ])
if pId == 1:
  time.sleep(5)
  call([rockstarComand, "-c", rockstarConf['OUTBASE'] + '/auto-rockstar.cfg' ])
