#!/bin/bash            # this line only there to enable syntax highlighting in this file

NTYPES=6                       # number of particle types
PERIODIC
VORONOI
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE
TREE_BASED_TIMESTEPS     # non-local timestep criterion (take 'signal speed' into account)
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
SELFGRAVITY                   # switch on for self-gravity     
ENFORCE_JEANS_STABILITY_OF_CELLS    # this imposes an adaptive floor for the temperature
EVALPOTENTIAL                 # computes gravitational potential
NSOFTTYPES=6                  # Number of different softening values to which particle types can be mapped.
MULTIPLE_NODE_SOFTENING       # If a tree node is to be used which is softened, this is done with the softenings of its different mass components 
INDIVIDUAL_GRAVITY_SOFTENING=32  # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
ADAPTIVE_HYDRO_SOFTENING
PMGRID=512
RCUT=5.5
PLACEHIGHRESREGION=2
ENLARGEREGION=1.2
GRIDBOOST=1
PM_ZOOM_OPTIMIZED
CHUNKING                 # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
OUTPUT_COORDINATES_IN_DOUBLEPRECISION    # will always output coordinates in double precision
FOF                                # enable FoF output
FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
SUBFIND                            # enables substructure finder
VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
ENLARGE_DYNAMIC_RANGE_IN_TIME  # This extends the dynamic range of the integer timeline from 32 to 64 bit
REDUCE_FLUSH
OUTPUT_CENTER_OF_MASS
OUTPUTPOTENTIAL
HAVE_HDF5                     # needed when HDF5 I/O support is desired
DEBUG                         # enables core-dumps 
HOST_MEMORY_REPORTING         # reports after start-up the available system memory by analyzing /proc/meminfo
LONGIDS
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
SIDM=2                                #activate and set types
SIDM_STATES=2                         #number of DM states (for inelastic models)
SIDM_REACTIONS=12                      #number of scatter reactions (for inelasitc models)
FOF_SECONDARY_LINK_TYPES=4   # 2^type for the types linked to nearest primaries
