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
FOF_SECONDARY_LINK_TYPES=1+4+16+32   # 2^type for the types linked to nearest primaries
MHD
MHD_POWELL
MHD_POWELL_LIMIT_TIMESTEP
MHD_SEEDFIELD
COOLING
UVB_SELF_SHIELDING            # gas is self-shielded from the cosmic background based on its density
USE_SFR
RIEMANN_HLLD
REFINEMENT_HIGH_RES_GAS
ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
ADD_MAGNETIC_GROUP_PROPERTIES
SAVE_HSML_IN_SNAPSHOT              # this will store hsml and density values in the snapshot files 
SOFTEREQS
BLACK_HOLES               # enables Black-Holes (master switch)
BH_THERMALFEEDBACK        # quasar-mode: couple a fraction of the BH luminosity into surrounding
DRAINGAS=3                # non-stochastic smooth accretion (1: on, 2: on + cell rho, 3: on + gas drained from all cells within hsml)
BH_EXACT_INTEGRATION      # integrates analytically mass accretion
BH_BONDI_DEFAULT          # default Bondi prescription
BH_DO_NOT_PREVENT_MERGERS # When this is enabled, BHs can merge irrespective of their relative velocity 
BH_USE_ALFVEN_SPEED_IN_BONDI  # when this is enabled the alfven speed is added to the gas sound speed in the Bondi rate and the total gas pressure around the BH includes the magnetic contribution when compared the the reference pressure in BH_PRESSURE_CRITERION (requires MHD)
BH_NEW_CENTERING          # an alternative to the BH_FRICTION and REPOSITION_ON_POTMIN switches
BH_PRESSURE_CRITERION
BH_ADIOS_WIND
BH_ADIOS_WIND_WITH_QUASARTHRESHOLD  # use a threshold value ("qusarthrehold") of bondi-rate over Eddington rate to decide about quasar mode vs. adios wind
BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD  # scales the threshold with black hole mass (with a factor (M_BH/M_ref)^2, where M_ref = 10^8 Msun)  
BH_ADIOS_RANDOMIZED        # inputs momentum along alternating random directions
BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY   # disable ADIOS wind if density around blackhole drops below a certain fraction of the star formation density
GENERATE_GAS_IN_ICS
SPLIT_PARTICLE_TYPE=2+4
GFM                                    #master switch
GFM_STELLAR_EVOLUTION=0                #stellar evolution: 0->default, 1->no mass loss (beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine 
GFM_CONST_IMF=0                        #0 for Chabrier (default), 1 for a pure power-law (requires parameter IMFslope, e.g. -2.35 for Salpeter)
GFM_PREENRICH                          #pre enrich gas at given redshift
GFM_WINDS                              #decoupled ISM winds 
GFM_WINDS_VARIABLE=1                   #decoupled ISM winds: 0->scale winds with halo mass, requires FoF, 1->sigma winds
GFM_WINDS_VARIABLE_HUBBLE              #add an additional H(z)^(-1/3) factor to the wind scaling, such that it scales with halo mass not halo velocity dispersion
GFM_WIND_ENERGY_METAL_DEPENDENCE       #this can be used to decrease the wind energy for high metallicity (mimicking higher cooling losses)
GFM_WINDS_STRIPPING                    #wind metal stripping
GFM_WINDS_THERMAL_NEWDEF               #with this switch, the thermal energy is specified as a fraction of the total energy
GFM_COOLING_METAL                      #metal line cooling
GFM_AGN_RADIATION                      #cooling suppression/heating due to AGN radiation field (proximity effect)
GFM_STELLAR_PHOTOMETRICS               #calculate stellar magnitudes for different filters based on GALAXEV/BC03 
GFM_OUTPUT_MASK=1+2+4+8+16+32+64+256   #which fields to output (see io_fields.c)
GFM_NORMALIZED_METAL_ADVECTION         #this introduces an additional pseudo element for all untracked metals and normalizes the extrapolated abundance vectors to unity
GFM_OUTPUT_BIRTH_POS                   #output BirthPos and BirthVel for all star particles
GFM_CHEMTAGS                           #see documentation/README_chemical_tagging
GFM_DISCRETE_ENRICHMENT                #allow stars to enrich nearby gas from stellar evolution only above some delta mass fraction threshold
GFM_SPLITFE                            #see documentation/README_splitfe
GFM_RPROCESS                           #see documentation/README_nsns, must have GFM_SPLITFE toggled as well
