import os
import numpy as np
import json
import warnings
from scipy.interpolate import interp1d

class cosmoSim:
    """
    Class for easily handling cosmological simulation
    data products.

    Args:
        run_name (str): Run name string as found in folder names

    Attributes:
        run_name (str): Run name string as found in folder names
        redshifts (list): List of snapshot redshifts
        boxsize (int): Periodic box size in kpc
        npart (int): Cube root of number of particles
        dm_type (str): Type of dark matter in simulation
        baryon_type (str): Type of baryonic treatment in simulation
        mass_resolution (float): Mass of a DM particle in 10^10 M_sun
        softening_length (float): Force resolution in kpc
        plot_label (str): Legend label to use on plots
    """

    def __init__(self, run_name, base_path="../../data_prods/", explicit_warn=False):
        
        self.__base_path = os.path.abspath(base_path)
        
        with open(os.path.join(self.__base_path, run_name, 'run_info.json')) as f:
            run_info = json.load(f)

        self.run_name = run_info["run name"]
        self.redshifts = run_info['redshifts']

        self.file_indices = run_info['file_indices']

        self.boxsize = run_info['BoxSize']
        self.npart = run_info['NPart']

        self.dm_type = run_info['DM type']
        self.baryon_type = run_info['baryon_type']

        self.mass_resolution = run_info["Mass Resolution"]
        self.softening_length = run_info["Softening Length"]

        self.plot_label = f'{self.dm_type}'

        if self.dm_type == '2cDM' or self.dm_type == 'SIDM':
            self.sigma0 = run_info['sigma0']
            self.powerLaws = run_info['powerLaws']
            self.plot_label += f' {self.powerLaws}, $\\sigma_0=$ {self.sigma0}'
        if self.dm_type == '2cDM':
            if 'Vkick' in run_info.keys():
                self.Vkick = run_info['Vkick']
            else:
                if explicit_warn:
                    warnings.warn(f'Vkick not explicitly set in run {self.run_name}! Assuming 100 km/s...')
                self.Vkick = 100.
        if self.baryon_type == 'HY':
            self.Omega0 = run_info['Omega0']
            self.OmegaB = run_info['OmegaB']
            self.OmegaStar = run_info['OmegaStar']
        self.run_info = run_info

    def __calculate_fourier_conversion(self, Boxsize):
        '''
        Takes the Boxsize in kpc and produces the conversion factors
        from Fourier units to physical units. For use with genPK output.
        '''
        Boxsize = Boxsize / 1000 # convert to Mpc

        k_conv = 2*np.pi / Boxsize
        p_conv = (Boxsize / (2 * np.pi))**3

        return k_conv, p_conv

    def __get_genPK_data(self, fpath, boxsize):

        k_conv, p_conv = self.__calculate_fourier_conversion(boxsize)

        genPK = np.loadtxt(fpath)

        bins = genPK[:, 0]
        pk = genPK[:, 1]

        dk = pk * (2 * np.pi)**3 * (4 * np.pi) * bins**3

        bins = genPK[:, 0] * k_conv

        pk = genPK[:, 1] * p_conv

        return bins, pk, dk
    
    def __get_pyliansPK_data(self, fpath):

        pyliansPK = np.loadtxt(fpath)

        bins = pyliansPK[:, 0]
        pk = pyliansPK[:, 1]

        dk = pk * (2 * np.pi)**3 * (4 * np.pi) * bins**3

        return bins, pk, dk

    def __get_redshift_table_index(self, redshift, tolerance=0.1):
        max_tolerance = 0.5
        redshift_arr = np.array(self.redshifts)
        diffs = np.abs(redshift_arr - redshift)

        # this indexes into the redshift array
        # self.file_indices[idx] gives the actual file number
        idx = diffs.argmin()
        diff = np.amin(diffs)

        if diff > tolerance:
            warnings.warn(f"WARNING: Requested redshift {redshift} is not within tolerance {tolerance} of snapshot redshift {self.redshifts[idx]} in run {self.run_name}!")
        if diff > max_tolerance:
            raise Exception(f"ERROR: Requested redshift {redshift} is not within {max_tolerance} of any snapshot in run {self.run_name}! Snapshot may not exist!")
        return idx

    def redshift_to_index(self, redshift, tolerance=0.1):
        """
        Converts a given redshift to list index

        Args:
            redshift (float): redshift of snapshot

        Returns:
            index (int): index of file with given redshift
        """

        return self.file_indices[self.__get_redshift_table_index(redshift, tolerance)]

    def __interpolate(self, domain, range):
        """
        Creates an interpolation function for range over domain

        Args:
            domain (np.array(float)): domain of function
            range (np.array(float)): range of function

        Returns:
            interpf (function): an interpolation function
            lims (n.array(float)): the bounds of validity for the interpolation function
                            packaged in the form [inf, sup]
        """

        inf = np.amin(domain)
        sup = np.amax(domain)

        return interp1d(domain, range, fill_value="extrapolate"), np.array([inf, sup])

    def get_nearest_redshift(self, r, tolerance=0.1):
        """
        Convenience function to return the nearest redshift to the requested value

        Args:
            r (float): requested redshift
            tolerance (float): tolerance within to search

        Returns:
            nearest (float): the redshift of the snapshot with closest redshift to r
        """

        return self.redshifts[self.__get_redshift_table_index(r, tolerance)]

    def load_subhalo_info(self, redshift):
        """
        Loads tabulated subhalo information

        Args:
            redshift (float): redshift of snapshot

        Returns:
            Vmax (np.array(float)): Vmax for all subhalos
            Rmax (np.array(float)): Radius at which Vmax is achieved
            subhaloMass (np.array(float)): Subhalo mass, defined by subfind SubhaloMass
            halfMassRad (np.array(float)): Subhalo half mass radius
            massInHalfRad (np.array(float)): Mass contained within half mass radius
            massInRad (np.array(float)): Mass contained within 2 * half mass radius
        """
        idx = self.redshift_to_index(redshift)

        Vmax, Rmax, subhaloMass, halfMasRad, massInHalfRad, massInRad = np.loadtxt(
            os.path.join(
            self.__base_path,
            self.run_name,
            f'subhalo_stats_{idx}.txt')
            )

        return Vmax, Rmax, subhaloMass, halfMasRad, massInHalfRad, massInRad
    
    def load_profile_info(self, redshift):
        """
        Loads tabulated subhalo information

        Args:
            redshift (float): redshift of snapshot

        Returns:
            cutoffs (np.array(float)): Array containing lower bounds for subhalo profiles
                                       with specified number of particles, i.e. stats above
                                       a certain cutoff.
                                       Format:
                                       [ mass_cutoff, vmax_cutoff, ...mass_cutoff_by_type...]
        """
        idx = self.redshift_to_index(redshift)

        cutoffs = np.loadtxt(
            os.path.join(
            self.__base_path,
            self.run_name,
            f"profile_cutoffs_{idx}.txt")
            )

        return cutoffs

    def load_power_spectra(self, redshift, part_type='DM', backend="pylians"):
        """
        Loads tabulated power spectra from disk for this run

        Args:
            redshift (float): redshift of snapshot
            part_type (str): particle type to load
                choices are ["DM", "by", "st"]
            backend (str): program used to compute power spectrum
                choices are ["pylians", "gen-pk"]

        Returns:
            bins (np.array(float)): power spectrum bins
            pk (np.array(float)): 1D power spectrum Mpc^3/h
            dk (np.array(float)): 1D dimensionless power spectrum
        """

        if backend == 'pylians':
            if part_type == "DM":
                part_type_string = "CDM"
            elif part_type == "by":
                part_type_string = "GAS"
            elif part_type == "st":
                part_type_string = "Stars"
            elif part_type == "All":
                if self.baryon_type == "HY":
                    part_type_string = "GAS+CDM+Stars"
                else:
                    part_type_string = "CDM"
            else:
                part_type_string = ""
                raise Exception("No Particle Types Specified")
            z = self.get_nearest_redshift(redshift)
            pk_file = os.path.join(
                self.__base_path,
                self.run_name,
                f'Pk_{part_type_string}_z={z:.3f}.dat')
            bins, pk, dk = self.__get_pyliansPK_data(pk_file)

        elif backend == 'gen-pk':
            idx = self.redshift_to_index(redshift)
            pk_file = os.path.join(
                self.__base_path,
                self.run_name,
                f'PK-{part_type}-snap_{idx:03}.hdf5')
            bins, pk, dk = self.__get_genPK_data(pk_file, self.boxsize)

        k_ny = self.npart * np.pi / (self.boxsize / 1000) # k_ny in Mpc^-1

        return bins, pk, dk, k_ny

    def interp_power_spectra(self, redshift, part_type='DM', backend="pylians"):
        """
        Loads tabulated power spectra from disk and interpolates
        the result

        Args:
            redshift (float): redshift of snapshot
            part_type (str): particle type to load
                choices are ["DM", "by", "st", "all"]
            backend (str): program used to compute power spectrum
                choices are ["pylians", "gen-pk"]

        Returns:
            lims (np.array(float)): the bounds of validity for the interpolation function
                                    packaged in the form [inf, sup]
            pk_interp (function): 1D power spectrum interpolation function Mpc^3/h
            dk_interp (function): 1D dimensionless power spectrum interpolation function
        """
        bins, pk, dk, k_ny = self.load_power_spectra(redshift, part_type, backend=backend)
        pk_interp, lims = self.__interpolate(bins, pk)
        dk_interp, lims = self.__interpolate(bins, dk)

        return lims, pk_interp, dk_interp, k_ny

    def load_combined_power_spectra(self, redshift, backend="pylians"):
        """
        Combines tabulated power spectra from DM and baryonic components
        weighting by contribution to Omega0

        Args:
            redshift (float): redshift of snapshot
            backend (str): program used to compute power spectrum
                choices are ["pylians", "gen-pk"]

        Returns:
            bins (np.array(float)): power spectrum bins
            pk (np.array(float)): 1D power spectrum Mpc^3/h
            dk (np.array(float)): 1D dimensionless power spectrum
        """
        if backend == "pylians":
            bins, pk, dk, k_ny = self.load_power_spectra(redshift, part_type='All', backend=backend)

        elif backend == "gen-pk":
            bins, pk_DM, dk_DM, k_ny = self.load_power_spectra(redshift, part_type='DM', backend=backend)
            if self.baryon_type == 'DM':
                return bins, pk_DM, dk_DM, k_ny
            ridx = self.__get_redshift_table_index(redshift)
            omegaM = self.Omega0[ridx]
            try:
                _, pk_by, dk_by, _ = self.load_power_spectra(redshift, part_type='by', backend=backend)
                omegaB = self.OmegaB[ridx]
            except:
                warnings.warn(f'No gas for redshift {redshift} in run {self.run_name}')
                pk_by = dk_by = 0
                omegaB = 0
            try:
                _, pk_st, dk_st, _ = self.load_power_spectra(redshift, part_type='st', backend=backend)
                omegaStar = self.OmegaStar[ridx]
            except:
                warnings.warn(f'No stars for redshift {redshift} in run {self.run_name}')
                pk_st = dk_st = 0
                omegaStar = 0
            # weight by contribution to omegaM
            pk = omegaB * pk_by + omegaStar * pk_st + (omegaM - (omegaB + omegaStar))/omegaM * pk_DM
            dk = omegaB * dk_by + omegaStar * dk_st + (omegaM - (omegaB + omegaStar))/omegaM * dk_DM

        return bins, pk, dk, k_ny

    def interp_combined_power_spectra(self, redshift, backend="pylians"):
        """
        Combines tabulated power spectra from DM and baryonic components
        and interpolates the result

        Args:
            redshift (float): redshift of snapshot

        Returns:
            lims (np.array(float)): the bounds of validity for the interpolation function
                                    packaged in the form [inf, sup]
            pk_interp (function): 1D power spectrum interpolation function Mpc^3/h
            dk_interp (function): 1D dimensionless power spectrum interpolation function
        """
        bins, pk, dk, k_ny = self.load_combined_power_spectra(redshift, backend=backend)
        pk_interp, lims = self.__interpolate(bins, pk)
        dk_interp, lims = self.__interpolate(bins, dk)

        return lims, pk_interp, dk_interp, k_ny

    def load_mass_profile(self, redshift, partType='all'):
        """
        Loads tabulated halo mass function from disk for this run

        Args:
            redshift (float): redshift of snapshot
            partType (str): Particle type to use for mass function
                choices are ["all", "DM", "gas", "stars", "by"]

        Returns:
            mbins (np.array(float)): logarithmic mass bins in M_sun
            m (np.array(float)): cumulative halo count in bins
        """

        idx = self.redshift_to_index(redshift)

        if partType == 'all':
            mbins, m = np.loadtxt(
                os.path.join(self.__base_path,
                            self.run_name,
                            f'mass_profile_{idx}.txt')
                )
        elif partType == 'DM':
            mbins, m = np.loadtxt(
                os.path.join(self.__base_path,
                            self.run_name,
                            f'part_type_{1}_mass_profile_{idx}.txt')
                )
        elif partType == 'gas':
            mbins, m = np.loadtxt(
                os.path.join(self.__base_path,
                            self.run_name,
                            f'part_type_{0}_mass_profile_{idx}.txt')
                )
        elif partType == 'stars':
            mbins, m = np.loadtxt(
                os.path.join(self.__base_path,
                            self.run_name,
                            f'part_type_{4}_mass_profile_{idx}.txt')
                )
        elif partType == 'by':
            mbins, m = np.loadtxt(
                os.path.join(self.__base_path,
                            self.run_name,
                            f'by_mass_profile_{idx}.txt')
                )


        return mbins, m

    def interp_mass_profile(self, redshift, partType='all'):
        """
        Loads tabulated halo mass function from disk
        and interpolates the result

        Args:
            redshift (float): redshift of snapshot

        Returns:
            lims (np.array(float)): the bounds of validity for the interpolation function
                                    packaged in the form [inf, sup]
            m_interp (function): mass function interpolation function
        """
        mbins, m = self.load_mass_profile(redshift, partType=partType)

        m_interp, lims = self.__interpolate(mbins, m)

        return lims, m_interp

    def load_vel_profile(self, redshift):
        """
        Loads tabulated circular velocity function from disk for this run

        Args:
            redshift (float): redshift of snapshot

        Returns:
            vbins (np.array(float)): logarithmic velocity bins in km/s
            v (np.array(float)): cumulative halo count in bins
        """

        idx = self.redshift_to_index(redshift)

        vbins, v = np.loadtxt(
            os.path.join(self.__base_path,
                         self.run_name,
                         f'vel_profile_{idx}.txt')
            )

        return vbins, v
        
    def interp_vels_profile(self, redshift):
        """
        Loads tabulated circular velocity function from disk
        and interpolates the result

        Args:
            redshift (float): redshift of snapshot

        Returns:
            lims (np.array(float)): the bounds of validity for the interpolation function
                                    packaged in the form [inf, sup]
            v_interp (function): velocity function interpolation function
        """
        vbins, v = self.load_vel_profile(redshift)

        v_interp, lims = self.__interpolate(vbins, v)

        return lims, v_interp
    
    # TODO: Unify these load/interp functions into just one with options

    def load_SFR_profile(self, redshift):
        """
        Loads tabulated SFR function from disk for this run

        Args:
            redshift (float): redshift of snapshot

        Returns:
            sfr_bins (np.array(float)): logarithmic SFR bins in M_sun/yr
            s (np.array(float)): cumulative halo count in bins
        """

        idx = self.redshift_to_index(redshift)

        sfr_bins, s = np.loadtxt(
            os.path.join(self.__base_path,
                         self.run_name,
                         f'SFR_profile_{idx}.txt')
            )

        return sfr_bins, s

    def interp_SFR_profile(self, redshift):
        """
        Loads tabulated SFR function from disk
        and interpolates the result

        Args:
            redshift (float): redshift of snapshot

        Returns:
            lims (np.array(float)): the bounds of validity for the interpolation function
                                    packaged in the form [inf, sup]
            s_interp (function): SFR function interpolation function
        """
        sfr_bins, s = self.load_SFR_profile(redshift)

        s_interp, lims = self.__interpolate(sfr_bins, s)

        return lims, s_interp

    def load_mass_density(self, redshift, subhalo_idx, partType='all'):
        """
        Loads tabulated halo mass density profile from disk

        Args:
            redshift (float): redshift of snapshot
            subhalo_idx (int): index of subhalo with 0 being largest,
                               1 second largest, etc.
            partType (str):
                choices are ["all", "gas", "by", "DM", "st", "bh"]
            Returns:
                rbins (np.array(float)): radius bins densities were calculated within
                densities (np.array(float)): Mass density within each bin in 10e10 M_sun/kpc^3
        """

        idx = self.redshift_to_index(redshift)

        X = np.loadtxt(
            os.path.join(self.__base_path,
                         self.run_name,
                         f'subhalo_densities_{idx}',
                         f'subhalo_{subhalo_idx}.txt')
            )
        rbins = X[0]
        if partType == 'all':
            densities = X[1]
        elif partType == 'gas':
            densities = X[2]
        elif partType == 'DM':
            densities = X[3]
        elif partType == 'st':
            densities = X[4]
        elif partType == 'bh':
            densities = X[5]
        elif partType == 'by':
            densities = X[2] + X[4]

        return rbins, densities
    
    def get_fake_spectra_savefile_path(self, redshift):
        """
        Outputs the folder path to `fake_spectra` generated savefiles
        
        Args:
            redshift (float):
                redshift of snapshot
        Returns:
            base (str):
                absolute path to the folder containing fake_spectra savefile folders
            num (int):
                index for this requested redshift
        """
        
        num = self.redshift_to_index(redshift)
        base = str(os.path.abspath(os.path.join(self.__base_path, self.run_name)))
        return num, base