import os
import numpy as np
import json

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
    __base_path = "../../data_prods/"

    def __init__(self, run_name):

        with open(os.path.join(self.__base_path, run_name, 'run_info.json')) as f:
            run_info = json.load(f)

        self.run_name = run_info["run name"]
        self.redshifts = run_info['redshifts']
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
            self.plot_label += f' {self.powerLaws}, $\sigma_0=$ {self.sigma0}'

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
    
    def load_power_spectra(self, redshift, part_type='DM'):
        """
        Loads tabulated power spectra from disk for this run

        Args:
            redshift (float): redshift of snapshot
            part_type (str): particle type to load
                choices are ["DM", "by"] 

        Returns:
            bins (np.array(float)): power spectrum bins
            pk (np.array(float)): 1D power spectrum Mpc^3/h
            dk (np.array(float)): 1D dimensionless power spectrum
        """
        pk_file = os.path.join(
            self.__base_path, 
            self.run_name, 
            f'PK-{part_type}-snap_{self.redshifts.index(redshift):03}.hdf5') 
        
        bins, pk, dk = self.__get_genPK_data(pk_file, self.boxsize)

        k_ny = self.npart * np.pi / (self.boxsize / 1000) # k_ny in Mpc^-1

        return bins, pk, dk, k_ny
    
    def load_mass_profile(self, redshift):
        """
        Loads tabulated halo mass function from disk for this run

        Args:
            redshift (float): redshift of snapshot

        Returns:
            mbins (np.array(float)): logarithmic mass bins in M_sun
            m (np.array(float)): cumulative halo count in bins
        """
        mbins, m = np.loadtxt( 
            os.path.join(self.__base_path, 
                         self.run_name, 
                         f'mass_profile_{self.redshifts.index(redshift)}.txt') 
            )
        
        return mbins, m
    
    def load_vel_profile(self, redshift):
        """
        Loads tabulated circular velocity function from disk for this run

        Args:
            redshift (float): redshift of snapshot

        Returns:
            vbins (np.array(float)): logarithmic velocity bins in km/s
            v (np.array(float)): cumulative halo count in bins
        """
        vbins, v = np.loadtxt( 
            os.path.join(self.__base_path, 
                         self.run_name, 
                         f'vel_profile_{self.redshifts.index(redshift)}.txt') 
            )
        
        return vbins, v