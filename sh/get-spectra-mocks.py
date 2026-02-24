import glob
import sys

from fake_spectra.randspectra import RandSpectra

IN_PATH = sys.argv[1]

for i in range(len(glob.glob(f"{IN_PATH}/snap*"))):
    rr = RandSpectra(i, IN_PATH, thresh=0.)
    rr.get_tau("H",1,1215)
    #Lyman-beta
    rr.get_tau("H",1,1025)
    rr.get_col_density("H",1)
    #Save spectra to file
    rr.save_file()