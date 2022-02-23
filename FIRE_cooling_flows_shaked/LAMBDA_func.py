import h5py, numpy as np
from numpy import log10 as log
import scipy
import scipy.stats
from scipy import interpolate

def LAMBDA(fn,Y=0.28):
    f=h5py.File(fn,'r')
    iHe = np.searchsorted(f['Metal_free']['Helium_mass_fraction_bins'][:],Y)    
    H_He_Cooling  = f['Metal_free']['Net_Cooling'][iHe,...]
    Tbins         = f['Metal_free']['Temperature_bins'][...]
    nHbins        = f['Metal_free']['Hydrogen_density_bins'][...]
    Metal_Cooling = f['Total_Metals']['Net_cooling'][...] 
    
    f_H_He = interpolate.RegularGridInterpolator((log(Tbins), log(nHbins)),
                                                    H_He_Cooling, 
                                                    bounds_error=False, fill_value=None)
    f_Z = interpolate.RegularGridInterpolator((log(Tbins), log(nHbins)),
                                                    Metal_Cooling, 
                                                    bounds_error=False, fill_value=None)
    return lambda T,nH,Z2Zsun,f_H_He=f_H_He,f_Z=f_Z: (
        f_H_He((log(T), log(nH))) + f_Z((log(T), log(nH))) * Z2Zsun )

