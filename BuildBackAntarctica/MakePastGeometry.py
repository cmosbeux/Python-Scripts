import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import h5py
from scipy.interpolate import griddata

# Charger le dataset netcdf
file_path =  '/Users/cmosbeux/Documents/Data_Antarctica.CLEAN/DEMs&dHdT/BedMachineAntarctica-v03.nc'
ds = xr.open_dataset(file_path)

# Récupérer les variables du dataset
thick = ds['thickness'].values
surf = ds['surface'].values
bed = ds['bed'].values
x_thick = ds['x'].values  # Coordonnées x
y_thick = ds['y'].values  # Coordonnées y

# Charger le fichier HDF5 contenant dhdt ground
with h5py.File('/Users/cmosbeux/Documents/Data/ICESat/AIS_mass_change.h5', 'r') as f_ground:
    dhdt_ground = f_ground['/dHdt'][()]
    x_data_g = f_ground['x'][:]
    y_data_g = f_ground['y'][:]

# Charger le fichier HDF5 contenant dhdt IS
with h5py.File('/Users/cmosbeux/Documents/Data_Antarctica.CLEAN/DEMs&dHdT/ICESat2/ICE1_ICE2_AnIS_dHdt_2003_2018_R209_05KM_FLOAT_MASS_F2.h5', 'r') as f_IS:
    dhdt_IS = f_IS['/dhdt'][()]
    x_data_IS = f_IS['x'][:]
    y_data_IS = f_IS['y'][:]
    
print('Files loaded...')

# Créer une grille 2D des coordonnées
lon_thick, lat_thick = np.meshgrid(x_thick, y_thick)

# Condition : glacier posé ou flottant
is_grounded = (thick == (surf - bed))  # Vrai si le glacier est posé
is_floating = ~is_grounded  # Le reste est flottant


print('Interpolation for grounded ice...')
# Interpolation pour les glaciers posés
thick_ground_interpolated = griddata(
    (lon_thick.flatten(), lat_thick.flatten()),
    thick.flatten(),
    (x_data_g, y_data_g),
    method='linear'
)

# Calculer h pour les glaciers posés
h_ground = thick_ground_interpolated - dhdt_ground

print('Interpolation for floating ice...')
# Interpolation pour les glaciers flottants
thick_float_interpolated = griddata(
    (lon_thick.flatten(), lat_thick.flatten()),
    thick.flatten(),
    (x_data_IS, y_data_IS),
    method='linear'
)

# Calculer h pour les glaciers flottants
h_float = thick_float_interpolated - dhdt_IS

# Masquer les zones posées et flottantes dans une nouvelle variable h
h_combined = np.where(is_grounded, h_ground, h_float)

# Afficher le résultat
plt.imshow(h_combined, cmap='viridis', extent=[x_thick.min(), x_thick.max(), y_thick.min(), y_thick.max()])
plt.colorbar(label='Épaisseur (m)')
plt.title('Différence d\'épaisseur (h)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()
plt.savefig('diff_thik')
