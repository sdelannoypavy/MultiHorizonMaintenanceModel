import netCDF4
import numpy as np
import os
from glob import glob
import shutil
#
# Répertoire contenant les fichiers NetCDF
directory = '/home/delannoypavysol/MultiHorizonMaintenanceModel/data/'
# Motif pour récupérer tous les fichiers NetCDF
files = glob(os.path.join(directory, '*.nc'))
output_file = 'data_merged.nc'

def merge_files(file_2023, file_2024):

    # Ouvrir les fichiers NetCDF en mode lecture
    file_2023_nc = netCDF4.Dataset(file_2023, 'r')
    file_2024_nc = netCDF4.Dataset(file_2024, 'r')

    # Vérification des dimensions 'latitude' et 'longitude' dans les deux fichiers
    lat_2023 = file_2023_nc.variables['latitude'][:]
    lon_2023 = file_2023_nc.variables['longitude'][:]
    lat_2024 = file_2024_nc.variables['latitude'][:]
    lon_2024 = file_2024_nc.variables['longitude'][:]

    valid_time_2023 = file_2023_nc.variables['valid_time'][:]
    valid_time_2024 = file_2024_nc.variables['valid_time'][:]

    var_name = 'swh'
    temperature_2023 = file_2023_nc.variables[var_name][:]
    temperature_2024 = file_2024_nc.variables[var_name][:]

    file_2023_nc.close()
    file_2024_nc.close()

    # Créer un fichier de sortie NetCDF
    output_nc = netCDF4.Dataset(output_file, 'w')

    # Copier les dimensions 'latitude', 'longitude' et 'valid_time'
    output_nc.createDimension('latitude', len(lat_2023))
    output_nc.createDimension('longitude', len(lon_2023))

    # Fusionner les dimensions 'valid_time' (la taille de cette dimension va augmenter)
    output_nc.createDimension('valid_time', len(valid_time_2023) + len(valid_time_2024))

    # Créer les variables dans le fichier de sortie pour 'latitude', 'longitude' et 'valid_time'
    output_lat = output_nc.createVariable('latitude', 'f4', ('latitude',))
    output_lon = output_nc.createVariable('longitude', 'f4', ('longitude',))
    output_valid_time = output_nc.createVariable('valid_time', 'f4', ('valid_time',))

    # Copier les données des dimensions
    output_lat[:] = lat_2023
    output_lon[:] = lon_2023
    output_valid_time[:] = np.concatenate((valid_time_2023, valid_time_2024), axis=0)

    # Copier les variables (par exemple, 'temperature' ou d'autres variables) en fusionnant sur 'valid_time'
    # !!!!!!!!!!!!!!!!!!!!!!!A changer si on change de variable!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Ouvrir et copier les données de la variable de chaque fichier
    # Créer la variable de sortie et la remplir
    output_temp = output_nc.createVariable(var_name, 'f4', ('valid_time', 'latitude', 'longitude'))
    output_temp[:] = np.concatenate((temperature_2023, temperature_2024), axis=0)

    # Fermer les fichiers
    output_nc.close()

    print(f"Les fichiers ont été fusionnés avec succès dans {output_file}.")

# Chemins des fichiers à fusionner

for i, file_path in enumerate(files):

    print(i)

    if i == 0:
        shutil.copy(file_path, output_file)

    if i > 0:
        merge_files(output_file, file_path)