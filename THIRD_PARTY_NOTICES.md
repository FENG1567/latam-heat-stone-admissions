# Third-party sources and software

Brazil administrative data originated from SIH-SUS/DATASUS; Mexican data from DGIS/SAEH. Weather inputs are derived from ECMWF/Copernicus ERA5-Land. Geographic mappings and source dictionaries remain attributable to their original providers. The repository's MIT and CC BY grants do not override those providers' terms. Original row-level records, monthly NetCDF fields and identifying administrative fields are not bundled.

Public source portals: https://datasus.saude.gov.br/ ; https://www.dgis.salud.gob.mx/ ; https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land . Source documentation and the paper identify selection years and limitations; portal contents may change. Historical Mexican retrieval requests were not retained, which limits exact raw-download reconstruction.

R, NumPy, Matplotlib, Pillow, data.table, jsonlite, sandwich, metafor and their dependencies retain their respective software licenses. They are installed through their package managers rather than vendored here.
