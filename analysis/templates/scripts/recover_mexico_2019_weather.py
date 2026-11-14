"""Recover matching 2019 state weather and validate against an existing month.

The spatial mask, cosine latitude weights, valid-land mask, daily UTC mean,
Kelvin conversion and valid-weight denominator match build_mex_state_day.py.
Shapely's vectorized contains_xy replaces repeated Point construction.
"""
from pathlib import Path
import calendar
import hashlib
import json
import warnings
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr

BASE=Path(__file__).resolve().parents[1]
RAW=Path(r'D:\数据库集合\新增数据\ERA5_crossnational\mexico')
GEO=Path(r'D:\数据库集合\新增数据\ERA5_crossnational\mex_states.geojson')
OLD=Path(r'D:\数据库集合\新增数据\ERA5_crossnational\mex_state_day_2015_2023.csv')
SOURCE_CODE=Path(__file__).resolve().parent / "build_mex_state_day.py"
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            h.update(b)
    return h.hexdigest()

states=gpd.read_file(GEO).sort_values('cve_ent').reset_index(drop=True)
assert len(states)==32
ref_path=sorted(RAW.glob('*.nc'))[0]
with xr.open_dataset(ref_path,engine='netcdf4') as ds:
    lat=ds.latitude.values.astype(float)
    lon=ds.longitude.values.astype(float)
    valid=np.isfinite(ds.t2m.isel(valid_time=0).values).ravel()
xx,yy=np.meshgrid(lon,lat)
area=np.cos(np.deg2rad(lat))[:,None]*np.ones((1,len(lon)),dtype=np.float32)
W=np.zeros((32,len(lat)*len(lon)),dtype=np.float32)
for i,geom in enumerate(states.geometry):
    inside=shapely.contains_xy(geom,xx,yy).ravel()&valid
    W[i,inside]=area.ravel()[inside]
s=W.sum(axis=1,keepdims=True)
if (s<=0).any():
    raise ValueError('a state has no valid grid cells')
W/=s
audit={'original_source_code':str(SOURCE_CODE),'source_code_sha256':sha(SOURCE_CODE),
       'geography_sha256':sha(GEO),'reference_nc':str(ref_path),'files':[],
       'spatial_rule':'grid-center within state, cosine latitude area weights, fixed reference valid-land mask, daily valid-weight renormalization'}
def month(year,mm):
    path=RAW/f'era5land_hourly_temp_dewpoint_mexico_{year}_{mm:02d}.nc'
    with xr.open_dataset(path,engine='netcdf4') as ds:
        if not np.array_equal(ds.latitude.values,lat) or not np.array_equal(ds.longitude.values,lon):
            raise ValueError(f'grid mismatch in {path.name}')
        dates=pd.to_datetime(ds.valid_time.values)
        expected=pd.date_range(f'{year}-{mm:02d}-01',periods=calendar.monthrange(year,mm)[1]*24,freq='h')
        if not dates.equals(expected):
            raise ValueError(f'incomplete/nonhourly time in {path.name}')
        t=ds.t2m.values-273.15
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        daily=np.nanmean(t.reshape(-1,24,len(lat),len(lon)),axis=1)
    rows=[]
    for day,field in zip(pd.date_range(f'{year}-{mm:02d}-01',periods=len(daily)),daily):
        flat=field.ravel()
        for i,st in enumerate(states.cve_ent):
            mask=W[i]>0
            numerator=np.nansum(flat[mask]*W[i,mask])
            denominator=np.nansum(W[i,mask]*np.isfinite(flat[mask]))
            value=float(numerator/denominator) if denominator>0 else np.nan
            rows.append((str(st).zfill(2),year,mm,day.day,value))
    audit['files'].append({'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path),'hours':len(dates)})
    print(f'{year}-{mm:02d} complete: {len(rows)} state-days',flush=True)
    return pd.DataFrame(rows,columns=['state','year','month','day','tmean'])

validation=month(2018,1)
old=pd.read_csv(OLD,dtype={'state':str})
old.state=old.state.str.zfill(2)
comparison=validation.merge(old,on=['state','year','month','day'],suffixes=('_new','_old'),validate='one_to_one')
maxdiff=float((comparison.tmean_new-comparison.tmean_old).abs().max())
if len(comparison)!=32*31 or maxdiff>1e-6:
    raise ValueError(f'weather validation failed {len(comparison)} {maxdiff}')
audit['validation']={'month':'2018-01','matched_days':len(comparison),'max_abs_difference':maxdiff}
frames=[]
for mm in range(1,13):
    cache=BASE/'data'/f'mx_weather_2019_{mm:02d}.csv'
    if cache.is_file():
        frames.append(pd.read_csv(cache,dtype={'state':str}))
    else:
        frame=month(2019,mm)
        frame.to_csv(cache,index=False)
        frames.append(frame)
out=pd.concat(frames,ignore_index=True)
if len(out)!=32*365 or out.tmean.isna().any() or out.duplicated(['state','year','month','day']).any():
    raise ValueError('2019 state weather coverage failed')
out.to_csv(BASE/'data'/'mx_weather_2019.csv',index=False)
audit['recovered_state_days']=len(out)
(BASE/'admin'/'mx_weather_2019_recovery.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
print(json.dumps({'validation':audit['validation'],'recovered_state_days':len(out)},indent=2))
