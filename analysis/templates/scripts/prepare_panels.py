"""Prepare unique daily panels and strict calendar weather histories for R."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / 'data' / 'source'
OUT = BASE / 'data'
MX_WEATHER = Path(r'D:\数据库集合\新增数据\ERA5_crossnational\mex_state_day_2015_2023.csv')
BR_WEATHER = Path(r'F:\分析专用\climate_attribution_20261106\work\expo_temp')

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def finish(d, label, start, end, history_start=None):
    d['date'] = pd.to_datetime(d.date)
    d['unit'] = d.unit.astype(str)
    if d.duplicated(['unit','date']).any():
        raise ValueError(f'{label}: duplicate unit-date')
    pieces = []
    audit = []
    for unit, s in d.groupby('unit', sort=True):
        s = s.set_index('date').sort_index()
        s = s.reindex(pd.date_range(history_start or start, end, freq='D'))
        s.index.name = 'date'
        s['unit'] = unit
        for k in [3,7,14,21]:
            s[f'E{k}'] = s.tmean.rolling(k, min_periods=k).mean()
        s['Emax3'] = s.tmax.rolling(3,min_periods=3).mean()
        s = s.loc[pd.Timestamp(start):pd.Timestamp(end)].copy()
        s['dow'] = s.index.dayofweek + 1
        s['month'] = s.index.month
        s['year'] = s.index.year
        # Population is an annual denominator. Existing unknown outcome days remain NA.
        if label == 'BR':
            s['population'] = s.groupby('year').population.transform(lambda x: x.ffill().bfill())
        s['country'] = label
        eligible = s[['y_K80','y_K81','y_urinary','E3']].notna().all(axis=1)
        audit.append({'unit':unit,'calendar_days':len(s),'weather_missing':int(s.tmean.isna().sum()),
                      'outcome_missing':int(s.y_K80.isna().sum()),'E3_valid':int(s.E3.notna().sum()),
                      'primary_eligible_days':int(eligible.sum()),
                      'P95':float(s.E3.quantile(.95)),
                      'median':float(s.E3.median())})
        pieces.append(s.reset_index())
    result = pd.concat(pieces,ignore_index=True)
    result.to_csv(OUT/f'{label.lower()}_calendar.csv',index=False)
    result.to_csv(OUT/f'{label.lower()}_calendar.csv.gz',index=False,compression='gzip')
    result.to_parquet(OUT/f'{label.lower()}_calendar.parquet',index=False)
    pd.DataFrame(audit).to_csv(OUT/f'{label.lower()}_unit_coverage.csv',index=False)
    return {'calendar_rows':len(result),'units':len(audit),'eligible_days':int(sum(a['primary_eligible_days'] for a in audit)),
            'missing_outcomes':int(result.y_K80.isna().sum()),'missing_weather':int(result.tmean.isna().sum()),
            'outcome_totals':{c:int(result[c].sum()) for c in ['y_K80','y_K81','y_urinary']}}

br_path=SOURCE/'brazil_panel.parquet'
mx_path=SOURCE/'mex_state_day_stone.csv'
br=pd.read_parquet(br_path).rename(columns={'hrid':'unit','K80':'y_K80','K81':'y_K81','N20_N23':'y_urinary',
                                          'sum3':'y_sum3','pop_region':'population'})
cols=['unit','date','tmean','tmax','population','y_K80','y_K81','y_urinary','y_sum3']
br=br[cols]
mx=pd.read_csv(mx_path,dtype={'state':str}).rename(columns={'state':'unit','k80':'y_K80','k81':'y_K81','n20':'y_urinary'})
mx['unit']=mx.unit.str.zfill(2)
mx['date']=pd.to_datetime(mx.date)
source_info={'br':{'sha256':sha(br_path)},'mx':{'sha256':sha(mx_path)}}
weather_files=[BR_WEATHER/f'region_day_temp_{y}.csv' for y in range(2005,2021)]
if all(p.is_file() for p in weather_files):
    wb=pd.concat([pd.read_csv(p,dtype={'hrid':str},parse_dates=['date']) for p in weather_files],ignore_index=True)
    wb=wb.rename(columns={'hrid':'unit'})[['unit','date','tmean','tmax']]
    wb=wb[wb.date.between('2005-12-12','2020-12-31')]
    br['date']=pd.to_datetime(br.date)
    br['unit']=br.unit.astype(str)
    check=br[['unit','date','tmean','tmax']].merge(wb,on=['unit','date'],suffixes=('_frozen','_weather'),validate='one_to_one')
    differences={c:float((check[c+'_frozen']-check[c+'_weather']).abs().max()) for c in ['tmean','tmax']}
    if len(check)!=len(br) or max(differences.values())>2e-5:
        raise ValueError(f'Brazil source weather does not match frozen panel: {len(check)}, {differences}')
    # Preserve original frozen exposure on existing rows; restore only genuinely available missing calendar weather.
    weather_only=wb.merge(br[['unit','date']],on=['unit','date'],how='left',indicator=True)
    weather_only=weather_only[weather_only['_merge']=='left_only'].drop(columns='_merge')
    br=pd.concat([br,weather_only],ignore_index=True)
    source_info['br']['weather_files']=[{'path':str(p),'sha256':sha(p)} for p in weather_files]
    source_info['br']['weather_agreement_max_abs']=differences
    source_info['br']['added_weather_only_days']=len(weather_only)
if MX_WEATHER.is_file():
    w=pd.read_csv(MX_WEATHER,dtype={'state':str})
    w['unit']=w.state.str.zfill(2)
    w['date']=pd.to_datetime(dict(year=w.year,month=w.month,day=w.day))
    w=w[['unit','date','tmean']]
    w=w[(w.unit!='00')&(w.date.between('2015-01-01','2022-12-31'))]
    if w.duplicated(['unit','date']).any():
        raise ValueError('Mexico full weather has duplicate dates')
    check=mx[['unit','date','tmean']].merge(w,on=['unit','date'],suffixes=('_frozen','_weather'),validate='one_to_one')
    delta=(check.tmean_frozen-check.tmean_weather).abs()
    if delta.max()>1e-6:
        raise ValueError(f'Mexico weather inconsistency {delta.max()}')
    source_info['mx']['weather_sha256']=sha(MX_WEATHER)
    source_info['mx']['weather_matched_frozen_days']=len(check)
    mx=w.merge(mx.drop(columns='tmean'),on=['unit','date'],how='outer',validate='one_to_one')
mx['tmax']=np.nan
mx['population']=np.nan
mx['y_sum3']=mx[['y_K80','y_K81','y_urinary']].sum(axis=1,min_count=3)
source_info['br']['prepared']=finish(br,'BR','2006-01-01','2020-12-31','2005-12-12' if 'weather_files' in source_info['br'] else None)
source_info['mx']['prepared']=finish(mx[cols],'MX','2015-01-01','2022-12-31')
(BASE/'admin'/'panel_preparation.json').write_text(json.dumps(source_info,indent=2),encoding='utf-8')
print(json.dumps(source_info,indent=2))
