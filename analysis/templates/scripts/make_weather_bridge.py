"""Export only restored weather days for preparing the same panel on server."""
from pathlib import Path
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
d=pd.read_parquet(BASE/'data'/'br_calendar.parquet',columns=['unit','date','tmean','tmax','population','y_K80'])
extra=d.loc[d.y_K80.isna(),['unit','date','tmean','tmax','population']]
h=pd.read_csv(r'F:\分析专用\climate_attribution_20261106\work\expo_temp\region_day_temp_2005.csv',dtype={'hrid':str},parse_dates=['date'])
h=h.loc[h.date>='2005-12-12',['hrid','date','tmean','tmax']].rename(columns={'hrid':'unit'})
bridge=pd.concat([extra,h],ignore_index=True)
assert not bridge.duplicated(['unit','date']).any()
bridge.to_csv(BASE/'data'/'br_weather_bridge.csv',index=False)
print('weather-only and prehistory rows',len(bridge))
