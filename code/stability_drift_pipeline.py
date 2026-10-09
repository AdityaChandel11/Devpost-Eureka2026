"""
Expired by Climate — stability-zone margin audit pipeline
=========================================================

Recomputes the climate statistics that pharmaceutical stability "climatic zones"
were built on (mean kinetic temperature, MKT, and mean water-vapour partial
pressure, Pd) from ERA5 reanalysis, year by year, and measures how much of the
margin between today's climate and each long-term stability TEST condition has
been consumed.

Test conditions audited (long-term):
    Zone II   25 C / 60 % RH   -> Pd 19.0 hPa
    Zone IVa  30 C / 65 % RH   -> Pd 27.6 hPa
    Zone IVb  30 C / 75 % RH   -> Pd 31.8 hPa
    (Pd values computed with the Buck (1981) saturation formula, see es_hpa)

Outputs per grid cell and year:
    mkt_c      MKT of open-air 2 m temperature (hourly diurnal cycle per month)
    mkt_grimm  MKT with the Grimm-style 19 C floor (values below 19 C set to 19 C,
               mimicking heated indoor storage) -- the comparator closest to how the
               zones were originally derived. VERIFY the exact Grimm procedure
               against Grimm (1998) before final reporting.
    pd_hpa     mean water-vapour partial pressure from 2 m dewpoint
    margin_*   test temperature minus MKT (deg C) and test Pd minus Pd (hPa);
               negative = the climate now exceeds the test condition.

Data: ERA5 "monthly averaged reanalysis by hour of day" (24 values per month),
which keeps the diurnal cycle that MKT is sensitive to while staying small.
Requires a free Copernicus CDS account and ~/.cdsapirc (see
https://cds.climate.copernicus.eu/how-to-api).  NOT run in the session that
wrote it -- start with one year and a small area to check variable names.

Usage:
    pip install cdsapi xarray netCDF4 numpy pandas
    python stability_drift_pipeline.py selftest
    python stability_drift_pipeline.py download 1980 2025 --area 38 66 6 98   # N W S E (India box)
    python stability_drift_pipeline.py analyse era5_hourofday_*.nc --out results.nc
"""
import sys
import argparse
import numpy as np

R = 8.3144598          # J/(mol K)
DH_DEFAULT = 83144.0   # J/mol, the default activation energy used for MKT (USP/ICH convention)

TEST_CONDITIONS = {    # name: (temperature C, RH %)
    "II": (25.0, 60.0),
    "IVa": (30.0, 65.0),
    "IVb": (30.0, 75.0),
}


def es_hpa(t_c):
    """Saturation vapour pressure over water, hPa (Buck 1981)."""
    t_c = np.asarray(t_c, dtype=float)
    return 6.1121 * np.exp((18.678 - t_c / 234.5) * (t_c / (257.14 + t_c)))


def mkt_c(temps_c, axis=-1, dh=DH_DEFAULT, floor_c=None):
    """Mean kinetic temperature (Haynes 1971) of equally weighted samples, deg C."""
    t = np.asarray(temps_c, dtype=float)
    if floor_c is not None:
        t = np.maximum(t, floor_c)
    tk = t + 273.15
    mean_exp = np.mean(np.exp(-dh / (R * tk)), axis=axis)
    return (dh / R) / (-np.log(mean_exp)) - 273.15


def rate_ratio(t_ref_c, t_c, dh=DH_DEFAULT):
    """Arrhenius degradation-rate ratio k(t)/k(t_ref)."""
    return np.exp(dh / R * (1.0 / (t_ref_c + 273.15) - 1.0 / (np.asarray(t_c) + 273.15)))


def effective_shelf_life(months_at_test, test_c, mkt_actual_c, dh=DH_DEFAULT):
    """Zero-order approximation: shelf life consumed in proportion to rate."""
    return months_at_test / rate_ratio(test_c, mkt_actual_c, dh)


def margins(mkt_value_c, pd_value_hpa):
    out = {}
    for zone, (t, rh) in TEST_CONDITIONS.items():
        out[f"margin_T_{zone}"] = t - mkt_value_c
        out[f"margin_Pd_{zone}"] = es_hpa(t) * rh / 100.0 - pd_value_hpa
    return out


# ----------------------------------------------------------------- download
def download(y0, y1, area):
    import cdsapi
    c = cdsapi.Client()
    for year in range(y0, y1 + 1):
        target = f"era5_hourofday_{year}.nc"
        c.retrieve(
            "reanalysis-era5-single-levels-monthly-means",
            {
                "product_type": ["monthly_averaged_reanalysis_by_hour_of_day"],
                "variable": ["2m_temperature", "2m_dewpoint_temperature"],
                "year": [str(year)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "time": [f"{h:02d}:00" for h in range(24)],
                "area": area,                     # [N, W, S, E]
                "data_format": "netcdf",
            },
            target,
        )
        print("saved", target)


# ----------------------------------------------------------------- analyse
def analyse(paths, out):
    import xarray as xr
    ds = xr.open_mfdataset(paths, combine="by_coords")
    tname = "valid_time" if "valid_time" in ds.dims else "time"
    t2m = ds["t2m"] - 273.15
    d2m = ds["d2m"] - 273.15
    year = ds[tname].dt.year
    res = []
    for y in np.unique(year.values):
        sel = {tname: year == y}
        t = t2m.sel(sel)
        d = d2m.sel(sel)
        mkt = xr.apply_ufunc(mkt_c, t, input_core_dims=[[tname]], kwargs={"axis": -1})
        mktg = xr.apply_ufunc(mkt_c, t, input_core_dims=[[tname]],
                              kwargs={"axis": -1, "floor_c": 19.0})
        pd_ = xr.apply_ufunc(es_hpa, d).mean(tname)
        yr = xr.Dataset({"mkt_c": mkt, "mkt_grimm": mktg, "pd_hpa": pd_})
        for k, v in margins(yr["mkt_grimm"], yr["pd_hpa"]).items():
            yr[k] = v
        res.append(yr.expand_dims(year=[int(y)]))
    xr.concat(res, "year").to_netcdf(out)
    print("wrote", out)
    # Next steps (see doc, Section 4): regrid population (GPWv4 / WorldPop) to the
    # ERA5 grid, compute population living where margin < 0, fit per-cell trends
    # (Theil-Sen + Mann-Kendall), and validate ERA5 against station data.


# ----------------------------------------------------------------- self-test
def selftest():
    # 1. Reproduce the worked example used in ICH training material: 20 C and 40 C -> MKT ~34.4 C
    m = mkt_c([20.0, 40.0])
    assert abs(m - 34.4) < 0.1, m
    # 2. Constant temperature -> MKT equals that temperature
    assert abs(mkt_c([30.0] * 24) - 30.0) < 1e-9
    # 3. Floor works
    assert abs(mkt_c([10.0, 10.0], floor_c=19.0) - 19.0) < 1e-9
    # 4. Test-condition partial pressures
    pd = {z: es_hpa(t) * rh / 100 for z, (t, rh) in TEST_CONDITIONS.items()}
    assert abs(pd["IVb"] - 31.8) < 0.1 and abs(pd["IVa"] - 27.6) < 0.1 and abs(pd["II"] - 19.0) < 0.1
    # 5. Synthetic hot-humid site: diurnal 27-37 C, dewpoint 25.5 C
    hours = np.arange(24)
    t = 32 + 5 * np.sin((hours - 9) / 24 * 2 * np.pi)
    m_site = mkt_c(t)
    p_site = es_hpa(25.5)
    mg = margins(m_site, p_site)
    print(f"synthetic site: MKT {m_site:.2f} C, Pd {p_site:.1f} hPa, "
          f"IVb margins T {mg['margin_T_IVb']:+.2f} C, Pd {mg['margin_Pd_IVb']:+.2f} hPa, "
          f"24-month label -> {effective_shelf_life(24, 30, m_site):.1f} months")
    print("selftest OK")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    d = sub.add_parser("download")
    d.add_argument("y0", type=int)
    d.add_argument("y1", type=int)
    d.add_argument("--area", type=float, nargs=4, default=[90, -180, -90, 180])
    a = sub.add_parser("analyse")
    a.add_argument("paths", nargs="+")
    a.add_argument("--out", default="stability_drift_results.nc")
    args = p.parse_args()
    if args.cmd == "selftest":
        selftest()
    elif args.cmd == "download":
        download(args.y0, args.y1, args.area)
    else:
        analyse(args.paths, args.out)
