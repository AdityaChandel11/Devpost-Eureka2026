"""Compute every number used in the paper from the merged CSV, and draw the data figures.
Run: python3 analyse.py  -> writes results.json and figs/fig2_margins.png, figs/fig3_humidity.png
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BIAS = 0.3          # daily-max/min reconstruction over-reads MKT by ~0.3 C vs hourly data (NOAA, 2024)
R = 8.3144598
df = pd.read_csv("data/expired_by_climate_results.csv")
df["is_era5"] = df.source.str.startswith("ERA5")

# Primary series per place: ERA5 where available, otherwise the station record.
primary = []
for city, g in df.groupby("city"):
    g = g[g.is_era5] if g.is_era5.any() else g
    primary.append(g)
prim = pd.concat(primary)

rng = np.random.default_rng(2026)
rows = []
for city, g in prim.groupby("city"):
    zone = g.zone.iloc[0]
    T = 30.0 if zone == "IV" else 25.0
    p1, p2 = g[g.year <= 2000], g[g.year >= 2006]
    if len(p1) < 15 or len(p2) < 15:
        continue
    d = p2.mkt_grimm_c.mean() - p1.mkt_grimm_c.mean()
    bs = [rng.choice(p2.mkt_grimm_c.values, len(p2)).mean() - rng.choice(p1.mkt_grimm_c.values, len(p1)).mean()
          for _ in range(5000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    # Theil-Sen slope per decade
    yrs, v = g.year.values, g.mkt_grimm_c.values
    sl = [(v[j] - v[i]) / (yrs[j] - yrs[i]) for i in range(len(v)) for j in range(i + 1, len(v)) if yrs[j] != yrs[i]]
    adj1, adj2 = p1.mkt_grimm_c - BIAS, p2.mkt_grimm_c - BIAS
    r = dict(city=city, country=g.country.iloc[0], zone=zone, test=T,
             source="ERA5" if g.is_era5.iloc[0] else "Station",
             n1=len(p1), n2=len(p2),
             p1=p1.mkt_grimm_c.mean(), p2=p2.mkt_grimm_c.mean(), d=d, lo=lo, hi=hi,
             slope=np.median(sl) * 10,
             p1c=adj1.mean(), p2c=adj2.mean(), margin_c=T - adj2.mean(),
             above1=int((adj1 > T).sum()), above2=int((adj2 > T).sum()))
    if g.is_era5.iloc[0]:
        r.update(pd1=p1.pd_hpa.mean(), pd2=p2.pd_hpa.mean(),
                 e50=(p2.mkt_ea50_c - BIAS).mean(), e120=(p2.mkt_ea120_c - BIAS).mean())
    def months(Tm, ea=83144):
        return 24 / np.exp(ea / R * (1 / (T + 273.15) - 1 / (Tm + 273.15)))
    r["life24"] = months(r["p2c"]) if r["p2c"] > T else None
    r["life24_50"] = months(r["p2c"], 50000) if r["p2c"] > T else None
    r["life24_120"] = months(r["p2c"], 120000) if r["p2c"] > T else None
    rows.append(r)
res = pd.DataFrame(rows).sort_values("margin_c")

# Validation: ERA5 minus station where both exist in the same year
val = []
for city, g in df.groupby("city"):
    e, s = g[g.is_era5], g[~g.is_era5]
    m = e.merge(s, on="year", suffixes=("_e", "_s"))
    if len(m):
        val.append(dict(city=city, n=len(m), years=f"{m.year.min()}–{m.year.max()}",
                        diff=(m.mkt_grimm_c_e - m.mkt_grimm_c_s).mean()))
val = pd.DataFrame(val)

out = dict(places=res.to_dict(orient="records"), validation=val.to_dict(orient="records"),
           n_places=len(res), n_era5=int((res.source == "ERA5").sum()), n_station=int((res.source == "Station").sum()),
           n_warmed=int((res.lo > 0).sum()), n_above=int((res.margin_c < 0).sum()),
           d_min=res.d.min(), d_max=res.d.max(), bias=BIAS)
json.dump(out, open("results.json", "w"), indent=1, default=float)

# ---------------------------------------------------------------- figures
plt.rcParams.update({"font.family": "Carlito", "font.size": 9, "axes.edgecolor": "#888",
                     "xtick.color": "#555", "ytick.color": "#555"})
INK, QUIET, GREY, ACCENT = "#1f1f1f", "#666666", "#b9b9b4", "#2563c9"

r2 = res.copy()
r2["x1"], r2["x2"] = r2.p1c - r2.test, r2.p2c - r2.test
r2 = r2.sort_values("x2", ascending=True)
h = 0.24 * len(r2) + 1.0
fig, ax = plt.subplots(figsize=(6.8, h), dpi=300)
for i, row in enumerate(r2.itertuples()):
    ax.plot([row.x1, row.x2], [i, i], color=GREY, lw=2, zorder=1)
    ax.plot(row.x1, i, "o", color=GREY, ms=5.5, zorder=2)
    ax.plot(row.x2, i, "o", color=ACCENT, ms=5.5, zorder=3)
    over = row.x2 > 0
    ax.text(r2.x2.max() + 0.45, i, f"{row.x2:+.2f}", va="center", fontsize=8.5,
            color=INK if over else QUIET, weight="bold" if over else "normal")
labels = [f"{c}{' (station)' if s == 'Station' else ''}" for c, s in zip(r2.city, r2.source)]
ax.set_yticks(range(len(r2))); ax.set_yticklabels(labels, fontsize=8.5)
for t, x in zip(ax.get_yticklabels(), r2.x2):
    if x > 0: t.set_fontweight("bold"); t.set_color(INK)
ax.axvline(0, color="#555", lw=1.2)
ax.text(0.08, len(r2) - 0.4, "test temperature", fontsize=8, color=QUIET, va="bottom")
ax.set_xlabel("Bias-corrected mean kinetic temperature minus long-term test temperature (°C)")
ax.set_ylim(-0.7, len(r2) - 0.1)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.grid(axis="x", color="#ececec", lw=0.6); ax.set_axisbelow(True)
ax.plot([], [], "o", color=GREY, label="1981–2000"); ax.plot([], [], "o", color=ACCENT, label="2006–2025")
ax.legend(loc="lower right", frameon=False, fontsize=8.5)
fig.savefig("figures/fig2_margins.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# humidity figure (ERA5 cities, Zone IV)
hz = res[(res.source == "ERA5") & (res.zone == "IV")].sort_values("pd2")
fig, ax = plt.subplots(figsize=(6.8, 0.24 * len(hz) + 1.0), dpi=300)
for i, row in enumerate(hz.itertuples()):
    ax.plot([row.pd1, row.pd2], [i, i], color=GREY, lw=2)
    ax.plot(row.pd1, i, "o", color=GREY, ms=5.5); ax.plot(row.pd2, i, "o", color=ACCENT, ms=5.5)
ax.set_yticks(range(len(hz))); ax.set_yticklabels(hz.city, fontsize=8.5)
for val_, lab in [(27.6, "30 °C / 65%"), (29.7, "30 °C / 70%"), (31.8, "30 °C / 75%")]:
    ax.axvline(val_, color="#777", lw=0.9, ls=(0, (4, 3)))
    ax.text(val_ + 0.08, len(hz) - 0.25, lab, fontsize=7.5, color=QUIET, ha="left", va="bottom")
ax.set_xlabel("Mean water-vapour partial pressure (hPa); dashed lines = test-chamber values")
ax.set_ylim(-0.7, len(hz) + 0.3)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.grid(axis="x", color="#ececec", lw=0.6); ax.set_axisbelow(True)
ax.plot([], [], "o", color=GREY, label="1981–2000"); ax.plot([], [], "o", color=ACCENT, label="2006–2025")
ax.legend(loc="lower left", bbox_to_anchor=(0.30, 0.02), frameon=False, fontsize=8.5)
fig.savefig("figures/fig3_humidity.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

pd.set_option("display.width", 220)
print(res[["city", "source", "test", "p1", "p2", "d", "lo", "hi", "p2c", "margin_c", "above1", "above2", "slope"]].round(2).to_string(index=False))
print(val.round(2).to_string(index=False))
print({k: v for k, v in out.items() if k not in ("places", "validation")})
