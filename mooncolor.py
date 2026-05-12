"""mooncolor — astronomical moonlight lamp.

A standalone HTTP service that drives a Philips Hue Color bulb to
match the *actual* moon at your location: real altitude, real phase,
real atmospheric color through the local aerosol optical depth (AOD),
clamped to your bulb's published color gamut.

What it does, plainly
---------------------
The moon's apparent color is not white. It depends on:
  - The sun's spectrum (moonlight is reflected sunlight)
  - The moon's spectral reflectance (slightly redder than the sun)
  - How much atmosphere the moonlight passes through to reach you
  - Aerosols / haze / smoke in that atmosphere (which redden it further)
  - Ozone in the upper atmosphere (which absorbs yellow at low altitude)
  - The lunar phase fraction (a thin crescent has a strong earthshine
    blue contribution from the un-lit limb)
  - Whether the moon is even above the horizon

This service computes all of that from first principles + live data,
projects the resulting xy chromaticity onto your bulb's published
gamut triangle, and PUTs `(xy, brightness)` to a Hue Bridge V1 API.
The lamp tracks the real moon in real time, all night, every night.

How you trigger it
------------------
HTTP webhook:

    POST http://<host>:<port>/lamp?state=on     # start tracking
    POST http://<host>:<port>/lamp?state=off    # stop + lamp off

Alternate enter/exit form (for callers that send state-machine events
rather than levels — e.g. a home-automation system delegating moonlight
mode to this service):

    POST http://<host>:<port>/event?event=enter
    POST http://<host>:<port>/event?event=exit

Either route is idempotent: re-asserting the current state is a no-op.

What's running while the lamp is "on"
-------------------------------------
A tick loop wakes every `tick.period_s` seconds (default 5 min). Each
tick:

  1. Compute moon altitude + illuminated fraction (astronomy-engine,
     pure-Python ephemeris — no API call needed).
  2. If altitude <= 0, turn the lamp off and skip.
  3. Read the latest atmospheric inputs (AOD, PM2.5, PM10, RH,
     pressure) from the in-memory poller cache.
  4. Run the physics pipeline → xy chromaticity → gamut-clamped xy.
  5. Compute brightness as a function of phase fraction (linear ramp
     between `min_brightness` at new moon and `max_brightness` at full).
  6. PUT `{"on": true, "xy": [...], "bri": N, "transitiontime": T}`
     to the Hue Bridge V1 API for the configured light_id.

Atmospheric pollers run continuously in the background, regardless of
whether the tick loop is active. This means the first tick after the
lamp is enabled uses already-warm data — no waiting for a fresh fetch.
If a source is unreachable, the last good value is used until it ages
out (`atmospheric.max_age_s`), then the configured default kicks in.

Physics (high level — derivations are inline below)
---------------------------------------------------
  Astronomy        astronomy-engine → moon altitude, phase fraction.
  Solar SPD        5778 K Planckian blackbody (≈ ASTM E490 AM0;
                   <1 % deviation across the visible band). Used
                   instead of D65 so we don't pre-bluify by including
                   Earth's atmospheric Rayleigh scattering — the moon
                   doesn't see that.
  Lunar reflectance  ROLO-derived B-V tilt, k=2.27 (reproduces the
                     observed lunar B-V index of 0.92 vs solar 0.65).
  Earthshine       At low phase fraction, the un-lit limb reflects
                   sunlight off Earth (oceans + clouds — bluer than
                   direct moonlight). Disk-integrated effective
                   reflectance: φ·refl_lunar + (1−φ)·α_E·refl_earth.
                   Meaningful only at φ < 0.15.
  Airmass          Kasten-Young 1989 (accurate down to low altitudes).
  Rayleigh         λ⁻⁴ vertical optical depth, pressure-scaled.
  Aerosol          AOD@550 from a CAMS-like source, Ångström-
                   extrapolated per-λ via α derived from PM10/PM2.5
                   mass ratio with hygroscopic-growth correction
                   (Hänel: high RH → swollen droplets → flatter
                   spectrum, models fog / marine layer correctly).
                   Multi-scatter forward-factor (1 − g·ω ≈ 0.36)
                   accounts for forward-scattered light reaching the
                   observer along the direct path. Optionally cross-
                   checks against a PurpleAir-derived AOD as a smoke-
                   event override when local sensors detect smoke the
                   upstream CAMS data hasn't ingested.
  Ozone            Chappuis-band absorption around 590-610 nm
                   (Bogumil 2003 cross-sections × 300 DU column ×
                   Kasten-Young airmass — removes yellow at low alt).
  Color rendering  spectrum_at_observer = solar_AM0 × lunar_refl ×
                   atm_trans, integrated against CIE 1931 2° standard
                   observer (x̄, ȳ, z̄) at 10 nm steps → X/Y/Z →
                   xy chromaticity → projected onto your bulb's
                   published gamut triangle if outside.

Resolution is limited by the granularity of public data sources
(CAMS air-quality is hourly; PurpleAir is ~10 min).

Configuration
-------------
YAML at `$MOONCOLOR_CONFIG` (default: `./config.yaml`). See
`config.example.yaml` for the schema and inline guidance.

Secrets via env:
  HUE_API_KEY        required — 40-char Hue Bridge V1 API key
  PURPLEAIR_API_KEY  optional — only if the PurpleAir cross-check is enabled

Limits
------
* Single bulb per service instance. To drive multiple bulbs from one
  instance, generalize `hue.light_id` to a list and dispatch the
  same (xy, bri) to all of them. If the bulbs are in different rooms
  with different gamuts, run multiple instances instead.
* Hue Bridge V1 API caps at ~10 commands/sec. We emit one per tick
  (default 5 min) — vastly under-limit.
* Hue V1 is plain HTTP; LAN-only.

License: see LICENSE.
"""

from __future__ import annotations

import argparse
import asyncio
import contextlib
import json
import logging
import math
import os
import signal
import sys
import time
from dataclasses import dataclass, field
from typing import Any

import aiohttp
import yaml
from aiohttp import web
from astronomy import (
    Body,
    Equator,
    Horizon,
    Illumination,
    Observer,
    Refraction,
    Time,
)


# ====================================================================
# PHYSICS.
# Derivations are inline with each function below. Do not modify any
# of these constants or routines without a corresponding physics
# review — small numerical changes here have outsized effects on the
# rendered color, and the values below are calibrated against
# established literature references (cited at each use site).
# ====================================================================

# Spectral integration grid: 380–780 nm, 10 nm steps (41 points).
LAMBDAS_NM = list(range(380, 781, 10))
LAMBDAS_UM = [l / 1000.0 for l in LAMBDAS_NM]

# CIE 1931 2° standard observer color-matching functions
# (x̄, ȳ, z̄), 10 nm steps 380–780 nm. Source: CIE 015:2018 Table T.4.
CIE_X_BAR = [
    0.001368, 0.004243, 0.013438, 0.034570, 0.077630, 0.137200, 0.230420, 0.302730,
    0.343230, 0.336200, 0.290800, 0.195360, 0.095640, 0.032010, 0.004900, 0.009300,
    0.063270, 0.165500, 0.290400, 0.433450, 0.594500, 0.762100, 0.916300, 1.026300,
    1.062200, 1.002600, 0.854450, 0.642400, 0.447900, 0.283500, 0.164900, 0.087400,
    0.046770, 0.022700, 0.011359, 0.005790, 0.002899, 0.001440, 0.000690, 0.000332, 0.000166,
]
CIE_Y_BAR = [
    0.000039, 0.000120, 0.000396, 0.001210, 0.002900, 0.006100, 0.010800, 0.018800,
    0.030100, 0.048000, 0.073900, 0.116200, 0.180800, 0.287000, 0.426400, 0.594800,
    0.762100, 0.875200, 0.954000, 0.994950, 0.995000, 0.952000, 0.870000, 0.757300,
    0.631000, 0.503000, 0.381000, 0.265000, 0.175000, 0.107000, 0.061000, 0.032000,
    0.017000, 0.008210, 0.004102, 0.002091, 0.001047, 0.000520, 0.000249, 0.000120, 0.000060,
]
CIE_Z_BAR = [
    0.006450, 0.020050, 0.064500, 0.168490, 0.385400, 0.656760, 1.082900, 1.391800,
    1.622960, 1.747060, 1.669200, 1.287640, 0.812950, 0.465180, 0.272000, 0.158200,
    0.078250, 0.042160, 0.020300, 0.008750, 0.003900, 0.002100, 0.001650, 0.001100,
    0.000800, 0.000340, 0.000190, 0.000050, 0.000020, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
]


def _planck_spd(lambdas_nm: list[int], temperature_k: float) -> list[float]:
    h = 6.62607015e-34
    c = 2.99792458e8
    k = 1.380649e-23
    raw = []
    for lam_nm in lambdas_nm:
        lam = lam_nm * 1e-9
        b = (2.0 * h * c * c / lam ** 5) / (math.exp(h * c / (lam * k * temperature_k)) - 1.0)
        raw.append(b)
    peak = max(raw)
    return [r / peak * 100.0 for r in raw]


SOLAR_AM0_SPD = _planck_spd(LAMBDAS_NM, 5778)
AOD_REF_UM = 0.550
OZONE_DU_DEFAULT = 300
DU_TO_MOLECULES_PER_CM2 = 2.687e16

# Ozone Chappuis-band absorption cross-sections at 380–780 nm in 10 nm
# steps. Units: 10⁻²¹ cm²/molecule. Source: Bogumil 2003.
OZONE_XSECT_E21_CM2 = [
    0.002, 0.003, 0.005, 0.008, 0.011, 0.014, 0.017, 0.022, 0.030, 0.040,
    0.055, 0.075, 0.105, 0.165, 0.255, 0.385, 0.555, 0.815, 1.180, 1.640,
    2.150, 2.650, 3.000, 3.080, 2.900, 2.500, 2.000, 1.480, 1.030, 0.700,
    0.480, 0.320, 0.225, 0.165, 0.120, 0.085, 0.060, 0.045, 0.035, 0.025, 0.020,
]

AEROSOL_ASYMMETRY = 0.7
AEROSOL_SSA = 0.92
AEROSOL_FORWARD_FACTOR = 1.0 - AEROSOL_ASYMMETRY * AEROSOL_SSA  # ≈ 0.356
HANEL_GROWTH_EXPONENT = 0.25
EARTHSHINE_INTENSITY_FACTOR = 0.02
K_LUNAR_REFL = 2.27
K_EARTH_REFL = -5.27


def _airmass_kasten_young(altitude_deg: float) -> float:
    a = max(altitude_deg, 0.1)
    return 1.0 / (
        math.sin(math.radians(a)) + 0.50572 * (a + 6.07995) ** -1.6364
    )


def _angstrom_alpha_humidity_corrected(pm10: float, pm25: float, rh: float) -> float:
    if pm25 <= 0:
        alpha_dry = 1.0
    else:
        ratio = max(pm10, pm25) / pm25
        alpha_dry = max(0.3, min(1.5, 1.5 / ratio))
    rh_clamped = max(0.0, min(99.0, rh))
    g_rh = (1.0 - rh_clamped / 100.0) ** -HANEL_GROWTH_EXPONENT
    return max(0.0, alpha_dry / (g_rh ** 1.5))


def _atmospheric_transmittance(
    alt_deg: float,
    aod_550: float,
    alpha: float,
    pressure_inhg: float,
    ozone_du: float = OZONE_DU_DEFAULT,
) -> list[float]:
    m = _airmass_kasten_young(alt_deg)
    p_factor = pressure_inhg / 29.92
    n_o3 = ozone_du * DU_TO_MOLECULES_PER_CM2
    out = []
    for i, l in enumerate(LAMBDAS_UM):
        tau_r = 0.008569 * p_factor * (l ** -4)
        tau_a_geom = aod_550 * (l / AOD_REF_UM) ** -alpha
        tau_a = AEROSOL_FORWARD_FACTOR * tau_a_geom
        tau_o = (OZONE_XSECT_E21_CM2[i] * 1e-21) * n_o3
        out.append(math.exp(-(tau_r + tau_a + tau_o) * m))
    return out


def _refl_total(phi: float) -> list[float]:
    refl_lunar = [1.0 + K_LUNAR_REFL * math.log10(l / 0.550) for l in LAMBDAS_UM]
    refl_earth = [
        max(0.001, 1.0 + K_EARTH_REFL * math.log10(l / 0.550)) for l in LAMBDAS_UM
    ]
    return [
        phi * rl + (1 - phi) * EARTHSHINE_INTENSITY_FACTOR * re
        for rl, re in zip(refl_lunar, refl_earth)
    ]


def _moon_xy(
    alt_deg: float,
    aod_550: float,
    alpha: float,
    pressure_inhg: float,
    phi: float,
    ozone_du: float = OZONE_DU_DEFAULT,
) -> list[float] | None:
    if alt_deg <= -0.8:
        return None
    trans = _atmospheric_transmittance(
        alt_deg, aod_550, alpha, pressure_inhg, ozone_du
    )
    refl = _refl_total(phi)
    X = Y = Z = 0.0
    for i in range(len(LAMBDAS_NM)):
        s = SOLAR_AM0_SPD[i] * refl[i] * trans[i]
        X += s * CIE_X_BAR[i]
        Y += s * CIE_Y_BAR[i]
        Z += s * CIE_Z_BAR[i]
    total = X + Y + Z
    if total <= 0:
        return None
    return [X / total, Y / total]


def _point_in_triangle(p: list[float], a: list[float], b: list[float], c: list[float]) -> bool:
    def sign(p1: list[float], p2: list[float], p3: list[float]) -> float:
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    d1, d2, d3 = sign(p, a, b), sign(p, b, c), sign(p, c, a)
    return not (
        ((d1 < 0) or (d2 < 0) or (d3 < 0))
        and ((d1 > 0) or (d2 > 0) or (d3 > 0))
    )


def _project_onto_segment(
    p: list[float], a: list[float], b: list[float]
) -> list[float]:
    apx, apy = p[0] - a[0], p[1] - a[1]
    abx, aby = b[0] - a[0], b[1] - a[1]
    denom = abx * abx + aby * aby
    if denom == 0:
        return list(a)
    t = max(0.0, min(1.0, (apx * abx + apy * aby) / denom))
    return [a[0] + abx * t, a[1] + aby * t]


def _clamp_to_gamut(
    xy: list[float], vertices: dict[str, list[float]]
) -> list[float]:
    R, G, B = vertices["R"], vertices["G"], vertices["B"]
    if _point_in_triangle(xy, R, G, B):
        return xy
    candidates = [
        _project_onto_segment(xy, R, G),
        _project_onto_segment(xy, G, B),
        _project_onto_segment(xy, B, R),
    ]
    best = min(
        candidates, key=lambda p: (p[0] - xy[0]) ** 2 + (p[1] - xy[1]) ** 2
    )
    return [round(best[0], 4), round(best[1], 4)]


def compute_moon_color_xy(
    alt_deg: float,
    aod_550: float,
    alpha: float,
    pressure_inhg: float,
    phi: float,
    bulb: dict[str, Any],
    ozone_du: float = OZONE_DU_DEFAULT,
) -> list[float] | None:
    xy = _moon_xy(alt_deg, aod_550, alpha, pressure_inhg, phi, ozone_du)
    if xy is None:
        return None
    return _clamp_to_gamut([round(xy[0], 4), round(xy[1], 4)], bulb["gamut_vertices"])


# ====================================================================
# CONFIG
# ====================================================================


@dataclass(frozen=True)
class AtmosphericConfig:
    """Defaults + freshness limits for each input source."""

    aod_default: float
    pm25_default: float
    pm10_default: float
    pressure_inhg_default: float
    rh_default: float
    max_age_s: int
    open_meteo_poll_s: int
    purpleair_poll_s: int
    purpleair_url: str | None
    purpleair_api_key: str | None
    purpleair_use_epa_correction: bool
    purpleair_pm25_to_aod_k: float
    purpleair_trigger_ratio: float


@dataclass(frozen=True)
class HueConfig:
    bridge_ip: str
    api_key: str
    light_id: int
    http_timeout_s: float


@dataclass(frozen=True)
class TickConfig:
    period_s: int
    transition_s: int
    max_brightness: int
    min_brightness: int


@dataclass(frozen=True)
class Config:
    listen_host: str
    listen_port: int
    latitude: float
    longitude: float
    elevation_m: float
    bulb_profile_path: str
    atmospheric: AtmosphericConfig
    hue: HueConfig
    tick: TickConfig
    log_level: str


def _require_env(name: str, allow_empty: bool = False) -> str:
    val = os.environ.get(name, "").strip()
    if not val and not allow_empty:
        logging.critical("missing required env var %s", name)
        sys.exit(1)
    return val


def load_config(path: str) -> Config:
    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    listen = raw.get("listen") or {}
    loc = raw.get("location") or {}
    atm = raw.get("atmospheric") or {}
    hue = raw.get("hue") or {}
    tick = raw.get("tick") or {}

    if "latitude" not in loc or "longitude" not in loc:
        logging.critical("config: location.latitude and location.longitude required")
        sys.exit(1)

    bulb_profile = raw.get("bulb_profile")
    if not bulb_profile:
        logging.critical("config: bulb_profile required")
        sys.exit(1)
    # Resolve relative to config file directory.
    config_dir = os.path.dirname(os.path.abspath(path))
    bulb_profile_path = (
        bulb_profile
        if os.path.isabs(bulb_profile)
        else os.path.join(config_dir, bulb_profile)
    )

    hue_key = _require_env(hue.get("api_key_env", "HUE_API_KEY"))
    bridge_ip = (hue.get("bridge_ip") or "").strip()
    if not bridge_ip:
        logging.critical("config: hue.bridge_ip required")
        sys.exit(1)
    light_id = hue.get("light_id")
    if light_id is None:
        logging.critical("config: hue.light_id required")
        sys.exit(1)

    pa_url = atm.get("purpleair_url")
    pa_key = None
    if pa_url:
        pa_key = _require_env(
            atm.get("purpleair_api_key_env", "PURPLEAIR_API_KEY"), allow_empty=False
        )

    return Config(
        listen_host=listen.get("host", "0.0.0.0"),
        listen_port=int(listen.get("port", 8081)),
        latitude=float(loc["latitude"]),
        longitude=float(loc["longitude"]),
        elevation_m=float(loc.get("elevation_m", 0.0)),
        bulb_profile_path=bulb_profile_path,
        atmospheric=AtmosphericConfig(
            aod_default=float(atm.get("aod_default", 0.10)),
            pm25_default=float(atm.get("pm25_default", 8.0)),
            pm10_default=float(atm.get("pm10_default", 12.0)),
            pressure_inhg_default=float(atm.get("pressure_inhg_default", 30.0)),
            rh_default=float(atm.get("rh_default", 50.0)),
            max_age_s=int(atm.get("max_age_s", 3 * 3600)),
            open_meteo_poll_s=int(atm.get("open_meteo_poll_s", 1800)),
            purpleair_poll_s=int(atm.get("purpleair_poll_s", 1800)),
            purpleair_url=pa_url,
            purpleair_api_key=pa_key,
            purpleair_use_epa_correction=bool(
                atm.get("purpleair_use_epa_correction", True)
            ),
            purpleair_pm25_to_aod_k=float(atm.get("purpleair_pm25_to_aod_k", 0.003)),
            purpleair_trigger_ratio=float(atm.get("purpleair_trigger_ratio", 1.5)),
        ),
        hue=HueConfig(
            bridge_ip=bridge_ip,
            api_key=hue_key,
            light_id=int(light_id),
            http_timeout_s=float(hue.get("http_timeout_s", 5.0)),
        ),
        tick=TickConfig(
            period_s=int(tick.get("period_s", 300)),
            transition_s=int(tick.get("transition_s", 30)),
            max_brightness=int(tick.get("max_brightness", 200)),
            min_brightness=int(tick.get("min_brightness", 3)),
        ),
        log_level=str(raw.get("log_level", "INFO")).upper(),
    )


# ====================================================================
# ATMOSPHERIC DATA — direct REST polling.
# Open-Meteo Air-Quality: AOD, PM2.5, PM10. No key needed.
# Open-Meteo Forecast: pressure, RH. No key needed.
# PurpleAir: local PM2.5 cross-check. Free API, key required.
# ====================================================================


@dataclass
class AtmosphericReading:
    """A timestamped sample of one input."""

    value: float | None = None
    updated_at: float = 0.0


@dataclass
class AtmosphericState:
    """In-memory cache of latest readings, updated by background pollers."""

    aod_550: AtmosphericReading = field(default_factory=AtmosphericReading)
    pm25_om: AtmosphericReading = field(default_factory=AtmosphericReading)
    pm10: AtmosphericReading = field(default_factory=AtmosphericReading)
    pressure_inhg: AtmosphericReading = field(default_factory=AtmosphericReading)
    rh: AtmosphericReading = field(default_factory=AtmosphericReading)
    pm25_purpleair: AtmosphericReading = field(default_factory=AtmosphericReading)


def _value_or_default(
    reading: AtmosphericReading, default: float, max_age_s: int, label: str
) -> float:
    """Return the cached value if fresh; else the default. Logs once on
    staleness transition would be nicer but state-tracking adds complexity
    we don't need yet — we just log every time we fall back."""
    if reading.value is None:
        logging.warning(
            "atmospheric: no reading for %s — using default %s", label, default
        )
        return default
    age = time.time() - reading.updated_at
    if age > max_age_s:
        logging.warning(
            "atmospheric: %s stale (%.0fs > %ds) — using default %s",
            label,
            age,
            max_age_s,
            default,
        )
        return default
    return reading.value


async def _poll_open_meteo_air_quality(
    session: aiohttp.ClientSession,
    config: Config,
    state: AtmosphericState,
) -> None:
    """One scrape of Open-Meteo Air-Quality API. Updates aod/pm25/pm10."""
    url = (
        "https://air-quality-api.open-meteo.com/v1/air-quality"
        f"?latitude={config.latitude}&longitude={config.longitude}"
        "&current=pm2_5,pm10,aerosol_optical_depth"
    )
    try:
        async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as resp:
            resp.raise_for_status()
            body = await resp.json()
    except Exception as e:
        logging.warning("open-meteo air-quality poll failed: %s", e)
        return
    cur = body.get("current") or {}
    now = time.time()
    aod = cur.get("aerosol_optical_depth")
    pm25 = cur.get("pm2_5")
    pm10 = cur.get("pm10")
    if aod is not None:
        state.aod_550 = AtmosphericReading(value=float(aod), updated_at=now)
    if pm25 is not None:
        state.pm25_om = AtmosphericReading(value=float(pm25), updated_at=now)
    if pm10 is not None:
        state.pm10 = AtmosphericReading(value=float(pm10), updated_at=now)
    logging.debug(
        "atmospheric/open-meteo-aq: aod=%s pm25=%s pm10=%s", aod, pm25, pm10
    )


async def _poll_open_meteo_forecast(
    session: aiohttp.ClientSession,
    config: Config,
    state: AtmosphericState,
) -> None:
    """One scrape of Open-Meteo Forecast API. Updates pressure + RH.

    Pressure unit returned by Open-Meteo (with no `pressure_unit=` query
    param) is hPa; we convert to inHg here so the cached value matches
    the physics module's expected unit.
    """
    url = (
        "https://api.open-meteo.com/v1/forecast"
        f"?latitude={config.latitude}&longitude={config.longitude}"
        "&current=pressure_msl,relative_humidity_2m"
    )
    try:
        async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as resp:
            resp.raise_for_status()
            body = await resp.json()
    except Exception as e:
        logging.warning("open-meteo forecast poll failed: %s", e)
        return
    cur = body.get("current") or {}
    now = time.time()
    p_hpa = cur.get("pressure_msl")
    rh = cur.get("relative_humidity_2m")
    if p_hpa is not None:
        state.pressure_inhg = AtmosphericReading(
            value=float(p_hpa) * 0.02953, updated_at=now
        )
    if rh is not None:
        state.rh = AtmosphericReading(value=float(rh), updated_at=now)
    logging.debug(
        "atmospheric/open-meteo-forecast: pressure=%s hPa rh=%s%%", p_hpa, rh
    )


async def _poll_purpleair(
    session: aiohttp.ClientSession,
    config: Config,
    state: AtmosphericState,
) -> None:
    """One scrape of the PurpleAir bounding-box API. Updates local PM2.5
    mean (used as a smoke override against CAMS-derived AOD)."""
    if (
        not config.atmospheric.purpleair_url
        or not config.atmospheric.purpleair_api_key
    ):
        return
    headers = {"X-API-Key": config.atmospheric.purpleair_api_key}
    try:
        async with session.get(
            config.atmospheric.purpleair_url,
            headers=headers,
            timeout=aiohttp.ClientTimeout(total=10),
        ) as resp:
            resp.raise_for_status()
            body = await resp.json()
    except Exception as e:
        logging.warning("purpleair poll failed: %s", e)
        return
    rows = body.get("data") or []
    total = 0.0
    count = 0
    for row in rows:
        if not isinstance(row, list) or len(row) < 2:
            continue
        try:
            v = float(row[1])
        except (TypeError, ValueError):
            continue
        if v > 0:
            total += v
            count += 1
    if count == 0:
        logging.warning("purpleair: no sensors with positive readings — keeping last value")
        return
    mean = total / count
    state.pm25_purpleair = AtmosphericReading(value=mean, updated_at=time.time())
    logging.debug(
        "atmospheric/purpleair: mean pm25=%.1f across %d sensors", mean, count
    )


async def _atmospheric_poller(
    app: web.Application, session: aiohttp.ClientSession
) -> None:
    """Background coroutine: polls each source on its own cadence
    until the app is shutting down."""
    config: Config = app["config"]
    state: AtmosphericState = app["atmospheric_state"]
    stop = app["stop_event"]

    # Initial fan-out (don't wait a full cycle to populate).
    await asyncio.gather(
        _poll_open_meteo_air_quality(session, config, state),
        _poll_open_meteo_forecast(session, config, state),
        _poll_purpleair(session, config, state),
        return_exceptions=True,
    )

    async def _loop_one(
        poll: Any, interval_s: int, label: str
    ) -> None:
        while not stop.is_set():
            try:
                await asyncio.wait_for(stop.wait(), timeout=interval_s)
                return  # stop set during sleep
            except asyncio.TimeoutError:
                pass
            try:
                await poll(session, config, state)
            except Exception:
                logging.exception("atmospheric %s poll error", label)

    tasks = [
        asyncio.create_task(
            _loop_one(
                _poll_open_meteo_air_quality,
                config.atmospheric.open_meteo_poll_s,
                "open-meteo-aq",
            )
        ),
        asyncio.create_task(
            _loop_one(
                _poll_open_meteo_forecast,
                config.atmospheric.open_meteo_poll_s,
                "open-meteo-forecast",
            )
        ),
        asyncio.create_task(
            _loop_one(
                _poll_purpleair,
                config.atmospheric.purpleair_poll_s,
                "purpleair",
            )
        ),
    ]
    try:
        await asyncio.gather(*tasks, return_exceptions=True)
    finally:
        for t in tasks:
            if not t.done():
                t.cancel()


def _effective_aod(
    config: Config, state: AtmosphericState
) -> float:
    """Resolve the effective AOD to use for this tick.

    Default is the CAMS-derived AOD from the air-quality API. If the
    PurpleAir cross-check is enabled AND a fresh PurpleAir reading is
    available AND the PurpleAir-derived AOD exceeds the CAMS value by
    `purpleair_trigger_ratio`, override with the PurpleAir value.

    Why: CAMS aerosol forecasts are produced upstream from satellite
    + chemistry-transport modeling and can lag local smoke events
    (wildfires, controlled burns) by hours. Local PurpleAir sensors
    detect smoke immediately. We want the lamp to redden when the
    actual atmosphere is smoky, not when the upstream model says it
    should be.

    The PurpleAir → AOD conversion uses the Barkjohn 2021 / EPA
    correction to translate raw PM2.5 to corrected PM2.5, then
    multiplies by `purpleair_pm25_to_aod_k` (a mass-extinction ×
    boundary-layer-height product). The trigger ratio prevents
    flapping when the two sources are close — only override when
    PurpleAir clearly disagrees in the smoky direction."""
    a = config.atmospheric
    aod_cams = _value_or_default(
        state.aod_550, a.aod_default, a.max_age_s, "aod_550"
    )
    if not a.purpleair_url:
        return aod_cams
    pa_raw_reading = state.pm25_purpleair
    if (
        pa_raw_reading.value is None
        or (time.time() - pa_raw_reading.updated_at) > a.max_age_s
    ):
        return aod_cams
    pa_raw = pa_raw_reading.value
    if a.purpleair_use_epa_correction:
        rh = _value_or_default(state.rh, a.rh_default, a.max_age_s, "rh")
        pa_corrected = max(0.0, 0.524 * pa_raw - 0.086 * rh + 5.75)
    else:
        pa_corrected = pa_raw
    aod_local = a.purpleair_pm25_to_aod_k * pa_corrected
    if aod_local > a.purpleair_trigger_ratio * aod_cams:
        logging.info(
            "smoke override — PurpleAir raw=%.1f corrected=%.1f → aod=%.3f "
            "> %sx CAMS (%.3f); using PurpleAir",
            pa_raw,
            pa_corrected,
            aod_local,
            a.purpleair_trigger_ratio,
            aod_cams,
        )
        return aod_local
    return aod_cams


# ====================================================================
# HUE V1 API CLIENT
# ====================================================================


class HueClient:
    """Minimal Hue V1 light-set client. Two operations: set color+brightness,
    turn off."""

    def __init__(self, bridge_ip: str, api_key: str, timeout_s: float):
        self._base = f"http://{bridge_ip}/api/{api_key}"
        self._timeout = aiohttp.ClientTimeout(total=timeout_s)
        self._session: aiohttp.ClientSession | None = None

    async def __aenter__(self) -> HueClient:
        self._session = aiohttp.ClientSession(timeout=self._timeout)
        return self

    async def __aexit__(self, *exc: Any) -> None:
        if self._session is not None:
            await self._session.close()
            self._session = None

    async def set_light(
        self,
        light_id: int,
        xy: list[float],
        brightness: int,
        transition_s: int,
    ) -> None:
        body = {
            "on": True,
            "xy": xy,
            "bri": brightness,
            # Hue's transitiontime is in units of 100 ms.
            "transitiontime": max(0, int(transition_s) * 10),
        }
        await self._put_light_state(light_id, body)

    async def turn_off(self, light_id: int, transition_s: int) -> None:
        body = {"on": False, "transitiontime": max(0, int(transition_s) * 10)}
        await self._put_light_state(light_id, body)

    async def _put_light_state(
        self, light_id: int, body: dict[str, Any]
    ) -> None:
        assert self._session is not None
        url = f"{self._base}/lights/{light_id}/state"
        # Path-only label for error messages — the full URL embeds the
        # API key, and these messages flow through journald → Vector →
        # VictoriaLogs, so the key must never appear there.
        log_path = f"/lights/{light_id}/state"
        async with self._session.put(url, json=body) as resp:
            text = await resp.text()
            if resp.status >= 400:
                raise HueError(
                    f"PUT {log_path} body={body!r} -> HTTP {resp.status} {text!r}"
                )
            try:
                parsed = json.loads(text)
            except json.JSONDecodeError:
                raise HueError(f"PUT {log_path} non-JSON: {text!r}")
            errors = [
                item["error"]
                for item in parsed
                if isinstance(item, dict) and "error" in item
            ]
            if errors:
                raise HueError(f"PUT {log_path} body={body!r} errors={errors!r}")


class HueError(RuntimeError):
    pass


# ====================================================================
# ASTRONOMY HELPER
# ====================================================================


def make_observer(config: Config) -> Observer:
    return Observer(config.latitude, config.longitude, config.elevation_m)


def moon_state_now(observer: Observer) -> tuple[float, float]:
    """(altitude_deg, illumination_fraction)."""
    t = Time.Now()
    equ = Equator(Body.Moon, t, observer, ofdate=True, aberration=True)
    hor = Horizon(t, observer, equ.ra, equ.dec, refraction=Refraction.Normal)
    illum = Illumination(Body.Moon, t).phase_fraction
    return hor.altitude, illum


# ====================================================================
# TICK LOOP — runs when state=on.
# ====================================================================


async def _tick_loop(app: web.Application) -> None:
    """Compute & push one moonlight frame every tick.period_s, until the
    stop signal is set (typically by a state=off webhook)."""
    config: Config = app["config"]
    state: AtmosphericState = app["atmospheric_state"]
    bulb: dict[str, Any] = app["bulb"]
    hue: HueClient = app["hue"]
    observer: Observer = app["observer"]
    tick_stop: asyncio.Event = app["tick_stop"]

    logging.info(
        "tick loop started — period=%ds transition=%ds",
        config.tick.period_s,
        config.tick.transition_s,
    )

    while not tick_stop.is_set():
        await _do_tick(app, config, state, bulb, hue, observer)
        try:
            await asyncio.wait_for(
                tick_stop.wait(), timeout=config.tick.period_s
            )
        except asyncio.TimeoutError:
            pass

    # Stop signal: turn the lamp off and exit cleanly.
    try:
        await hue.turn_off(config.hue.light_id, config.tick.transition_s)
        logging.info("tick loop stopped — lamp off")
    except Exception:
        logging.exception("tick loop: failed to turn off lamp on shutdown")


async def _do_tick(
    app: web.Application,
    config: Config,
    state: AtmosphericState,
    bulb: dict[str, Any],
    hue: HueClient,
    observer: Observer,
) -> None:
    try:
        alt, illum = moon_state_now(observer)
    except Exception:
        logging.exception("astronomy failed; skipping tick")
        return

    if alt <= 0:
        try:
            await hue.turn_off(config.hue.light_id, config.tick.transition_s)
            logging.info(
                "alt=%.1f° below horizon — lamp off (illum=%.2f)", alt, illum
            )
        except Exception:
            logging.exception("hue off (below-horizon) failed")
            app["state"]["dispatch_errors"] += 1
        app["state"]["last_tick"] = {
            "ts": time.time(),
            "alt": alt,
            "illum": illum,
            "lamp": "off (below horizon)",
        }
        return

    a = config.atmospheric
    aod_eff = _effective_aod(config, state)
    pm25 = _value_or_default(state.pm25_om, a.pm25_default, a.max_age_s, "pm25")
    pm10 = _value_or_default(state.pm10, a.pm10_default, a.max_age_s, "pm10")
    rh = _value_or_default(state.rh, a.rh_default, a.max_age_s, "rh")
    pressure = _value_or_default(
        state.pressure_inhg,
        a.pressure_inhg_default,
        a.max_age_s,
        "pressure",
    )
    alpha = _angstrom_alpha_humidity_corrected(pm10, pm25, rh)

    xy = compute_moon_color_xy(alt, aod_eff, alpha, pressure, illum, bulb)
    if xy is None:
        try:
            await hue.turn_off(config.hue.light_id, config.tick.transition_s)
        except Exception:
            logging.exception("hue off (no-xy) failed")
            app["state"]["dispatch_errors"] += 1
        app["state"]["last_tick"] = {
            "ts": time.time(),
            "alt": alt,
            "illum": illum,
            "lamp": "off (xy unresolvable)",
        }
        return

    bri = int(
        config.tick.min_brightness
        + illum * (config.tick.max_brightness - config.tick.min_brightness)
    )
    try:
        await hue.set_light(
            config.hue.light_id, xy, bri, config.tick.transition_s
        )
        app["state"]["dispatches"] += 1
        logging.info(
            "tick alt=%.1f° illum=%.2f aod=%.3f α=%.2f rh=%.0f%% "
            "P=%.2finHg → xy=%s bri=%d",
            alt,
            illum,
            aod_eff,
            alpha,
            rh,
            pressure,
            xy,
            bri,
        )
    except Exception:
        logging.exception("hue set_light failed")
        app["state"]["dispatch_errors"] += 1
    app["state"]["last_tick"] = {
        "ts": time.time(),
        "alt": alt,
        "illum": illum,
        "aod": aod_eff,
        "alpha": alpha,
        "rh": rh,
        "pressure_inhg": pressure,
        "xy": xy,
        "bri": bri,
        "lamp": "on",
    }


# ====================================================================
# HTTP HANDLERS
# ====================================================================


def _tick_task_done(task: asyncio.Task[None]) -> None:
    """Surface unhandled exceptions from the tick loop. asyncio
    otherwise silently swallows them until the task is garbage-collected
    (at which point a 'Task exception was never retrieved' warning fires
    to stderr — much harder to debug than a logged error)."""
    if task.cancelled():
        return
    exc = task.exception()
    if exc is not None:
        logging.error("tick loop terminated with exception", exc_info=exc)


def _parse_state(qs: dict[str, str] | Any) -> bool | None:
    """Return True for on, False for off, None for ambiguous/missing.

    Accepts two input shapes so both webhook styles work:

      - ?state=on|off|true|false|yes|no|1|0   (level-style)
      - ?event=enter|exit                     (event-style, for callers
                                              that delegate state-machine
                                              transitions to this service)
    """
    event = (qs.get("event") or "").strip().lower()
    if event == "enter":
        return True
    if event == "exit":
        return False
    state = (qs.get("state") or "").strip().lower()
    if state in ("on", "true", "1", "yes"):
        return True
    if state in ("off", "false", "0", "no"):
        return False
    return None


async def _handle_lamp(request: web.Request) -> web.Response:
    desired = _parse_state(request.query)
    if desired is None:
        return web.json_response(
            {"error": "missing or invalid state — pass ?state=on|off or ?event=enter|exit"},
            status=400,
        )

    app = request.app
    tick_task = app["state"]["tick_task"]
    currently_on = tick_task is not None and not tick_task.done()

    if desired and not currently_on:
        # Transition off→on. Start the tick loop. add_done_callback
        # surfaces unhandled exceptions (otherwise asyncio drops them
        # to stderr on garbage collection without context).
        app["tick_stop"].clear()
        task = asyncio.create_task(_tick_loop(app))
        task.add_done_callback(_tick_task_done)
        app["state"]["tick_task"] = task
        app["state"]["transitions"] += 1
        logging.info("trigger: state=on — tick loop started")
        return web.json_response({"state": "on", "transition": "started"})

    if not desired and currently_on:
        # Transition on→off. Signal the tick loop to stop; it turns off
        # the lamp itself.
        app["tick_stop"].set()
        app["state"]["transitions"] += 1
        logging.info("trigger: state=off — tick loop stop signal sent")
        return web.json_response({"state": "off", "transition": "stopping"})

    # Idempotent no-op.
    logging.debug(
        "trigger: state=%s already in desired state; no-op",
        "on" if desired else "off",
    )
    return web.json_response(
        {"state": "on" if desired else "off", "transition": "noop"}
    )


async def _handle_healthz(request: web.Request) -> web.Response:
    app = request.app
    state: AtmosphericState = app["atmospheric_state"]
    now = time.time()

    def freshness(r: AtmosphericReading) -> dict[str, Any]:
        return {
            "value": r.value,
            "age_s": (now - r.updated_at) if r.updated_at else None,
        }

    s = app["state"]
    tick_task = s["tick_task"]
    return web.json_response(
        {
            "status": "ok",
            "uptime_s": int(now - app["started_at"]),
            "tick": {
                "running": tick_task is not None and not tick_task.done(),
                "last": s["last_tick"],
            },
            "transitions": s["transitions"],
            "dispatches": s["dispatches"],
            "dispatch_errors": s["dispatch_errors"],
            "atmospheric": {
                "aod_550": freshness(state.aod_550),
                "pm25_om": freshness(state.pm25_om),
                "pm10": freshness(state.pm10),
                "rh": freshness(state.rh),
                "pressure_inhg": freshness(state.pressure_inhg),
                "pm25_purpleair": freshness(state.pm25_purpleair),
            },
        }
    )


# ====================================================================
# APP WIRING
# ====================================================================


def make_app(config: Config) -> web.Application:
    app = web.Application()
    app["config"] = config
    app["started_at"] = time.time()
    app["atmospheric_state"] = AtmosphericState()
    app["tick_stop"] = asyncio.Event()
    # aiohttp 3.10+ deprecates re-assigning app[key] after startup, but
    # mutating the contents of a value set at startup is fine. Stash
    # the counters + last_tick + the live tick_task ref inside a single
    # mutable dict so the dict reference stays constant in app[].
    app["state"] = {
        "tick_task": None,  # type: asyncio.Task[None] | None
        "transitions": 0,
        "dispatches": 0,
        "dispatch_errors": 0,
        "last_tick": None,
    }
    app["stop_event"] = asyncio.Event()  # shutdown signal for the atmospheric poller

    # Two routes accept identical query-shape semantics — `/lamp` is
    # the simple/canonical name, `/event` is for upstream callers that
    # think in state-machine terms (enter/exit). Both call the same
    # handler; _parse_state recognizes both `state=` and `event=`.
    app.router.add_post("/lamp", _handle_lamp)
    app.router.add_post("/event", _handle_lamp)
    app.router.add_get("/healthz", _handle_healthz)

    async def _on_startup(app: web.Application) -> None:
        # Bulb profile.
        with open(config.bulb_profile_path) as f:
            app["bulb"] = json.load(f)
        # Astronomy observer.
        app["observer"] = make_observer(config)
        # Hue client.
        client = HueClient(
            config.hue.bridge_ip, config.hue.api_key, config.hue.http_timeout_s
        )
        await client.__aenter__()
        app["hue"] = client
        # Atmospheric poller HTTP session (lifetime = app lifetime).
        app["atm_session"] = aiohttp.ClientSession()
        app["atm_task"] = asyncio.create_task(
            _atmospheric_poller(app, app["atm_session"])
        )
        logging.info(
            "mooncolor: ready listen=%s:%d bridge=%s light_id=%d bulb=%s",
            config.listen_host,
            config.listen_port,
            config.hue.bridge_ip,
            config.hue.light_id,
            app["bulb"].get("model", "?"),
        )

    async def _on_cleanup(app: web.Application) -> None:
        # Stop the tick loop if running. The loop turns the lamp off
        # itself on stop signal. If it doesn't drain cleanly within
        # the timeout, cancel it explicitly — wait_for() alone leaves
        # the task running, which would then hit a closed Hue session.
        tick_task = app["state"]["tick_task"]
        if tick_task is not None and not tick_task.done():
            app["tick_stop"].set()
            try:
                await asyncio.wait_for(tick_task, timeout=10)
            except asyncio.TimeoutError:
                logging.warning(
                    "tick loop did not stop within 10s — cancelling"
                )
                tick_task.cancel()
                with contextlib.suppress(Exception):
                    await tick_task
            except Exception:
                logging.exception("tick loop teardown")
        # Stop atmospheric pollers. Same wait-then-cancel pattern so a
        # poll mid-flight can't outlive the aiohttp session below.
        app["stop_event"].set()
        atm_task = app.get("atm_task")
        if atm_task is not None and not atm_task.done():
            try:
                await asyncio.wait_for(atm_task, timeout=5)
            except asyncio.TimeoutError:
                logging.warning(
                    "atmospheric poller did not stop within 5s — cancelling"
                )
                atm_task.cancel()
                with contextlib.suppress(Exception):
                    await atm_task
            except Exception:
                logging.exception("atmospheric poller teardown")
        atm_session = app.get("atm_session")
        if atm_session is not None:
            await atm_session.close()
        # Close Hue client.
        hue = app.get("hue")
        if hue is not None:
            await hue.__aexit__(None, None, None)
        logging.info("mooncolor: stopped")

    app.on_startup.append(_on_startup)
    app.on_cleanup.append(_on_cleanup)
    return app


# ====================================================================
# ENTRYPOINT
# ====================================================================


def _configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level, logging.INFO),
        format="%(asctime)s.%(msecs)03dZ %(levelname)s %(name)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        stream=sys.stdout,
    )
    logging.Formatter.converter = time.gmtime


def main() -> None:
    parser = argparse.ArgumentParser(description="Astronomical moonlight lamp service.")
    parser.add_argument(
        "--config",
        default=os.environ.get(
            "MOONCOLOR_CONFIG",
            "config.yaml",
        ),
    )
    args = parser.parse_args()

    config = load_config(args.config)
    _configure_logging(config.log_level)

    app = make_app(config)
    runner = web.AppRunner(app)

    async def _serve() -> None:
        await runner.setup()
        site = web.TCPSite(runner, config.listen_host, config.listen_port)
        await site.start()
        loop = asyncio.get_running_loop()
        stop = asyncio.Event()
        for sig in (signal.SIGTERM, signal.SIGINT):
            loop.add_signal_handler(sig, stop.set)
        await stop.wait()
        await runner.cleanup()

    with contextlib.suppress(KeyboardInterrupt):
        asyncio.run(_serve())


if __name__ == "__main__":
    main()
