# mooncolor

### An attempt at accurate mooncolor simulation for Hue bulbs.

> *Warning: This is 100% AI-coded between Claude Opus 4.7 and Gemini Pro.
> It is for a fun little thing so I did not sweat it. User beware.*

I 3D printed a moon and put a Hue bulb in it.

And it was good.

Time went by. And suddenly I had a thought... *That's no moon!* It doesn't match the moon at all.

So I fixed that by offloading all my critical thought into AI.

And it was good.

So here we are. A small program that should accurately model the color of the moon where you are, using a Hue bulb. The output is not a vibe — it's the spectrum that would enter your eye if the bulb were the moon: sunlight reflected off lunar regolith, attenuated by Rayleigh + aerosol + ozone, with earthshine blended in at low phase fractions, integrated against the CIE 1931 color-matching functions, and projected onto your bulb's color triangle.

I'm including a profile for the bulb I'm using (Philips Hue Color A19, LCA003). If you're using a different one, you may need to figure out your own gamut-vertex mapping — Philips publishes these for every Hue model.

The program is documented inline. There is no support here. Godspeed.

The next step is mechanical — to build a rotating shield over the bulb to physically render moon phases.

---

## Quick demo

  - Full moon, clear summer night → cool slightly-greenish white
  - Crescent moon → deep blue-tinted (earthshine dominates)
  - Wildfire smoke moves in → bulb shifts visibly red within ~30 min
  - Moon below horizon → bulb off
  - Service restart → resumes from cold within one poll cycle

Once it's running you stop thinking about it and just notice the room.

## What you need

  - Python 3.11+
  - A Philips Hue Color bulb (LCA003 / Hue Color A19, or anything
    with a published gamut triangle — drop the gamut JSON into
    `bulbs/` and point the config at it)
  - A Hue Bridge on your LAN
  - Your latitude, longitude, and approximate elevation
  - Optionally: a [PurpleAir](https://api.purpleair.com/) API key, if
    you want local smoke detection to override the upstream
    CAMS-derived aerosol baseline

No paid services, no cloud, no SaaS. Open-Meteo's free tier is the
default atmospheric source and doesn't require a key.

## Install

```sh
git clone https://github.com/tknx/mooncolor.git
cd mooncolor
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

## Get a Hue Bridge API key

One time, on the same LAN as the bridge:

```sh
# Press the physical link button on the bridge, then within 30 s:
curl -X POST http://<bridge_ip>/api -d '{"devicetype":"mooncolor#mylamp"}'
```

Response includes `{"success":{"username":"<40-char-hex>"}}`. Save that string as `HUE_API_KEY` — it's your API key.

Then find the numeric V1 light ID for the bulb you want to drive:

```sh
curl http://<bridge_ip>/api/<key>/lights | jq 'to_entries | map({id: .key, name: .value.name, model: .value.modelid})'
```

The dict key (e.g. `"5"`) is your `light_id`.

## Configure

```sh
cp config.example.yaml config.yaml
$EDITOR config.yaml
```

At minimum: set `location.latitude`, `location.longitude`, `location.elevation_m`, and `hue.bridge_ip` + `hue.light_id`. The example file documents every other knob inline.

## Run

### Foreground (development / testing)

```sh
export HUE_API_KEY=<your 40-char key>
# export PURPLEAIR_API_KEY=<your key>  # only if PurpleAir is enabled
python3 mooncolor.py --config config.yaml
```

### systemd (production)

```ini
[Unit]
Description=mooncolor — astronomical moonlight lamp
After=network-online.target
Wants=network-online.target

[Service]
Type=simple
User=mooncolor
Group=mooncolor
WorkingDirectory=/opt/mooncolor
Environment=PYTHONUNBUFFERED=1
EnvironmentFile=/etc/mooncolor/env
ExecStart=/opt/mooncolor/.venv/bin/python -u /opt/mooncolor/mooncolor.py --config /opt/mooncolor/config.yaml
Restart=always
RestartSec=5s
StandardOutput=journal
StandardError=journal
NoNewPrivileges=true
ProtectSystem=strict
ProtectHome=true
PrivateTmp=true

[Install]
WantedBy=multi-user.target
```

…and place `HUE_API_KEY=...` (plus optional `PURPLEAIR_API_KEY=...`) in `/etc/mooncolor/env` (mode 0600, owned by the service user).

### Docker

No `Dockerfile` is bundled but it's trivial since the dependencies are pure-Python:

```dockerfile
FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY mooncolor.py config.yaml ./
COPY bulbs/ ./bulbs/
EXPOSE 8080
CMD ["python", "-u", "mooncolor.py"]
```

Pass `HUE_API_KEY` (and optional `PURPLEAIR_API_KEY`) via `-e` / `--env-file`. Use `network_mode: host` (or an explicit LAN bridge) so the container can reach your Hue Bridge.

## Trigger the lamp

The service is dormant until something tells it to start. Two equivalent HTTP shapes — pick whichever fits your trigger source:

```sh
# Level-style
curl -X POST 'http://localhost:8080/lamp?state=on'    # start tracking
curl -X POST 'http://localhost:8080/lamp?state=off'   # stop + lamp off

# Event-style (for callers that emit state-machine transitions)
curl -X POST 'http://localhost:8080/event?event=enter'
curl -X POST 'http://localhost:8080/event?event=exit'
```

Both routes are idempotent — re-asserting the current state is a no-op. The `state=` form accepts `on/off/true/false/yes/no/1/0`.

Trigger possibilities: a cron entry ("on at sunset, off at sunrise"), a home-automation scene change ("when bedroom mood = Moon, trigger on"), a smart switch firing a webhook, a button on the wall fronted by an ESPHome device — anything that can issue an HTTP POST.

## Health endpoint

```sh
curl http://localhost:8080/healthz | jq
```

Returns uptime, current tick state, last tick result (altitude, phase, computed xy, dispatched brightness), atmospheric data freshness for every source, and Hue dispatch counters.

## How it works

The module-level docstring in `mooncolor.py` walks the physics in detail, with each function carrying a derivation comment in front of it. The short version:

```
   ┌──────────────────────┐
   │  astronomy-engine    │   moon altitude + phase fraction
   └──────────┬───────────┘
              │
   ┌──────────▼───────────┐
   │  atmospheric pollers │   AOD, PM2.5, PM10, RH, pressure
   │  Open-Meteo +        │   (each on its own background cadence;
   │  optional PurpleAir  │    in-memory cache with staleness guard)
   └──────────┬───────────┘
              │
   ┌──────────▼────────────────┐
   │  physics                  │
   │   solar SPD (Planck 5778K)│
   │   × lunar reflectance     │   (φ-blended with earthshine)
   │   × atmospheric trans     │   (Rayleigh × aerosol × ozone,
   │                           │    Kasten-Young airmass)
   │   integrate × CIE 1931    │   → XYZ → xy
   │   clamp to bulb gamut     │   → final xy
   └──────────┬────────────────┘
              │
   ┌──────────▼───────────┐
   │  Hue Bridge V1 PUT   │   {on, xy, bri, transitiontime}
   └──────────────────────┘
```

## Tuning notes

The physics constants are calibrated against literature references — they should not be touched without good reason. The values you might legitimately want to tune:

  - `tick.max_brightness` / `tick.min_brightness` — taste / use-case
  - `tick.period_s` — faster ticks give smoother color motion; slower
    saves API calls. 5 min is a fine sweet spot.
  - `tick.transition_s` — how long the bulb fades between frames
  - `atmospheric.purpleair_trigger_ratio` — how aggressively local
    smoke detection overrides the upstream baseline
  - `atmospheric.*_default` — what conservative values to fall back
    to if a source is unreachable. The defaults assume a clean / clear
    night; bias toward your local climate.

## Limits

  - Single bulb per service instance. To drive multiple identical
    bulbs in the same room, generalize `hue.light_id` to a list and
    fan-out the PUT. Different bulbs with different gamuts → run
    multiple instances.
  - Hue Bridge V1 only. The V2 CLIP API is HTTPS-only and many home
    setups have certificate-handling friction; V1 is plain HTTP on
    your LAN, no TLS, fast and reliable for this use case.
  - 5-min tick cadence implies up to 5-min lag in atmospheric color
    when conditions change rapidly. Drop `tick.period_s` to 60 or 120
    if you want faster response — Hue V1 caps at ~10 cmds/sec, you're
    vastly under-limit.

## Acknowledgements

The physics is built on standard atmospheric optics + astronomical ephemeris work. Specific calibrations cited inline in `mooncolor.py`:

  - **CIE 1931** color-matching functions, 10 nm tabulation
    (CIE 015:2018, Table T.4)
  - **ASTM E490 / Planck blackbody at 5778 K** for the solar SPD
  - **ROLO** lunar disk reflectance model (Stone & Kieffer)
  - **Kasten & Young 1989** airmass formula
  - **Bogumil 2003** ozone Chappuis-band cross-sections
  - **Hänel growth** correction for hygroscopic aerosols
  - **Barkjohn 2021 / EPA** PurpleAir PM2.5 correction
  - **astronomy-engine** for ephemeris (https://github.com/cosinekitty/astronomy)
  - **Open-Meteo** air-quality and forecast APIs (https://open-meteo.com/)
  - **PurpleAir** sensor network API (https://api.purpleair.com/)

## License

[Parity Public License 7.0.0](LICENSE.md) — free use + share, but derivatives must be open-source under the same terms.
