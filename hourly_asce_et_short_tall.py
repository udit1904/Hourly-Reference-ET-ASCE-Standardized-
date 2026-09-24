# -*- coding: utf-8 -*-

# =====================================================================
# Hourly standardized reference ET (ETsz) for the short crop (ETos, grass) and the tall crop (ETrs, alfalfa).
#
# Source: ASCE-EWRI (2005). The ASCE Standardized Reference Evapotranspiration
#         Equation. Allen, R.G. et al. (eds.), ASCE, Reston, VA.
#         "Calculations required for hourly time-steps", Eqs. 34-67.
#         Eq. 1 and Table 1 (Cn, Cd) are from the start of the same report. https://doi.org/10.1061/9780784408056.ch04
#
# Variable names follow the report (P, gamma, Delta, es, ea, Rns, Rnl, Rn,
# fcd, Rso, Ra, dr, delta, J, omega, omega1, omega2, omega_s, Sc, b, beta, G, u2).
# =====================================================================

import math
import re
from datetime import datetime, timedelta
import pandas as pd

# ---------------- Input / output files ----------------
CSV_IN = "/content/Book3.csv"                                                                               # SET
CSV_OUT = "Book1_ETsz_hourly.csv"                                                                           # SET
COLUMNS = {
    "timestamp": "Timestamp",
    "T":   "AirTemperature",           # deg C
    "RH":  "RelativeHumidity",         # %
    "u_z": "WindSpeed2m",              # m/s
    "Rs":  "SolarRadiation",           # W/m2
}

# ---------------- Site ----------------

LATITUDE = """38°02'09.6"N"""       # station location as N and W                                           # SET
LONGITUDE = """101°01'51.6"W"""                                                                             # SET
TIME_ZONE_MERIDIAN = "90 W"   # centre of your standard time zone:
                              # 75 W Eastern, 90 W Central,
                              # 105 W Mountain, 120 W Pacific


def to_deg(text, positive, negative):
    """'101.031 W' or 101°01'51.6"W -> decimal degrees with the sign set by the letter.
    positive: the letter that gives a + value; negative: the letter that gives a - value."""
    s = text.strip()
    if "-" in s:
        raise ValueError(f"Write '{text}' without a minus sign. The letter sets the direction.")
    letter = s[-1].upper()
    numbers = [float(n) for n in re.findall(r"\d+(?:\.\d+)?", s)]
    if not 1 <= len(numbers) <= 3:
        raise ValueError(f"Cannot read '{text}'.")
    numbers += [0.0] * (3 - len(numbers))
    value = numbers[0] + numbers[1] / 60 + numbers[2] / 3600     # degrees + minutes/60 + seconds/3600
    if letter == positive:
        return value
    if letter == negative:
        return -value
    raise ValueError(f"'{text}': the letter must be {positive} or {negative}.")


LAT_DEG = to_deg(LATITUDE, "N", "S")            # + north (Eq. 49 uses this)
L_m = to_deg(LONGITUDE, "W", "E")               # + west, as Eq. 55 defines L_m
L_z = to_deg(TIME_ZONE_MERIDIAN, "W", "E")      # + west, as Eq. 55 defines L_z
z = 895.8               # station elevation above sea level, m                                               # SET
z_w = 2.0               # height of the wind measurement, m                                                  # SET

# ---------------- Time convention of the logger ----------------
# Eq. 55: t = standard clock time at the MIDPOINT of the period, after correcting for any daylight savings shift.

TIMESTAMP_IS_END_OF_HOUR = True   # True: stamp 13:00 = period 12:00-13:00                                   # SET
CLOCK_IS_DAYLIGHT_SAVING = False  # True: logger clock is CDT, not CST                                       # SET
t_1 = 1.0                         # length of the period, h (Eq. 53-54)                                      # SET

# ---------------- fcd before the first beta >= 0.3 in the file ----------------

#   None  -> use fcd from the first period of the morning with beta >= 0.3
#            (the report, p. 36, allows the following-morning value for night hours) 0.05-1.0 -> use this value

FCD_START = None                                                                                             # SET

SET_NEGATIVE_TO_ZERO = False      # report (p. 45): negative night values                                    # SET
                                  # may be kept; "it may be appropriate to
                                  # retain the negative values"

# =====================================================================
# Constants given in the report
# =====================================================================
lambda_ = 2.45          # latent heat of vaporization, MJ/kg (inverse = 0.408)
G_sc = 4.92             # solar constant, MJ m-2 h-1 (Eq. 48)
sigma = 2.042e-10       # Stefan-Boltzmann constant, MJ K-4 m-2 h-1 (Eq. 44)
alpha = 0.23            # albedo (Eq. 43)

# Eq. 1 / Table 1, hourly time step: Cn, Cd for daytime and nighttime
# Eq. 65-66: G as a fraction of Rn. Nighttime = Rn < 0 (p. 44).
REFERENCE = {
    "ETos": {"Cn": 37.0, "Cd_day": 0.24, "Cd_night": 0.96, "G_day": 0.10, "G_night": 0.50},  # Eq. 65a, 65b
    "ETrs": {"Cn": 66.0, "Cd_day": 0.25, "Cd_night": 1.70, "G_day": 0.04, "G_night": 0.20},  # Eq. 66a, 66b
}

# =====================================================================
# Psychrometric and atmospheric variables
# =====================================================================
def eq34_P(z):
    """Eq. 34: mean atmospheric pressure, kPa."""
    return 101.3 * ((293 - 0.0065 * z) / 293) ** 5.26


def eq35_gamma(P):
    """Eq. 35: psychrometric constant, kPa/C."""
    return 0.000665 * P


def eq36_Delta(T):
    """Eq. 36: slope of the saturation vapor pressure curve, kPa/C."""
    return 2503 * math.exp(17.27 * T / (T + 237.3)) / (T + 237.3) ** 2


def eq37_es(T):
    """Eq. 37: saturation vapor pressure e0(T), kPa."""
    return 0.6108 * math.exp(17.27 * T / (T + 237.3))


def eq41_ea(RH, T):
    """Eq. 41: actual vapor pressure from mean RH and mean T for the hour, kPa."""
    return RH / 100.0 * eq37_es(T)

# =====================================================================
# Extraterrestrial radiation and sun position
# =====================================================================
def eq49_radians(deg):
    """Eq. 49: decimal degrees to radians."""
    return math.pi / 180.0 * deg


def eq52_J(D_M, M, Y):
    """Eq. 52: day of the year J."""
    return (D_M - 32 + int(275 * M / 9) + 2 * int(3 / (M + 1))
            + int(M / 100 - (Y % 4) / 4 + 0.975))


def eq50_dr(J):
    """Eq. 50: inverse relative Earth-Sun distance."""
    return 1 + 0.033 * math.cos(2 * math.pi / 365 * J)


def eq51_delta(J):
    """Eq. 51: solar declination, rad."""
    return 0.409 * math.sin(2 * math.pi / 365 * J - 1.39)


def eq58_b(J):
    """Eq. 58: b, rad."""
    return 2 * math.pi * (J - 81) / 364


def eq57_Sc(b):
    """Eq. 57: seasonal correction for solar time, h."""
    return 0.1645 * math.sin(2 * b) - 0.1255 * math.cos(b) - 0.025 * math.sin(b)


def eq55_omega(t, L_z, L_m, Sc):
    """Eq. 55: solar time angle at the midpoint of the period, rad."""
    return math.pi / 12 * ((t + 0.06667 * (L_z - L_m) + Sc) - 12)


def eq53_omega1(omega, t_1):
    """Eq. 53: solar time angle at the beginning of the period, rad."""
    return omega - math.pi * t_1 / 24


def eq54_omega2(omega, t_1):
    """Eq. 54: solar time angle at the end of the period, rad."""
    return omega + math.pi * t_1 / 24


def eq59_omega_s(phi, delta):
    """Eq. 59: sunset hour angle, rad."""
    return math.acos(-math.tan(phi) * math.tan(delta))


def eq56_limits(omega1, omega2, omega_s):
    """Eq. 56: limit omega1 and omega2 to sunrise (-omega_s) and sunset (omega_s)."""
    if omega1 < -omega_s:
        omega1 = -omega_s
    if omega2 < -omega_s:
        omega2 = -omega_s
    if omega1 > omega_s:
        omega1 = omega_s
    if omega2 > omega_s:
        omega2 = omega_s
    if omega1 > omega2:
        omega1 = omega2
    return omega1, omega2


def eq48_Ra(dr, omega1, omega2, phi, delta):
    """Eq. 48: extraterrestrial radiation for the period, MJ m-2 h-1."""
    return 12 / math.pi * G_sc * dr * ((omega2 - omega1) * math.sin(phi) * math.sin(delta)
                                       + math.cos(phi) * math.cos(delta)
                                       * (math.sin(omega2) - math.sin(omega1)))


def eq62_beta(phi, delta, omega):
    """Eq. 62: sun angle above the horizon at the midpoint of the period, rad."""
    return math.asin(math.sin(phi) * math.sin(delta) + math.cos(phi) * math.cos(delta) * math.cos(omega))


def eq47_Rso(z, Ra):
    """Eq. 47: clear-sky solar radiation, MJ m-2 h-1."""
    return (0.75 + 2e-5 * z) * Ra


# =====================================================================
# Net radiation
# =====================================================================
def eq45_fcd(Rs, Rso):
    """Eq. 45: cloudiness function, with 0.3 <= Rs/Rso <= 1.0 so 0.05 <= fcd <= 1.0."""
    ratio = min(max(Rs / Rso, 0.3), 1.0)
    return 1.35 * ratio - 0.35


def eq43_Rns(Rs):
    """Eq. 43: net short-wave radiation, MJ m-2 h-1."""
    return (1 - alpha) * Rs


def eq44_Rnl(fcd, ea, T):
    """Eq. 44: net long-wave radiation, MJ m-2 h-1. T_K = T + 273.16."""
    T_K_hr = T + 273.16
    return sigma * fcd * (0.34 - 0.14 * math.sqrt(ea)) * T_K_hr ** 4


def eq42_Rn(Rns, Rnl):
    """Eq. 42: net radiation, MJ m-2 h-1."""
    return Rns - Rnl


# =====================================================================
# Wind and the standardized equation
# =====================================================================
def eq67_u2(u_z, z_w):
    """Eq. 67: wind speed adjusted to 2 m, m/s."""
    return u_z * 4.87 / math.log(67.8 * z_w - 5.42)


def eq1_ETsz(Delta, Rn, G, gamma, Cn, Cd, T, u2, es, ea):
    """Eq. 1: standardized reference ET, mm/h."""
    return ((0.408 * Delta * (Rn - G) + gamma * Cn / (T + 273) * u2 * (es - ea))
            / (Delta + gamma * (1 + Cd * u2)))


# =====================================================================
# Driver
# =====================================================================
def parse_timestamp(raw):
    """Reads '2026-06-02+13:00' (the '+HH:MM' is the hour of day) and normal formats."""
    s = str(raw).strip()
    m = re.match(r"^(\d{4}-\d{2}-\d{2})\+(\d{1,2}):(\d{2})$", s)
    if m:
        d, hh, mm = m.groups()
        return datetime.strptime(f"{d} {int(hh):02d}:{mm}", "%Y-%m-%d %H:%M")
    return pd.to_datetime(s).to_pydatetime().replace(tzinfo=None)


def period_midpoint(stamp):
    """Standard clock time t at the midpoint of the period (Eq. 55 definition) and its date."""
    start = stamp - timedelta(hours=t_1) if TIMESTAMP_IS_END_OF_HOUR else stamp
    if CLOCK_IS_DAYLIGHT_SAVING:
        start = start - timedelta(hours=1)
    mid = start + timedelta(hours=t_1 / 2)
    t = mid.hour + mid.minute / 60 + mid.second / 3600
    return t, mid.date()


def solar_terms(t, day, phi):
    """Eqs. 47-62 for one period."""
    J = eq52_J(day.day, day.month, day.year)       # Eq. 52
    dr = eq50_dr(J)                                # Eq. 50
    delta = eq51_delta(J)                          # Eq. 51
    b = eq58_b(J)                                  # Eq. 58
    Sc = eq57_Sc(b)                                # Eq. 57
    omega = eq55_omega(t, L_z, L_m, Sc)            # Eq. 55
    omega_s = eq59_omega_s(phi, delta)             # Eq. 59
    omega1, omega2 = eq56_limits(eq53_omega1(omega, t_1),   # Eq. 53
                                 eq54_omega2(omega, t_1),   # Eq. 54
                                 omega_s)                   # Eq. 56
    Ra = eq48_Ra(dr, omega1, omega2, phi, delta)   # Eq. 48
    beta = eq62_beta(phi, delta, omega)            # Eq. 62
    Rso = eq47_Rso(z, Ra)                          # Eq. 47
    return dict(J=J, dr=dr, delta=delta, b=b, Sc=Sc, t_mid=t, omega=omega, omega_s=omega_s,
                omega1=omega1, omega2=omega2, Ra=Ra, Rso=Rso, beta=beta)


def compute(rows):
    """rows: list of dicts with stamp, T, ea, Rs (MJ m-2 h-1), u_z. Returns list of result dicts."""
    phi = eq49_radians(LAT_DEG)                    # Eq. 49
    P = eq34_P(z)                                  # Eq. 34
    gamma = eq35_gamma(P)                          # Eq. 35

    # Pass 1: sun position and fcd where beta >= 0.3 rad (Eq. 45)
    recs = []
    for r in rows:
        t, day = period_midpoint(r["stamp"])
        s = solar_terms(t, day, phi)
        s["fcd_calc"] = eq45_fcd(r["Rs"], s["Rso"]) if s["beta"] >= 0.3 and s["Rso"] > 0 else None
        recs.append({**r, **s})

    # Value for the night at the start of the file (see FCD_START)
    first_calc = next((x["fcd_calc"] for x in recs if x["fcd_calc"] is not None), 1.0)
    fcd_carry = FCD_START if FCD_START is not None else first_calc

    out = []
    for x in recs:
        # Eq. 45 when beta >= 0.3, otherwise Eq. 46: fcd = fcd(beta>0.3) from the prior period
        if x["fcd_calc"] is not None:
            fcd, src = x["fcd_calc"], "Eq.45"
            fcd_carry = fcd
        else:
            fcd, src = fcd_carry, "Eq.46"

        T, ea, Rs = x["T"], x["ea"], x["Rs"]
        Delta = eq36_Delta(T)                      # Eq. 36
        es = eq37_es(T)                            # Eq. 37
        Rns = eq43_Rns(Rs)                         # Eq. 43
        Rnl = eq44_Rnl(fcd, ea, T)                 # Eq. 44
        Rn = eq42_Rn(Rns, Rnl)                     # Eq. 42
        u2 = eq67_u2(x["u_z"], z_w)                # Eq. 67
        night = Rn < 0                             # p. 44: nighttime is Rn < 0

        res = dict(x)
        res.pop("fcd_calc")
        res.update(P=P, gamma=gamma, Delta=Delta, es=es, VPD=es - ea, fcd=fcd, fcd_eq=src,
                   Rns=Rns, Rnl=Rnl, Rn=Rn, u2=u2, period="night" if night else "day")
        for name, c in REFERENCE.items():
            G = (c["G_night"] if night else c["G_day"]) * Rn          # Eq. 65 / 66
            Cd = c["Cd_night"] if night else c["Cd_day"]              # Table 1
            et = eq1_ETsz(Delta, Rn, G, gamma, c["Cn"], Cd, T, u2, es, ea) * t_1   # Eq. 1
            if SET_NEGATIVE_TO_ZERO:
                et = max(et, 0.0)
            res[f"G_{name[2:]}"] = G
            res[f"Cd_{name[2:]}"] = Cd
            res[f"{name}_mm"] = et
        out.append(res)
    return out


if __name__ == "__main__":
    df = pd.read_csv(CSV_IN)
    df.columns = [c.strip() for c in df.columns]
    rows = []
    for _, r in df.iterrows():
        T = float(r[COLUMNS["T"]])
        rows.append(dict(
            stamp=parse_timestamp(r[COLUMNS["timestamp"]]),
            T=T,
            RH=float(r[COLUMNS["RH"]]),
            ea=eq41_ea(float(r[COLUMNS["RH"]]), T),                  # Eq. 41
            u_z=float(r[COLUMNS["u_z"]]),
            Rs=float(r[COLUMNS["Rs"]]) * 0.0036,                     # W/m2 -> MJ m-2 h-1 (x 3600 s / 1e6)
        ))
    rows.sort(key=lambda x: x["stamp"])
    res = pd.DataFrame(compute(rows))
    res.insert(0, "Timestamp", df.sort_values(
        COLUMNS["timestamp"], key=lambda s: s.map(parse_timestamp))[COLUMNS["timestamp"]].values)
    res = res.drop(columns="stamp")

    cols = (["Timestamp", "period", "T", "RH", "ea", "es", "VPD", "Delta", "P", "gamma", "u2",
             "J", "t_mid", "Sc", "omega", "omega1", "omega2", "omega_s", "delta", "dr", "beta",
             "Ra", "Rso", "Rs", "fcd", "fcd_eq", "Rns", "Rnl", "Rn",
             "G_os", "Cd_os", "ETos_mm", "G_rs", "Cd_rs", "ETrs_mm"])
    res = res[cols]
    res.to_csv(CSV_OUT, index=False)

    print(res[["Timestamp", "period", "beta", "fcd", "fcd_eq", "Rn", "ETos_mm", "ETrs_mm"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    for p in ["day", "night"]:
        sub = res[res["period"] == p]
        print(f"{p:5s}: ETos {sub['ETos_mm'].sum():.2f} mm, ETrs {sub['ETrs_mm'].sum():.2f} mm ({len(sub)} h)")
    print(f"total: ETos {res['ETos_mm'].sum():.2f} mm, ETrs {res['ETrs_mm'].sum():.2f} mm ({len(res)} h)")
    print(f"Saved: {CSV_OUT}")

