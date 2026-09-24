# Hourly Reference ET (ASCE Standardized)

Python script that calculates hourly standardized reference evapotranspiration for:

- **ETos**: short crop (clipped grass)
- **ETrs**: tall crop (alfalfa)

## Input

A CSV file with one row per hour:

| Column | Unit |
|---|---|
| Timestamp | `2026-06-02+13:00` or a standard date-time |
| AirTemperature | °C |
| RelativeHumidity | % |
| WindSpeed2m | m/s |
| SolarRadiation | W/m² |

Actual vapor pressure is computed from RH and air temperature (Eq. 41).

## Settings

Edit the lines marked `# SET` at the top of the script:

- `CSV_IN`, `CSV_OUT`: file names
- `LATITUDE`, `LONGITUDE`: for example `38°02'09.6"N` and `101°01'51.6"W`
- `TIME_ZONE_MERIDIAN`: `90 W` set for Central time
- `z`: station elevation (m)
- `z_w`: wind sensor height (m)

## Output

A CSV with every intermediate term, a day or night label (night is Rn < 0), and ETos and ETrs in mm for each hour with daily totals.

## Run

```bash
pip install pandas
python hourly_ETsz_ASCE.py
```

## Check and validation

The script reproduces the worked example in Appendix C (Table C-4, Greeley, CO) of the ASCE report.

## Notes

- At night and when the sun is low (β < 0.3 rad), fcd is carried forward from the last afternoon hour (Eq. 46).
- For the first night in a file, fcd is taken from the next morning unless `FCD_START` is set. Include the previous day in the file to avoid this.
- Negative night values are kept by default, as the report allows.

## Reference

1. Allen, R.G., Walter, I.A., Elliott, R.L., Howell, T.A., Itenfisu, D., Jensen, M.E., and Snyder, R.L. (2005). *The ASCE Standardized Reference Evapotranspiration Equation*. American Society of Civil Engineers, Reston, VA. https://doi.org/10.1061/9780784408056
2. ASCE-EWRI (2005). The ASCE Standardized Reference Evapotranspiration Equation. Allen, R.G. et al. (eds.), ASCE, Reston, VA. "Calculations required for hourly time-steps", Eqs. 34-67. Eq. 1 and Table 1 (Cn, Cd) are from the start of the same report. https://doi.org/10.1061/9780784408056.ch04

## Author

Udit Debangshi, PhD Scholar, Department of Agronomy, Kansas State University

Contact: uditdebangshi9251@gmail.com
