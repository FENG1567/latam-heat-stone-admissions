# Main tables

## Table 1

Coverage, diagnostic counts and exposure descriptors for the two primary cohorts. Calendar and exposure-eligible unit-days are distinct; missing health unit-days and zero-admission eligible days are shown separately. Admission totals are full-health-calendar records, not unique people; they include all calendar unit-days rather than only exposure-eligible days. E3 percentiles describe pooled eligible unit-days. The target increment summarizes the already-defined unit-specific P95-minus-median contrasts; its across-unit median and quartiles use empirical type-7 quantiles. 

| Characteristic | Brazil, primary 2006–2020 | Mexico, preferred 2015–2022 |
| --- | --- | --- |
| Geographic units | 450 | 32 |
| Calendar unit-days | 2,465,550 | 93,504 |
| Calendar dates | 5,479 | 2,922 |
| Exposure-eligible unit-days | 2,464,415 | 93,440 |
| Missing-health unit-days | 1,135 | 0 |
| Zero-admission unit-days (health calendar) | 990,206 | 7,300 |
| Zero-admission eligible unit-days | 990,206 | 7,289 |
| Mean E3 temperature (°C) | 23.34 | 20.60 |
| Median E3 temperature (°C) | 24.04 | 20.83 |
| P05 E3 temperature (°C) | 16.17 | 12.34 |
| P95 E3 temperature (°C) | 28.49 | 28.26 |
| Unit-specific P95-minus-median increment, median (IQR; °C) | 3.51 (2.84–3.96) | 3.68 (2.96–4.51) |
| K80 admission records | 2,173,811 | 318,331 |
| K81 admission records | 1,270,715 | 230,020 |
| Urinary N20–N23 admission records | 1,108,174 | 70,911 |

## Table 2

Aligned country-specific marginal RRs and directly estimated urinary-versus-K80 RRRs for the same primary joint-fit unit sets, using each unit’s P95 versus its shared median. Marginal RRs retain endpoint-specific pooling weights; dividing them does not recover the directly pooled RRR. Intervals are the reported calendar-synchronized Bartlett HAC14 inference with the modeled DL between-unit component, conditional on observed weights and heterogeneity. k, I² and positive-unit direction are reported unit-level summaries. 

| Country / endpoint | RR or direct RRR (95% CI) | k | I² (%) | Direction >0 (%) |
| --- | --- | --- | --- | --- |
| BR &#124; K80 | 1.007 (0.985–1.029) | 450 | 33.7 | 53.1% |
| BR &#124; Urinary N20–N23 | 1.105 (1.089–1.122) | 450 | 29.1 | 82.2% |
| BR &#124; Direct urinary-versus-K80 RRR | 1.098 (1.080–1.116) | 450 | 21.1 | 76.4% |
| MX &#124; K80 | 0.940 (0.898–0.984) | 32 | 0.0 | 12.5% |
| MX &#124; Urinary N20–N23 | 1.018 (0.961–1.079) | 32 | 21.8 | 46.9% |
| MX &#124; Direct urinary-versus-K80 RRR | 1.083 (1.033–1.136) | 32 | 30.3 | 62.5% |

**Additional panel: standardized Brazilian daily absolute-count scenario.** These reported panel summaries compare simultaneous local P95 settings with local medians; they are not observed national heat-event totals.

| Endpoint | Daily count difference (95% CI) | Included regions |
| --- | --- | --- |
| K80 | -0.2 (-8.5–8.1) | 450 |
| Urinary N20–N23 | 19.3 (16.1–22.5) | 450 |
| Urinary minus K80 | 19.5 (13.1–25.8) | 450 |

## Table 3

Sensitivity estimates for the eight-cell Brazilian grid (panel a) and rolling exposure windows (panel b). All are direct urinary-versus-K80 RRRs from DL/HAC14 pooling. Panel a factors the eight combinations by exposure, calendar specification and period. Shared calendar means shared weekday, month and smooth-trend coefficients across diagnoses; exposure coefficients remain diagnosis-specific. Year–month fixed effects are diagnosis-specific. The full-period grid is distinct from the primary diagnosis-specific smooth-trend model. Panel b reports four rolling windows for each country; Mexico uses the comparison cohort excluding 2019. The distributed-lag cumulative contrast is a distinct estimand.

### a. Eight-cell Brazilian grid

| Exposure | Calendar | Period | Direct RRR (95% CI) | k | I² (%) | Positive units (%) |
| --- | --- | --- | --- | --- | --- | --- |
| E3 | Shared calendar | 2006–2020 | 1.105 (1.092–1.118) | 450 | 36.3 | 72.2% |
| E3 | Shared calendar | 2008–2015 | 1.131 (1.115–1.148) | 446 | 28.7 | 75.1% |
| E3 | Year–month fixed effects | 2006–2020 | 1.090 (1.078–1.102) | 450 | 16.3 | 75.1% |
| E3 | Year–month fixed effects | 2008–2015 | 1.087 (1.070–1.104) | 446 | 15.4 | 63.9% |
| Emax3 | Shared calendar | 2006–2020 | 1.089 (1.076–1.101) | 450 | 32.4 | 71.8% |
| Emax3 | Shared calendar | 2008–2015 | 1.107 (1.092–1.123) | 446 | 25.4 | 70.0% |
| Emax3 | Year–month fixed effects | 2006–2020 | 1.077 (1.065–1.089) | 450 | 23.4 | 70.4% |
| Emax3 | Year–month fixed effects | 2008–2015 | 1.072 (1.056–1.089) | 446 | 17.5 | 64.6% |

### b. Rolling exposure windows

| Window | Brazil, primary cohort (2006–2020) | Mexico, comparison cohort (excluding 2019) |
| --- | --- | --- |
| E3 | RRR 1.098 (1.086–1.110)<br>k 450; I² 21.1%; positive units 76.4% | RRR 1.090 (1.043–1.140)<br>k 32; I² 18.5%; positive units 71.9% |
| E7 | RRR 1.111 (1.097–1.125)<br>k 450; I² 31.9%; positive units 75.3% | RRR 1.089 (1.038–1.142)<br>k 32; I² 19.9%; positive units 68.8% |
| E14 | RRR 1.106 (1.090–1.122)<br>k 450; I² 36.6%; positive units 71.1% | RRR 1.076 (1.030–1.124)<br>k 32; I² 0.0%; positive units 65.6% |
| E21 | RRR 1.105 (1.088–1.123)<br>k 450; I² 39.0%; positive units 68.2% | RRR 1.059 (1.011–1.109)<br>k 32; I² 0.0%; positive units 56.2% |

## Table 4

Targeted calendar, phenotype, adult-bridge, urgency and moisture comparisons. Calendar, phenotype, adult and urgency rows use unit-HAC14 contrasts with independent-unit DL pooling; moisture rows retain their labeled source inference. These differ from the synchronized primary intervals in Table 2. Confidence intervals correspond to the stated model. The urgency interaction is estimated jointly, not by dividing separately pooled strata. Table S8 gives the attempted, fitted, failed and ineligible counts, common-success denominators and selected-subset sizes for all 21 rows. The 2023 estimate is a period refit, not prediction validation.

| Scenario / contrast | Effect (95% CI) | k | I² (%) | Direction >0 (%) |
| --- | --- | --- | --- | --- |
| Brazil: calendar + holidays | 1.088 (1.076–1.101) | 422 | 24.3 | 75.1% |
| Mexico: year–month effects | 1.086 (1.034–1.141) | 32 | 30.3 | 68.8% |
| Mexico: calendar + holidays | 1.086 (1.032–1.143) | 32 | 39.5 | 62.5% |
| Mexico: 2023 period refit | 1.096 (0.970–1.238) | 30 | 50.0 | 63.3% |
| Brazil: urgent N20/K80 | 1.139 (1.111–1.166) | 411 | 25.1 | 67.4% |
| Brazil: elective N20/K80 | 0.998 (0.951–1.048) | 411 | 68.2 | 48.2% |
| Brazil: urgent/elective interaction | 1.156 (1.086–1.232) | 411 | 67.4 | 62.5% |
| Brazil: interaction, flexible calendar | 1.153 (1.102–1.206) | 242 | 34.2 | 66.9% |
| Brazil: exclude 2020, fixed targets | 1.184 (1.106–1.266) | 412 | 64.8 | 61.7% |
| Mexico: N20/K80 | 1.091 (1.036–1.150) | 32 | 33.8 | 71.9% |
| Mexico: N20–N21/K80 | 1.078 (1.028–1.131) | 32 | 29.0 | 59.4% |
| Mexico: N23/K80 | 1.150 (1.034–1.280) | 31 | 0.0 | 64.5% |
| Brazil adults: unstandardized | 1.113 (1.095–1.131) | 208 | 23.5 | 82.7% |
| Brazil adults: case-weight standardized | 1.099 (1.079–1.119) | 208 | 25.9 | 77.9% |
| Mexico adults: unstandardized | 1.065 (0.987–1.148) | 16 | 47.9 | 68.8% |
| Mexico adults: case-weight standardized | 1.055 (0.964–1.155) | 16 | 55.8 | 62.5% |
| Brazil: dew-point adjustment | 1.085 (1.073–1.098) | 450 | 19.5 | 70.0% |
| Brazil: humidity adjustment | 1.087 (1.073–1.101) | 450 | 18.6 | 68.0% |
| Mexico: dew-point adjustment | 1.077 (1.022–1.136) | 32 | 40.9 | 62.5% |
| Mexico: humidity, DL | 1.063 (1.001–1.129) | 32 | 42.9 | 56.2% |
| Mexico: humidity, REML–HK | 1.063 (0.999–1.131) | 32 | 42.9 | 56.2% |

