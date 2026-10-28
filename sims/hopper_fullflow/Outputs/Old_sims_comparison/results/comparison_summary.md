# vehicle_sim.py vs FullFlow model (1dof thrust profile)

Statistics over profile time 1.0-29.5 s, on vehicle_sim.py's 0.1 s grid.
"Plateau" columns use only the last sample of each profile hold (30 points), where both models have settled after a step.

| Quantity | Units | Old mean | New mean | Mean difference (new - old) | Mean |difference| [%] | Max |difference| [%] | Plateau mean |difference| [%] | Plateau max |difference| [%] |
|---|---|---|---|---|---|---|---|---|
| Thrust | N | 1041 | 1041 | -0.3867 | 0.2053 | 2.458 | 0.04585 | 0.3736 |
| Chamber pressure | psia | 180.2 | 178.8 | -1.404 | 1.318 | 10.92 | 0.6728 | 1.352 |
| Mixture ratio | - | 1.8 | 1.8 | 2.803e-16 | 2.325e-14 | 7.401e-14 | 2.549e-14 | 4.934e-14 |
| Ox mass flow | kg/s | 0.3608 | 0.3593 | -0.00149 | 0.4125 | 0.4213 | 0.4126 | 0.4213 |
| Fuel mass flow | kg/s | 0.2004 | 0.1996 | -0.0008279 | 0.4125 | 0.4213 | 0.4126 | 0.4213 |
| Isp | s | 188.6 | 189.4 | 0.7601 | 0.4884 | 2.919 | 0.457 | 0.8257 |
| C* | m/s | 1486 | 1483 | -3.039 | 0.2046 | 0.2516 | 0.2033 | 0.2176 |
| LOX tank pressure | psia | 501.1 | 501.8 | 0.7175 | 1.289 | 6.232 | 1.187 | 3.629 |
| Fuel tank pressure | psia | 501 | 501.3 | 0.271 | 1.054 | 4.389 | 1.035 | 3.453 |
| COPV pressure | psia | 3193 | 3714 | 520.7 | 18.56 | 41.87 | 18.89 | 41.42 |
| COPV temperature | K | 265.2 | 277.5 | 12.28 | 4.791 | 10.41 | 4.876 | 10.31 |
| LOX ullage temperature | K | 150.5 | 151.3 | 0.7388 | 0.9432 | 3.395 | 0.9515 | 2.692 |
| Fuel ullage temperature | K | 266.2 | 270.4 | 4.256 | 1.689 | 4.109 | 1.717 | 3.835 |
| LOX injector dP | psi | 35.59 | 35.42 | -0.1695 | 0.4767 | 2.06 | 0.4568 | 0.5726 |
| Fuel injector dP | psi | 34.98 | 34.75 | -0.2284 | 0.6507 | 1.953 | 0.6392 | 0.7065 |
| Regen dP | psi | 6.288 | 6.206 | -0.08249 | 2.912 | 94.53 | 0.8232 | 0.8447 |
| OTV CdA | mm^2 | 5.546 | 5.499 | -0.04685 | 1.668 | 8.998 | 1.671 | 5.518 |
| FTV CdA | mm^2 | 4.051 | 4.023 | -0.02799 | 1.617 | 9.592 | 1.059 | 4.737 |
| LOX remaining | kg | 13.74 | 13.75 | 0.01334 | 0.1216 | 0.3773 | 0.1238 | 0.3668 |
| Fuel remaining | kg | 6.539 | 6.547 | 0.007412 | 0.1489 | 0.4878 | 0.1513 | 0.4716 |
| COPV gas remaining | kg | 4.099 | 7.696 | 3.597 | 88.57 | 106.6 | 88.82 | 106.3 |

## Thrust vs commanded thrust (most recent profile point)

| Model | RMS error [N] | Max error [N] |
|---|---|---|
| vehicle_sim.py | 4.5 | 25.8 |
| FullFlow | 0.0 | 0.0 |

## Values at profile time 29.5 s

| Quantity | vehicle_sim.py | FullFlow |
|---|---|---|
| Thrust [N] | 988.2 | 987.6 |
| Chamber pressure [psia] | 172.8 | 171.3 |
| Mixture ratio [-] | 1.8 | 1.8 |
| Ox mass flow [kg/s] | 0.3457 | 0.3443 |
| Fuel mass flow [kg/s] | 0.1921 | 0.1913 |
| Isp [s] | 187.3 | 188 |
| C* [m/s] | 1485 | 1482 |
| LOX tank pressure [psia] | 491.7 | 503.8 |
| Fuel tank pressure [psia] | 497.9 | 497.2 |
| COPV pressure [psia] | 2165 | 3056 |
| COPV temperature [K] | 238.7 | 263.2 |
| LOX ullage temperature [K] | 146.3 | 150.1 |
| Fuel ullage temperature [K] | 248.1 | 257.6 |
| LOX injector dP [psi] | 32.21 | 32.07 |
| Fuel injector dP [psi] | 31.66 | 31.46 |
| Regen dP [psi] | 5.692 | 5.645 |
| OTV CdA [mm^2] | 5.187 | 5.073 |
| FTV CdA [mm^2] | 3.703 | 3.694 |
| LOX remaining [kg] | 8.694 | 8.726 |
| Fuel remaining [kg] | 3.735 | 3.754 |
| COPV gas remaining [kg] | 3.425 | 7.061 |
