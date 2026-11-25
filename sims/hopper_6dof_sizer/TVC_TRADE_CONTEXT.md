# TVC Actuator Trade — Context for Claude Code

## Goal
Confirm the E+L **AG 9113** is acceptable for the TVC gimbal (Avionics + Structures prefer it; the previous owner called the 9313 overkill). Its angular accel per the TVC team is below GNC's desired 300–450 °/s². Run stability analysis in the sims to settle it.

## Actuator specs (E+L datasheet)
| | AG 9113 | AG 9313 |
|---|---|---|
| Stroke | ±50 mm | ±50 mm |
| Nominal force | 1000 N | 3000 N |
| Nominal speed | 60 mm/s | 30 mm/s |
| Mid-stroke length M | 360 mm | 390 mm |
| Weight | 3.9 kg | 5.2 kg |
| Current @24 V | 5.6 A | 7.7 A |
| Mech. position accuracy | <0.05 mm | <0.05 mm |
| Encoder | 12 bit/rev | 12 bit/rev |

The drive has configurable position/speed/**acceleration**/force limits (cascaded position→speed→torque control).

## TVC team sizing (method/layout not yet confirmed)
AG 9113: ±8.0° → 20.5 °/s, 235 °/s² | ±6.0° → 21.0 °/s, 250 °/s² | ±2.5° → 22.0 °/s, 275 °/s²
AG 9313: ±8.0° → 13.5 °/s, 640 °/s² | ±5.0° → 14.0 °/s, 700 °/s² | ±3.5° → 14.0 °/s, 725 °/s²
TVC says it can't reach 900 °/s².

## GNC requirement (from GNC team)
- 200 °/s² is the hard minimum. Below that, tracking breaks down (lag near/past −150° at 2 Hz with a 2° sine).
- Desired: 300–450 °/s².

## Existing trade doc (src/analysis/tvc/tvcRequirements.m and the trade write-up)
- Gimbal ±8° per axis, simultaneous; hard stop ±9.5° (TBR)
- Rate req ≥25 °/s (= 2° × 12.25 rad/s pitch crossover). The 6DOF currently caps rate at 15 °/s.
- Actuator model the LQT is tuned against: ωn = 36 rad/s, 50 ms rise, ~−27° phase at 2 Hz
- Resolution ≤0.05°, latency ≤10 ms
- Recommended layout: AG 9113, r_eff 122.5 mm, tail stock at 271 mm, mid-stroke at 0° → 29.8 °/s at load, 1.58× force margin at the stress load, ~5,900 °/s² capability
- MIL (324 runs, TVC_MODEL = ag9113/ag9313): flight demand never exceeded 1.05°, 7.9 °/s, 250 °/s² (≤100 °/s² after removing quantizer steps). Landing within 1–2 mm median of the PDR model.

## Known issues to resolve
1. **The accel numbers disagree ~25×.** The trade doc gives ~5,900 °/s² and TVC gives 235. Force can't be the limit (1000 N × 0.12 m ≈ 120 N·m vs 0.042 kg·m² engine). The likely source is the drive's accel-limit setting, reflected rotor inertia, or an assumption. Need TVC's method, arm, and single- vs dual-axis basis.
2. **The requirements table is inconsistent.** The accel row says ≥300 °/s², but its basis (25 °/s × 36 rad/s) = 900 °/s².
3. **TVC's 9113 rate (20.5 °/s) fails the 25 °/s requirement.** Either they're using the old sizer layout or the requirement needs revisiting.
4. **MIL raw peak accel (250) is above TVC's 9113 capability (235).** Confirm that's quantizer artifact.

## Key argument
What limits tracking at crossover is the linear amplitude = min(a/ω², r/ω), not accel alone. At ω = 12.25 rad/s:
- AG 9113 (TVC #s): min(1.57°, 1.67°) = **1.57°** (accel-bound)
- AG 9313 (TVC #s): min(4.26°, 1.10°) = **1.10°** (rate-bound)

Flight demand is ~0.65° amplitude at ~2 Hz (≤7.9 °/s, ~100 °/s²), so the 9113 has ~2.3× accel / ~2.6× rate margin.

## Analysis plan
1. **Describing-function sweep**: amplitude 0.25–3°, freq 0.5–5 Hz, actuator with rate + accel clipping. Find where 2 Hz phase hits −150° and where gain rolls off.
2. **Margins vs amplitude**: L(jω)·N(A,ω). Find the amplitude where PM < ~30°, and compare to flight demand.
3. **Monte Carlo with limits applied**: the previous MIL used ~5,900 °/s² capability, so a 235 °/s² limit was never exercised. Add rate/accel limits as parameters to the ag9113/ag9313 variants, set them to TVC's numbers, and rerun nominal + 50 seeds + stress for both. Track time in saturation, peak demand vs capability, attitude error, limit cycling, and landing offset.
4. **Stress cases** with 9113 limits (gusts, thrust misalignment, FDIR F1 hard-over).

## Quick DF sanity script
```matlab
A = [0.25 0.5 1 1.5 2 3];  f = 1.95;
acts = {'AG9113 TVC',20.5,235; 'AG9313 TVC',13.5,640; 'GNC floor',25,200};
for i = 1:size(acts,1)
    fprintf('\n%s\n  A(deg)  gain   phase(deg)\n', acts{i,1});
    for a = A
        [G,ph] = dfPoint(a, f, acts{i,2}, acts{i,3});
        fprintf('  %5.2f  %5.3f  %7.1f\n', a, G, ph);
    end
end
function [G,ph] = dfPoint(A, f, rmax, amax)
    wn = 36; z = 0.7; dt = 1e-4; w = 2*pi*f; T = 1/f;
    t = 0:dt:12*T-dt; th = 0; thd = 0; y = zeros(size(t));
    for k = 1:numel(t)
        thdd = wn^2*(A*sin(w*t(k)) - th) - 2*z*wn*thd;
        thdd = max(min(thdd, amax), -amax);
        thd  = max(min(thd + thdd*dt, rmax), -rmax);
        th = th + thd*dt;  y(k) = th;
    end
    m = t >= 6*T;
    s = 2*mean(y(m).*sin(w*t(m)));  c = 2*mean(y(m).*cos(w*t(m)));
    G = hypot(s,c)/A;  ph = rad2deg(atan2(c,s));
end
```
