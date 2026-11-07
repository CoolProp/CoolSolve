# Solid phases, freezing and melting: design study

> **Status: design study, nothing is implemented.** Written 2026-11-07 on `fix/library-gaps-5`. It answers one
> question of the maintainer: *CoolProp has no solid phase (`enthalpy(Water, T=-5, P=1E5)` is NaN); how do we get
> solid / freezing / melting properties, coherently with what CoolProp computes?* Register items concerned:
> `CS-BUG-WATER-NEAR-FREEZING`, `CS-GAP-ENTHALPY-FUSION`, `CS-GAP-FREEZINGPT`, `CS-GAP-INCOMP-EES-NAMES`
> (see [model_library_support.md](model_library_support.md)); `CS-GAP-HUMIDAIR-HIGHT` is the same kind of
> CoolProp-limit gap (clear error now, model later).

**Evidence tags.** **[V]** verified for this study (I read the pinned CoolProp source, ran code, or ran the
`coolsolve` binary); **[W]** read on a web page (URL in §8), not re-derived; **[A]** assumed, from memory or
inferred, to be checked before relying on it. The numerical checks used the CoolProp 8.0.0 Python wheel
(CoolSolve pins 8.1.0dev, commit `8eb3b778`: its source was read for the guards, the melting-line data are
the same), the Python `iapws` 1.5.5 package (GPLv3, used as a **test oracle only**, no code copied) and `gsw` 3.6.23
(TEOS-10, wraps GSW-C), in a scratch virtualenv outside the repository.

---

## 0. Recommendation in 12 lines

1. **Water ice first, as a CoolSolve-side Gibbs module, not a CoolProp patch.** CoolProp already contains the
   IAPWS R10-06 ice Gibbs function (`src/Ice.cpp`, `g_Ice` and derivatives, used by its humid-air code) **[V]**, but
   exposes it neither through `AbstractState` nor `PropsSI`. Re-implement it in CoolSolve (closed-form derivatives
   to 3rd order, ~60 lines), as `nh3h2o` was.
2. **It is thermodynamically coherent with CoolProp's liquid without any fudge**: along CoolProp's own ice-Ih
   melting line, `g_liquid(CoolProp) - g_ice(R10-06)` is below 0.6 J/kg up to 200 MPa (below 0.5 mK in `T`), the
   slopes agree with Clausius-Clapeyron to 1e-5 relative, and at the triple point they are equal to 1e-6 J/kg **[V]** (§3.3).
3. **Architecture: a phase router in front of the existing AbstractState block** of `evaluateCoolPropFunction`, returning
   the same triple (value, d/d input1, d/d input2); liquid and vapour states go through the unchanged CoolProp path, so
   units, reference-state offsets, clamps, error handling and the AD assembly are reused (§3.9). No regression
   outside the new regions.
4. **Boundary convention that the library needs**: `T = 0 degC, P = 101325 Pa` should be *liquid* (the native
   `CSL-0041` calls `VOLUME(Water, T=0, P=P)` during melting, and `T=-0.01` for ice) although CoolProp's melting line is at
   0.00252 degC **[V]**: a liquid band of a few mK below `T_melt(p)` (§3.8).
5. **Quick wins first (1-2 days)**: `T_Melting`, `P_Melting` (CoolProp melting lines, 30 fluids; they reproduce
   the EES help examples) and `FreezingPt` (CoolProp `T_freeze`, equal to the EES help values) **[V]**; a clear error
   for any state below the melting line.
6. **CO2 second** (Jaeger & Span 2012 Gibbs EOS; coefficients only in the paper, not read here), **PCMs third**
   as EES-language library functions on CSV tables (no C++), **seawater/ice slurry last** (TEOS-10, BSD GSW-C).
7. Total for water + CO2 + PCM: about 4-6 weeks; the complete plan 6-10 weeks (§6).

---

## 1. Problem and scope

Models that leave the CoolProp envelope on the cold side:

| Application | State range | What is needed |
|:--|:--|:--|
| Ice storage tank, `CSL-0041` (`ice_storage_tank_discharge_phase_change`) | water/ice at 0 degC, 1 atm | `u`, `v` of ice at -0.01 degC and of liquid at 0 degC; latent heat; ice fraction |
| Freezing / defrost of evaporators, frost on heat pumps (`CSL-0155`) | ice at -20..0 degC, humid air | `Enthalpy_Fusion(Water)`, `h`, `cp` of ice; ice-vapour (frost) saturation |
| Dry ice, CO2 below the triple point (216.59 K, 5.18 bar) | CO2 solid + vapour at 1 atm (194.7 K) | `h`, `s`, `v` of solid CO2, sublimation pressure/enthalpy |
| Ice slurry, secondary refrigerants below their freezing point | brine/glycol + ice | freezing curve, ice fraction, enthalpy of the slurry |
| Seawater freezing/desalination, sea ice | `SA` 0-120 g/kg, -10..0 degC | freezing temperature, brine salinity vs `T`, ice fraction |
| PCM storage (paraffins, salt hydrates) | melting range of a few K | `h(T)`, liquid fraction, densities, hysteresis |

Humid air already works below 0 degC: CoolProp's humid-air code has its own ice branch (`IceProps`, `psub_Ice`);
`humrat(AirH2O, T=-10, P=101325, R=1)` = 1.606e-3 and `dewpoint` evaluate in the shipped binary **[V]**. What is
missing is the *condensed* solid phase of the pure fluid and of its solutions.

What happens today **[V]** (shipped binary, CoolProp 8.1.0dev `8eb3b778`): `enthalpy(Water, T=-5, P=1E5)`,
`enthalpy(Water, T=0.0001, P=101325)` and `density(CarbonDioxide, T=-60, P=1E5)` all end with *"CoolProp returned
invalid result (NaN or Inf)"*. The informative message exists in CoolProp (*"For now, we don't support T [268.15 K]
below Tmelt(p) [273.153 K]"*, `HelmholtzEOSMixtureBackend.cpp`, guards at lines ~1847, 1917, 2338) but is lost: the
AbstractState path throws, `evaluateCoolPropFunction` falls back to `PropsSI`, which returns Inf (the mixture path
re-throws with the reason; the pure-fluid path does not). The short-term "clear error" is therefore one line of
plumbing: keep the AbstractState message.

---

## 2. Survey of models and libraries

| Substance | Phases / content | Formulation, range, accuracy | Licence, language | Fit with CoolSolve |
|:--|:--|:--|:--|:--|
| **Water ice Ih**: IAPWS R10-06(2009), Feistel & Wagner 2006 | solid Ih, all properties | Gibbs `g(T,p)`, 0-273.16 K, 0-210 MPa **[W]**; designed to be used with IAPWS-95 | IAPWS release (free to use, attribution **[A]**); C++ in CoolProp `Ice.cpp` (MIT) **[V]**; Python `iapws` GPLv3 **[V]**; GSW-C `gsw_gibbs_ice` BSD-3 **[V]** | **Best.** Same reference state as IAPWS-95 (CoolProp's Water). Closed-form derivatives |
| **Water melting + sublimation**: IAPWS R14-08(2011) | `p_sub(T)` 50-273.16 K; `p_melt(T)` Ih 251.165-273.16 K, ices III, V, VI, VII to 715 K **[V]** (`iapws` docstrings) | correlations consistent with IAPWS-95 and R10-06 **[W]** | in CoolProp: `psub_Ice`, `melting_line` (Water: 611.657 Pa-2.18 GPa) **[V]** | Direct: boundaries come for free |
| **Supercooled water**: IAPWS G12-15 | metastable liquid to the nucleation limit, 0-400 MPa | guideline, does not replace IAPWS-95 **[W]** | public; `iapws._Supercooled` | Optional. CoolProp's IAPWS-95 extrapolates it already: density within 4e-5 (relative) of G12-15 at 263 K, `cp` within 0.2 % **[V]** |
| **Water, whole diagram** (Ih, II, III, V, VI, VII/X + liquid to 2300 MPa): SeaFreeze, Journaux et al. 2020 | high-pressure ices, NaCl(aq) | Gibbs representations | **GPL-3**, Python **[W]** | Not embeddable in an MIT code. Use as an oracle; coefficients are in the papers. Low priority (HVAC never needs ice III+) |
| **Seawater + ice + sea ice**: TEOS-10 (Feistel 2008, IAPWS-08/09/06), GSW | Gibbs of seawater `g(SA,t,p)`, ice, freezing `t_freezing`, brine `SA_freezing_from_t`, `ice_fraction_to_freeze_seawater`, `melting_ice_into_seawater`, frazil, sea ice **[V]** (function list of `gsw`) | `t_freezing(35.165 g/kg, p=0)` = -1.9191 degC **[V]**; McDougall et al. 2014 **[W]** | GSW-C **BSD-3** **[V]**; C; MATLAB/Python/R wrappers | Good: C, no dependency, exact first derivatives; own `Jet<3>` port possible (as `nh3h2o`). Pure-water part is IAPWS-09, not IAPWS-95: g equal to 1e-3 J/kg up to 25 degC, 17 J/kg in `h` at 60 degC **[V]** |
| **CO2 solid** (dry ice): Jaeger & Span 2012 | solid CO2, sublimation and melting pressures predicted with Span-Wagner fluid | Gibbs-explicit, fitted to `cp`, `v`, alpha, kappa of the solid; equilibrium pressures "very sensitive" to g and to the fluid EOS **[W]** (abstract) | paper only (J. Chem. Eng. Data 57, 590); coefficients **not read** **[A]** | Good fit (CoolProp CO2 *is* Span-Wagner **[V]**). Needs the paper's coefficient table |
| **Argon solid**: Maltby et al. 2024 | solid Ar, quasi-harmonic Helmholtz EOS, `T` <= 300 K, `p` <= 16 GPa | JPCRD 53(4) **[W]** (blog) | paper | Out of scope; shows the pattern (solid EOS per substance) |
| **Other pure fluids** | melting line in CoolProp for 30 of 136 fluids **[V]**; EES has `P_Melting`/`T_Melting`/`P_Sublimation` for "selected fluids" **[W]** (N2: melting Span et al. 2000, sublimation Fray & Schmitt 2009) | correlations only, no solid caloric data | CoolProp MIT | Lines yes (immediate); solid properties no |
| **REFPROP** | melting-line state flag; "does not calculate solid properties" **[W]** | - | NIST, closed | Not a source |
| **Cantera** | `StoichSubstance`: incompressible fixed-composition solid, thermo from polynomials **[W]** (partly) | no ice EOS found | BSD, C++ | Pattern for a *generic* solid (`h(T)`, `rho`), not for the water phase diagram |
| **EES** | `Ice` fluid: Hyland & Wexler (1983, NBS/NRC steam tables) below 0 degC, ice `k` Coles 1954, sound speed Vogt 2008, melting/sublimation Wagner et al. 1994; `Phase$` gives `solid` only for ice **[W]** | reference: h = s = 0 for saturated liquid at 0 degC (for `Water` too) **[W]** | closed | The compatibility target (§5) |
| **Modelica MSL** | `IF97_Utilities.BaseIF97.IceBoundaries`: melting/sublimation pressure only **[W]** (search result) | - | open | Not a source |
| **PCM data**: Rubitherm RT, SP; RAL-GZ 896; ASTM C1784 | `h(T)` per heating/cooling, melting/congealing ranges, storage capacity with tolerance (RT 28 HC: 250 kJ/kg, +-7.5 % **[W]**, garbled extract) | tables, hysteresis models (Barz & Sommer 2018) **[W]** | data sheets | Tabulated layer (§4) |
| **Brines / glycols**: Melinder (IIR 2010) in CoolProp INCOMP | `T_freeze(x)` (a function of concentration only) **[V]**; below it CoolProp **raises** (*"Your temperature 253.15 is below the freezing point of 258.574"*) **[V]**; `MITSW` (sea water) has no `T_freeze` and `Tmin` = 273.15 K **[V]** | `T_freeze`: CaCl2 20 % = -18.260 degC, MEG 20 % = -7.949 degC: **identical to the EES help values** **[V]** | CoolProp MIT | Freezing curve for free; no ice fraction |

Findings that shape the plan: (a) for **water**, everything needed is public, already consistent with CoolProp
and partly inside CoolProp; (b) for **other substances** there is **no open library that gives a whole phase
diagram**: only per-substance papers (CO2, Ar), the GPL water-only SeaFreeze, and TEOS-10 (water-salt); (c) the
generic substitute is **line + data** (melting line from CoolProp, solid `cp`/`rho`/`dh_fus` from tables).

---

## 3. Coupling with CoolProp, made rigorous

### 3.1 What CoolProp has and has not **[V]**

- HEOS has **no solid EOS**. For 30 fluids it has a **melting line** (`melting_line(iT|iP, ...)`; Water 611.657 Pa-2.18 GPa,
  CO2 5.18e5-8.23e8 Pa, N2, Ar, CH4, ...) used only as a guard: PT flashes refuse `T < T_melt(p) - 1 mK`.
  `has_melting_line()` makes the list discoverable at run time. CoolProp 8.0.0 adds `MeltingCaloric` (fluid properties
  along the line, to seed H,S flashes) **[W]**: still liquid-side only.
- **No sublimation line** for pure fluids (only `psub_Ice` for water in the humid-air code). Below the triple
  point, `p < p_triple` and `T < Tmin` is refused (*"we don't support p [...] below ptriple [...] when T [...] is less than Tmin"*).
- Those refusals are **guards, not EOS limits**: `update(DmassT_INPUTS, rho, T)` (`DHSU_T_flash`, no melting or `Tmin`
  test in the pinned source) evaluates the vapour at 200-270 K for water and 180-216 K for CO2, and the
  *metastable liquid* of IAPWS-95 down to 248 K **[V]** (a small Newton on `rho` to meet `p`). The config key
  `DONT_CHECK_PROPERTY_LIMITS` also lifts the PT guards (PR #3031 **[W]**) but is global and would weaken every error
  path: a dedicated bypass is safer. Risk: it depends on CoolProp internals, so a test must pin the behaviour at each bump
  (§6).
- Reference state: default IAPWS (u = s = 0 for saturated liquid at 273.16 K): `enthalpy(Water, T=0.01, x=0)` = 0.6118 J/kg
  in the binary **[V]**. **EES references `Water` to 0 for saturated liquid at 0 degC [W]**, so liquid enthalpies
  already differ from EES by about 41.6 J/kg (u of the saturated liquid at 273.15 K is -42.2 J/kg in the IAPWS
  reference **[V]**). Nothing below changes that; it matters only for absolute comparisons with EES (§6).

### 3.2 Gibbs framework and reference anchoring

For every phase `phi` use the specific Gibbs energy `g_phi(T,p)` (mass basis) as the potential:

```
v = g_p        s = -g_T        h = g - T g_T        u = h - p v        cp = -T g_TT
alpha = g_Tp/g_p    kappa_T = -g_pp/g_p    c_sound = g_p * sqrt( g_TT / (g_Tp^2 - g_TT g_pp) )
d(h)/dT|p = cp    d(h)/dp|T = g_p - T g_Tp    d(s)/dp|T = -g_Tp    d(v)/dp|T = g_pp    d(v)/dT|p = g_Tp
```

For the fluid, CoolProp gives `h, s, v, cp, alpha, kappa`, hence `g_f = h - T s` and all the derivatives above. A solid
EOS gives its own `g_S`. Equilibrium needs equal chemical potentials, and **`g` is defined up to `a - b T`**
(a reference-state change) *per phase*. Two rules follow:

1. Offsets that CoolSolve applies to results (`refState.h_offset/s_offset` in `postProcess`) apply to all phases after
   the property call: the solid must be returned in CoolProp's *raw* reference, and the router stays in front of
   `postProcess`.
2. **Anchor the solid at the triple point.** Choose `(a_S, b_S)` so that (i) `g_S(Tt,pt) = g_L(Tt,pt)` and (ii)
   `h_L - h_S = dh_fus` at `(Tt,pt)`. When no calorimetric `dh_fus` is trusted, take it from the **Clapeyron relation on
   CoolProp's own melting line**, `dh_fus = Tt (v_L - v_S) (dp_m/dT)_t`, so the solid is tied to what CoolProp
   already says. Caveat: `dh_fus` then inherits the error of `v_S` amplified by `v_L/(v_L - v_S)`: about 1 % on
   `dh_fus` per 0.1 % on `rho_S` for water, 0.5 % for CO2 **[V]** (CO2, CoolProp slope 0.2138 K/MPa at the triple point: with the
   EES help volumes `v_L` = 8.485e-4 and `v_S` = 6.471e-4 m3/kg, i.e. `rho_S` = 1545 kg/m3, Clapeyron gives `dh_fus` = 204.0 kJ/kg,
   EES's `h_L - h_S` is 204.9; with `rho_S` = 1510 kg/m3 it would be about 189 kJ/kg: calorimetric values for dry ice of that
   size are in the literature **[A]**). For water the anchor is practically exact: the R10-06 `g00` of CoolProp's
   `Ice.cpp` differs from the 2009 revision (the `iapws` constant) by 1.1e-4 J/kg **[V]** (that CoolProp's is the pre-2009 value is **[A]**);
at the triple point `g_L` (CoolProp) and `g_ice` (2009 constants) are both 0.611782 J/kg **[V]**.

### 3.3 Are CoolProp's melting line and the solid EOS consistent? Measured for water **[V]**

Ice from R10-06 (2009, `iapws` oracle; CoolProp's own `g_Ice` is the same function with the older `g00`), liquid
from CoolProp HEOS at `T_m(p)` of CoolProp's melting line (Wagner et al., R14-08):

| `p` [MPa] | `T_m` [K] | `g_L - g_ice` [J/kg] | `dT_m/dp` CoolProp [K/MPa] | `dv/ds` (Clapeyron) [K/MPa] |
|--:|--:|--:|--:|--:|
| 0.0007 (near triple pt) | 273.16000 | < 1e-5 | -0.074263 | -0.074262 |
| 0.101325 | 273.15252 | 0.0002 | -0.074295 | -0.074293 |
| 10 | 272.40165 | -0.0068 | -0.077408 | -0.077410 |
| 50 | 269.05946 | 0.052 | -0.089609 | -0.089610 |
| 100 | 264.20875 | -0.189 | -0.104377 | -0.104376 |
| 200 | 252.31670 | -0.505 | -0.133243 | -0.133246 |

`|dg|` below 0.6 J/kg is below 0.5 mK in `T` (divide by `ds` = 1.22 kJ/kg-K). Other facts: `dh_fus` = 333.4449 kJ/kg at
the triple point and 333.4265 kJ/kg at `T_m(1 atm)` = 273.15252 K (TEOS-10 `latentheat_melting` gives the same
333426.52 J/kg); `v_ice` = 1.090842e-3 m3/kg at 273.14 K. **Sublimation line**: with the vapour from CoolProp's EOS evaluated below
the triple point (DmassT bypass, §3.1) and `p_sub` of R14-08, `g_V - g_ice` = +0.74, -5.4, -6.9, +4.9 J/kg at 270, 250,
230, 200 K, i.e. `|dg|/(R T)` <= 7e-5: **the sublimation pressure is consistent to 0.007 %** (the IAPWS-95 ideal part is
only extrapolated below 273.16 K). The sublimation correlation in CoolProp (`psub_Ice`) uses the R14-08 coefficients.
**CO2 [V for the line, A for the solid]**: the melting line is Span-Wagner's quadratic in `Theta = T/216.592 - 1` (216.592-330 K), slope
at the triple point 0.2138 K/MPa; with EES's `h_L - h_S` = 204.9 kJ/kg it implies `dv_fus` = 2.023e-4 m3/kg, EES's `DeltaV_fusion`
help value is 2.014e-4 (0.4 %): CoolProp's line and EES's are the same physics. The Clapeyron test of the CO2 *solid* cannot be run
before the Jaeger-Span coefficients are available.

### 3.4 When the solid EOS and CoolProp's line disagree: decision rules

1. **The phase boundary is CoolProp's melting correlation** (it is what CoolProp itself guards with, it was fitted to
   the melting data, it is explicit and monotone per branch), not a root of `g_S = g_L`. Properties come from each EOS;
   `g` jumps at the boundary by `eps(p)` (water: < 0.6 J/kg), `h`, `s`, `v` by the physical latent amounts. The jump
   `h_L - h_S` computed from the two EOS is the latent heat, and a unit test compares it with the Clapeyron value
   (`T dv dp_m/dT`; criterion 1e-3 relative for water, 1e-2 per new substance).
2. If `eps` or the Clapeyron residual of a new substance exceeds the criterion, first **re-anchor** (rule of 3.2);
   if still above, document the substance as "solid properties approximate near melting" and keep the CoolProp
   line for the boundary. Do **not** hide it by adding a smooth correction to `g_S` (it would be absorbed into `v_S`).
3. **Sublimation boundary** (no CoolProp line): water uses R14-08 (`psub_Ice`); another substance uses the root of
   `g_S(T,p) = g_V(T,p)` (1D Newton from the ideal-gas start), which is how Jaeger & Span validated CO2 **[W]**.
   Both curves must meet CoolProp's melting line at `(Tt, pt)`: tested to 1e-9 relative on `p`.
4. No solid EOS (30 fluids with a line): the line is exposed (`T_Melting`, `P_Melting`), solid properties are an error.

### 3.5 Phase determination

By minimum `g` in principle, by the three lines in practice (identical to within `eps/ds`, no solve):

| Region (pure substance with a solid model) | Condition | Properties from |
|:--|:--|:--|
| liquid, vapour, supercritical | `p > pt`: `T >= T_m(p) - band`; or `T > Tt` and not solid | CoolProp, **unchanged path** |
| solid | `p >= pt`: `T < T_m(p) - band`; `T <= Tt`, `p >= p_sub(T)` | solid EOS |
| solid + liquid | `T = T_m(p)`, `p > pt` | lever rule between both EOS |
| solid + vapour | `p = p_sub(T)`, `T < Tt` | lever rule: solid EOS + vapour |
| triple point | `(Tt, pt)` | two phase fractions needed: error unless given (out of scope) |
| other solids | water `p` above about 209.9 MPa **[A]** and `T < T_m(p)`: ice III, V, VI, VII | error naming the unimplemented phase |

The Ih branch of water ends at the Ih-III-liquid point (251.165 K **[V]**, the lower limit of the `iapws` Ih melting correlation; about 209.9 MPa **[A]**);
`melting_line(iT, iP, p)` above it describes other ices (non-monotone `T_m`: 252.3 K at 200 MPa, 254.96 K at 300 MPa **[V]**).

### 3.6 Two-phase regions and the meaning of `quality`

- `quality` / `x` **stays the vapour mass fraction** (CoolProp and EES). **EES-confirmed**: the help example
  `temperature(Ice, v=1.088e-3, x=0)` = -16.85 degC; ice on the R14-08 sublimation line with that `v` is at **-16.843 degC [V]**,
  so in EES `x` is the vapour fraction in the solid-vapour region and `x = 0` the saturated solid. Same meaning here.
- Add the **solid mass fraction** `xs` (property `SolidFraction`, input designator `xs=`): 1 for the solid, 0 for
  liquid and vapour, in S+L `xs`, in S+V `xs = 1 - x`. The liquid fraction is `1 - xs - x`. (EES has no such property;
  the name is ours.) `quality` is 0 for a single-phase solid or liquid (as the existing clamp does for subcooled states).
- Extensive properties are linear in the fractions: `phi = xs phi_S + x phi_V + (1-xs-x) phi_L`, evaluated at the
  boundary state. Input pairs: S+L: `(P,h)`, `(P,s)`, `(P,u)`, `(P,v)`, `(P,xs)`, `(T,h)`, `(T,xs)` (the Ih branch
  `p_m(T)` is single-valued on 251.165-273.16 K); S+V: `(T,x)`, `(P,x)`, `(T,h)`, `(P,h)` ... The pair `(T,P)` **on** the line is ambiguous and
  resolves to the liquid/vapour side (§3.8), with a warning.
- Single-phase solid: any two of `T, P, h, s, u, v/rho` by a 2D Newton on `(T,p)` with the analytic Hessian above.

### 3.7 Derivatives, kinks, solver

- Inside a region everything is analytic (formulas of 3.2): no finite difference. The evaluator's
  "analytic vs forward-FD" consistency check of the CoolProp path is **not applied** to the router (a FD step across the
  boundary would produce `dh_fus/dT`).
- Across the boundary: for `(T,P)` inputs `h, s, v, u` **jump** (`dh_fus` = 333 kJ/kg for water) and `cp` is a delta; for
  `(P,h)` inputs `T(p,h)` is continuous and piecewise smooth with a flat segment (`dT/dh` = `1/cp_S`, 0, `1/cp_L`;
  `dxs/dh = -1/dh_fus` inside; `dT_m/dp` = -0.074 K/MPa): a kink of the kind the vapour dome already
  produces with CoolProp `(P,h)` inputs. Models must use `(P,h)` (or the state `h` with an `IF`-free formulation), never
  `h = h(T,P)` solved for `T` inside the latent range: that equation has no solution. A diagnostic should say so.
- The AD gradient is assembled exactly as today (`dResult/dInput1`, `dResult/dInput2` times the input gradients);
  `xs`-dependent outputs add `dResult/dxs`. Because T_m depends on `p` the S+L derivatives w.r.t. `p` carry the
  chain-rule term `dT_m/dp` (closed form from the correlation).
- Trial-point clamps (`P_MIN_SI` = 1000 Pa, `T_MIN_SI` = 50 K, and the fluid's own minimum in `evaluator.cpp`) must be lowered
  for fluids with a solid model: the S+V region of water lives at `p` < 611 Pa (R14-08 reaches 50 K). Otherwise a
  frost state is silently moved to 611.65 Pa.

### 3.8 The boundary band (a requirement found in the library)

CoolProp's melting temperature at 101325 Pa is 273.15252 K (0.0025 degC, = TEOS-10 `t_freezing(0,0,0)` **[V]**), so
the textbook state `T = 0 degC, P = 101325 Pa` is, strictly, **ice by 2.5 mK** and CoolProp refuses it (guard at `T_m - 1 mK`,
**[V]**: `enthalpy(Water, T=0.0001, P=101325)` fails today). The native `CSL-0041` evaluates
`v_liq = VOLUME(Water, T=T_liq_ref, P=P)` with `t_liq_ref = max(t_liq_0, t_liq)` = 0 degC during the whole melting
phase (`t_liq` < 0 there), and `v_ice` / `u_ice_0` at -0.01 degC (12.5 mK below `T_m`) **[V]** (native file, lines 57, 82-84). The
library variant replaces them by "ice at 0 degC" (0.001091) and "liquid water at 0 degC" (0.001000), and takes `u_ice_0` from the
EES stored solution of the sibling file (-333451.5211 J/kg); that EES returns *liquid* for the 0 degC call is the intention of the model
and the variant, **not** something read in an EES output **[A]**. The convention that satisfies both calls:

> the solid region is `T < T_m(p) - band`, with a default **band of 5 mK** (any value in 2.7-12.5 mK works for these two
> calls; configurable); inside the band the state is the *liquid* (metastable by at most 5 mK, evaluated through the
> DmassT bypass: the error is `cp * 5 mK` = 21 J/kg at most, 0.006 % of the latent heat), and a note is logged.

How EES places its own boundary for `p` != 1 atm is **[A] unknown**; the band should be a documented parameter.

### 3.9 Integration in the code (precedent: `nh3h2o`)

`nh3h2o` duplicates the whole evaluator because it has three inputs and a composition. A solid phase has the
**same two-input, one-output shape as CoolProp**, so the cheapest hook is **inside the AbstractState block** of
`evaluateCoolPropFunction` (`src/evaluator.cpp`, around the `state->update(ipair, ...)` / `first_partial_deriv` code):

```
fluid has a SolidModel?  (registry: Water/Steam/R718/Ice -> IceIh, CarbonDioxide/R744 -> DryIce, later others)
  -> router(T,p | pair, SolidModel, bands) -> {Solid, SolidLiquid, SolidVapour, Fluid}
  -> Fluid: fall through to the CoolProp code that exists today (bit-identical results)
  -> else: value, d/dinput1, d/dinput2 from the module  -> postProcess(...)   (units, offsets, v = 1/rho, AD)
```

New files in the style of `nh3h2o`: `include/coolsolve/solid_phase.h` (the `SolidModel` interface: `g` with derivatives to
3rd order, validity box, anchor constants), `src/solid_ice.cpp`, `src/solid_router.cpp`, hooks in `evaluator.cpp` for
property calls, `phase$` and the pure-fluid constants (`Enthalpy_Fusion`, `T_Melting`, ...). No `taylor_jet.h` is needed for
ice (the derivatives of the complex-log terms are closed forms: with `f(th) = (t-th)ln(t-th) + (t+th)ln(t+th) - 2t ln t - th^2/t`,
`f' = ln((t+th)/(t-th)) - 2th/t`, `f'' = 1/(t-th) + 1/(t+th) - 2/t`, `f''' = 1/(t-th)^2 - 1/(t+th)^2`); it would be used for TEOS-10
seawater (`Jet<3>` in `S, T, p`).

---

## 4. Generic layers: PCMs and solutions

### 4.1 PCMs (tabulated, EES-language library, no C++)

A PCM is not a fluid with an EOS: it is a pair of condensed phases with a melting *range*, so the model is the **enthalpy
method** `h(T)` (equivalently the apparent heat capacity `cp_app = dh/dT`), as in the data sheets (heat storage capacity
within temperature ranges, separate heating and cooling curves, RAL-GZ 896, ASTM C1784 **[W]**). Proposal, delivered through the
existing machinery (`$INCLUDE` of a library model; tables `<stem>-<table>.csv`; linear `INTERPOLATE`, §11 of the
language reference, with the slope of the active interval as exact derivative):

- one CSV per product, `pcm_rt28hc-hT.csv` with `T_C, h_heat_J_per_kg, h_cool_J_per_kg` (cooling column optional) and a
  one-row metadata table `rho_s, rho_l, cp_s, cp_l, T_solidus, T_liquidus, k_s, k_l, source` (`LOOKUP$ROW` key);
- FUNCTIONs `PCM_h(Name$, T, Branch$)`, `PCM_T(Name$, h, Branch$)` (inversion = interpolation with swapped columns, strictly
  monotone `h` required; checked), `PCM_xl(Name$, T)` (liquid fraction `(h - h_s(T))/(h_l(T) - h_s(T))`, with `h_s`, `h_l` the sensible
  extrapolations of the solid and liquid branches), `PCM_s`, `PCM_rho(Name$, T)`. Entropy `s = s0 + integral(dh/T)` evaluated on the
  table (needed for exergy, otherwise optional). Pressure dependence is ignored (Clapeyron estimate for a paraffin:
  about 0.2 K/MPa **[A]**, 0.02 K at 1 bar).
- Parametric fall-back (no table): `T_s, T_l, L, cp_s, cp_l` with a linear or smoothed (sine, Gaussian) melting fraction.
- **Hysteresis is a state, not a function**: it needs the history of the cell, which an algebraic solve does not have. Offer
  only the two branches and, for dynamic models, the *static hysteresis* formulation of Barz & Sommer (2018, IJHMT 127, 701 **[W]**)
  with an `INTEGRAL` state for the liquid fraction. Steady-state models choose a branch explicitly.
- Smoothness: the tables are C0 with piecewise-constant `cp_app`; a user-set `smoothing` width (K) rounds corners for Newton.
  The pure-substance case is **not** smoothed (the exact physics has a jump, §3.7).
- Data: only open or redistributable data; Rubitherm and others are not copied into the repository unless their terms allow it
  **[A: not checked]**: the format and a generator from a user's DSC curve are the deliverable.

### 4.2 Solutions (brines, glycols, seawater)

- **Freezing curve**: CoolProp INCOMP `T_freeze(x)` is available for the Melinder solutions (values equal to the EES help)
  **[V]**: `FreezingPt` is a wiring job. Its inverse `x_b(T)` (brine concentration in equilibrium with ice at `T`) is monotone and is
  solved by 1D Newton/bisection. Lever rule: ice mass fraction `w_ice = 1 - x_0/x_b(T)` for `T <= T_freeze(x_0)` and `T >= T_eutectic`.
- **Enthalpy of an ice slurry**: `h = w_ice h_ice(T) + (1-w_ice) h_brine(T, x_b)` needs a brine enthalpy that includes the
  **enthalpy of dilution**. CoolProp's incompressible solutions integrate `cp` at fixed concentration: *no enthalpy of mixing*
  (already documented for `LiBr`/`SeaWater`, code `C004`). So this formula is not an energy balance without a correction;
  the engineering literature builds enthalpy-phase diagrams for the purpose (Melinder & Granryd 2005 **[W]**). Verdict: ice fraction and
  freezing temperature yes; slurry *enthalpy* only with a stated accuracy (to be assessed against that paper), tagged `C004`-style.
- **Seawater via TEOS-10** is the rigorous case: `g(SA,T,p)` and `g_ice` are one consistent set, the equilibrium
  `mu_w(SA,t,p) = g_ice(t,p)` fixes the freezing curve (check at `SA` = 35.165 g/kg: `mu_w` = -2.2587 J/g, `mu_ice` = -2258.70 J/kg
  **[V]**), the brine salinity at `T`, the ice fraction and the heat of melting ice into seawater (329 928.65 J/kg at 35.165 g/kg
  against 333 426.52 for pure water **[V]**). Plan: a `SeaWaterIce` solution built on a `Jet<3>` port of the Feistel 2008
  function (BSD GSW-C as source of the coefficients and as oracle), valid 0-120 g/kg **[A]**. It is **not** the CoolProp `MITSW` model
  (Sharqawy et al., no freezing, `T` >= 0 degC): both would coexist, with the TEOS-10 route as the one that supports ice; a unification
  is a decision for the maintainer (§7). Note the two reference states: TEOS-10 seawater has h = s = 0 at 35.165 g/kg, 0 degC, 1 atm
  **[V]** (pure water in that scale: h(0 degC) = 61.0 J/kg), so enthalpies must not be mixed with CoolProp's without the (constant)
  offset.
- Below the eutectic (NaCl -21.2 degC, CaCl2 -49.8 degC **[A]**) salt hydrates appear: error "below the eutectic", not modelled.

---

## 5. API proposal (EES-compatible)

| Call | Behaviour | EES basis |
|:--|:--|:--|
| `enthalpy(Water, T=-5, P=1E5)` (any property, any input pair) | ice Ih value, IAPWS reference like the liquid (`h` = -343.75 kJ/kg and `v` = 1.0900e-3 m3/kg at -5 degC, R10-06 **[V]**) | EES evaluates its `Water` below 0 degC with ice (stored `INTENERGY(Water, T=-0.01, P=1E5)` = -333451.5211 J/kg, `CS-BUG-WATER-NEAR-FREEZING`) |
| `enthalpy(Water, T=0, P=101325)` | **liquid** (band, §3.8) | `CSL-0041` |
| `enthalpy(Ice, ...)`, `Phase$(Ice, ...)` | alias of Water; with ambiguous `(v, x)` inputs the lower-temperature solution (EES help) | EES `Ice` fluid |
| `temperature(Ice, v=1.088e-3, x=0)` | -16.843 degC (EES help -16.85) | EES help |
| `Enthalpy_Fusion(Fluid)` | `h_L - h_S` at the normal melting point (`T_m` at 1 atm; at the triple point if `pt` > 1 atm, CO2); Water 333 426.5 J/kg | EES help (EES: 333 605.9169, +0.054 %, other EOS) |
| `T_Melting(Fluid, P=)`, `P_Melting(Fluid, T=)`, `P_Sublimation(Fluid, T=)`, `DeltaV_Fusion(Fluid)` | CoolProp melting line (30 fluids), sublimation: water R14-08, CO2 Span-Wagner correlation or `g_S = g_V` root; `DeltaV_Fusion` by Clapeyron | EES help; names already in the parser list |
| `FreezingPt(Fluid, C=)` | CoolProp `T_freeze` (Melinder solutions); a `T` input accepted and ignored | EES help |
| `SolidFraction(Fluid, ...)` and input `xs=` | §3.6 (new, no EES counterpart) | - |
| `Phase$` | adds `solid`, `solid_liquid`, `solid_gas` (existing style: `supercritical_gas`); inside the band `liquid` | EES: `solid` for ice only |
| `Quality(...)` | vapour mass fraction; 0 for single-phase solid/liquid | EES |

Checks of the EES help examples **[V]**: `P_Melting(Water, T=270 K)` = 39.313 MPa (EES 39.31); `T_Melting(ice, P=20 MPa)` = 271.612 K
(EES 271.6); `P_Sublimation(CarbonDioxide, T=-78.5 C)` with the (recalled **[A]**) Span-Wagner sublimation form
`ln(p/pt) = (Tt/T)(a1 th + a2 th^1.9 + a3 th^2.9)`, `th = 1 - T/Tt`, `a` = -14.740846, 2.4327015, -5.3061778 = 101.03 kPa (EES 101);
the same form gives 101.33 kPa at 194.686 K, the sublimation point at 1 atm.

**Errors** (never NaN, names the call as `CS-BUG-INCOMP-SOLIDS` does): a fluid with a melting line but no solid model:
*"enthalpy(Nitrogen, T=-230, P=1E5): below the melting/sublimation temperature (63.15 K at this pressure); CoolSolve has no
solid-phase properties for Nitrogen (available: Water, CarbonDioxide); `T_Melting(Nitrogen, P=...)` gives the boundary"*; ice outside the
validity (*"p = 300 MPa: ice III/V region, not implemented"*); `T` < 50 K or `p` above the box; the triple point without two fractions;
a state inside the band logs one note (`C005`-style) per run. A switch in `coolsolve.conf` (`solidPhases = off|on`, `meltingBand = 0.005`)
restores the plain error: the change turns an error into a value, which may hide a modelling mistake.

---

## 6. Options, order, validation

| Option | Content | Size | Risks |
|:--|:--|:--|:--|
| **A. Minimal** | Clear error for every state below the melting line (keep the AbstractState message; use `has_melting_line()`); `T_Melting`, `P_Melting`, `FreezingPt` (INCOMP), `T_Triple` kept; `Enthalpy_Fusion(Water)` as a constant | ~150 lines C++, ~150 lines tests, docs; 1-2 days | none; no solid values |
| **B. Water** | A + `SolidModel` interface, ice Ih (R10-06), sublimation line, router, S+L and S+V with all input pairs, `xs`, `Ice` alias, `phase$`, `Enthalpy_Fusion`, `DeltaV_Fusion`, band | ~1000 lines src, ~700 tests (for scale: `nh3h2o` is ~2800 src + 1250 tests **[V]**); 2-3 weeks | CoolProp internals (DmassT bypass); band convention; trial-point clamps; Newton on kinks; EES reference offset |
| **C. B + CO2 + other substances** | CO2 via Jaeger-Span (needs the paper's table), sublimation by root, anchor by Clapeyron; optional per-fluid "line + `cp`/`rho` table" solids for other fluids | +400 src, +300 tests, 1-2 weeks for CO2; each other fluid is a data job | Coefficients transcribed from a paywalled paper; dry-ice literature spread (§3.2); fluid EOS vintage |
| **D. Complete** | C + PCM library (EES-language, CSV; ~1 week, independent of C++) + brine/ice slurry (~1 week) + TEOS-10 seawater-ice with `Jet<3>` (~2 weeks) + optional metastable liquid (`phase` selection) | 6-10 weeks in all | `C004`-type enthalpy inconsistency of CoolProp brines; two seawater models; data licensing of PCM sheets |

**Recommended order**: (0) short-term clear error; (1) option A quick wins; (2) water ice (B); (3) CO2 (C); (4) PCM tables, which can start
**in parallel** since it is library work; (5) seawater/ice slurry (D); (6) other fluids only on demand.

**Validation strategy.**
1. *Unit/oracle*: R10-06 Table 6 check values (the release PDF was **not** read here, the IAPWS page gave only the landing page: take them from the
   PDF `https://iapws.org/technical-guidance/release/Ice-2009.download`), tolerance 1e-9 relative; three independent oracles available offline: CoolProp's own `g_Ice` /
   `h_Ice` / `rho_Ice` (header `CoolProp/fluids/Ice.h`, linkable in the test binary **[V]**), GSW-C `gsw_gibbs_ice`, and `iapws` (offline script that writes the numbers into
   the test file; the three agree to the printed digits in the probes **[V]**).
2. *Thermodynamic identities* (property-based): Maxwell relations, `h = g - T g_T`, analytic derivative vs central FD (1e-8), `(T,P)` -> `(P,h)` -> `T` round
   trips for all pairs, lever-rule linearity, continuity at the boundaries (jump = `dh_fus` to 1e-12), Clapeyron along the melting line (<= 1e-3 relative; the table of 3.3
   shows 1e-5), `|g_S - g_L|` < 1 J/kg on the melting line up to 200 MPa, < 10 J/kg on the sublimation line, triple-point equality 1e-9.
3. *EES references* (help examples of §5 and stored solutions): TM-0095 `MSTH051024exercice4.EES`, `u_ice_0` = -333451.5211 J/kg: the R10-06 value on the EES reference
   (liquid at 0 degC = 0) is -333449.39 J/kg **[V]**, i.e. 2.1 J/kg (6e-6 relative) from EES (Hyland-Wexler vs R10-06); in CoolSolve's reference it is -333491.59 J/kg, 40 J/kg lower, so either compare
   differences or shift both phases by the 41.6 J/kg reference offset (§3.1). `Enthalpy_Fusion(Water)` within 0.1 % of 333 605.9.
4. *Library model `CSL-0041`*: run the native file (its `INTENERGY`/`VOLUME` calls as written) and compare with the verified variant: `v_ice` 1.0908e-3 m3/kg at 273.14 K (variant 0.001091), `v_liq` 1.0002e-3
   at 0 degC (variant 0.001000), melting end at 30 460 s, final `u` 20 934.2 J/kg (the README trajectory; `u_ice_0` differs by 40 J/kg through the reference, 0.012 % on the melting time, about 4 s); the
   0 degC liquid call is the band test. Its `c_ice` = 2200 J/kg-K is an exercise constant (R10-06: 2096.6 at 273.14 K).
5. *Regression*: the whole suite unchanged (router is a no-op above `T_m` + band); a test of the CoolProp bypass and of the melting-line data to catch a pin bump.
6. *Human check*: GUI/`phase$` sweep of a `T-s` or `p-h` diagram of water across 0 degC.

---

## 7. Open questions for the maintainer

1. Band default (5 mK) and the name of the switch; should `Water` at `T = 0, P = 101325` really be liquid (needed by `CSL-0041`)?
2. Align the Water reference state with EES (0 at saturated liquid at 0 degC, 41.6 J/kg), or keep CoolProp's? Independent of this study but it decides how EES stored values are compared.
3. Name of the solid fraction (`xs` / `SolidFraction` proposed) and whether `quality` of a solid should be 0 (proposed) or an error.
4. Supercooled liquid as an explicit option (the machinery is the same bypass; IAPWS G12-15 as oracle)? Default stays "stable phase = ice".
5. One seawater model or two (`MITSW` without ice, TEOS-10 with ice)?
6. Ice III-VII: not planned (SeaFreeze is GPL; HVAC does not need it).

## 8. Sources

- IAPWS R10-06(2009) ice Ih: https://iapws.org/relguide/Ice-2009.html ; Feistel & Wagner, JPCRD 35, 1021 (2006), doi:10.1063/1.2183324 ; PDF https://iapws.org/technical-guidance/release/Ice-2009.download
- IAPWS R14-08(2011) melting and sublimation: https://iapws.org/relguide/MeltSub.html ; doi:10.1063/1.3657937 ; IAPWS G12-15 supercooled water: https://iapws.org/relguide/Supercooled.html
- TEOS-10 / GSW-C: https://github.com/TEOS-10/GSW-C (licence: https://raw.githubusercontent.com/TEOS-10/GSW-C/master/LICENSE) ; `gsw_ice_fraction_to_freeze_seawater`: https://teos-10.org/pubs/gsw/html/gsw_ice_fraction_to_freeze_seawater.html ; McDougall et al., JPO 44, 1751 (2014)
- Python `iapws` (GPLv3): https://github.com/jjgomera/iapws ; SeaFreeze (GPL-3): https://pypi.org/project/seafreeze/
- Jaeger & Span, J. Chem. Eng. Data 57, 590 (2012), doi:10.1021/je2011677 : https://fis.tu-dresden.de/portal/en/publications/equation-of-state-for-solid-carbon-dioxide-based-on-the-gibbs-free-energy(96592c23-b9d9-401b-8d93-4804f35c5369).html ; Maltby et al., JPCRD 53(4) (2024): https://oivindw.com/2024/12/25/new-article-presenting-a-solid-phase-equation-of-state-for-argon-with-state-of-the-art-accuracy/
- CoolProp changelog (8.0.0: melting-line items): https://coolprop.org/coolprop/changelog.html ; `Ice.cpp`: https://raw.githubusercontent.com/CoolProp/CoolProp/master/src/Ice.cpp ; `MeltingCaloric`: https://coolprop.org/_static/doxygen/html/_melting_caloric_8cpp_source.html ; incompressible/Melinder: https://secondarycoolantprops.readthedocs.io/
- EES help: https://fchartsoftware.com/ees/eeshelp/steam_nbs.htm (Ice), `steam_iapws.htm` (reference state), `enthalpy_fusion.htm`, `hs804.htm` (P_SUBLIMATION), `hs805.htm` (P_MELTING), `hs821.htm` (T_MELTING), `hs967.htm` (DELTAV_fusion), `hs158.htm` (PHASE$), `hs651.htm` (FREEZINGPT), `nitrogen.htm`, `hs731.htm` (sea water) — all under `https://fchartsoftware.com/ees/eeshelp/`
- REFPROP: https://refprop-docs.readthedocs.io/en/latest/GUI/Menu%20Commands/Calculate%20Menu/specifiedstatepoints.html ; Cantera `StoichSubstance`: https://cantera.org/science/phase-thermo.html ; Modelica `IceBoundaries`: https://doc.modelica.org/Modelica%204.1.0/Resources/helpDymola/Modelica_Media_Water_IF97_Utilities_BaseIF97_IceBoundaries.html
- PCM: https://www.rubitherm.eu/media/products/datasheets/Techdata_-RT28HC_EN_21012026.PDF ; Barz & Sommer, IJHMT 127, 701 (2018), doi:10.1016/j.ijheatmasstransfer.2018.08.032 ; ASTM C1784 ; RAL-GZ 896 (https://pcm-ral.org, not read) ; Melinder & Granryd, Int. J. Refrig. 28, 13 (2005): https://iifiir.org/en/fridoc/using-property-values-of-aqueous-solutions-and-ice-to-estimate-ice-123737

*Not verified in this study*: the R10-06 and R14-08 PDFs (equations and Table 6), the Jaeger-Span coefficients, EES's boundary rule for `p` != 1 atm and
its CO2 solid values beyond the help page, GSW validity ranges, Rubitherm numbers (garbled extract), the licence terms of IAPWS releases and PCM data sheets.
