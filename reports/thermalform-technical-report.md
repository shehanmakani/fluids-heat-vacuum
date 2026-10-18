# ThermalForm: designing PFAS-free single-phase immersion cooling fluids

*An IntelliForm module*

**Shehan Makani · ChemeNova LLC · October 2026**

## Abstract

ThermalForm designs non-fluorinated dielectric fluids for single-phase immersion cooling of AI servers. Blends of published base stocks (polyalphaolefins, synthetic esters and a renewable isoparaffin) are scored with open property models, chip-level heat-sink physics, a 1 MW facility loop and a signal-integrity model, all under the Open Compute Project (OCP) immersion-fluid specification. Three lead blends pass the specification at worst-case permittivity and beat a PAO benchmark at the chip and/or the pump. Lead A runs a 1 kW chip 0.40 °C cooler and saves about $3,100 per MW a year in pumping; Lead G is fully readily biodegradable and saves about $5,500. The governing unknown is ester permittivity at 20–40 GHz, which is unpublished and sets how much ester a blend can carry.

## 1. Problem and specification

High-density AI racks are moving to immersion cooling. Many high-performance dielectric fluids are fluorinated, and PFAS restrictions are tightening. Hydrocarbon and ester fluids avoid PFAS but trade off heat transfer, flammability and dielectric behavior.

Design constraints come from the OCP Base Specification for immersion fluids:

| Requirement | Limit |
| --- | --- |
| Figure of merit FOM1 = k(βc<sub>p</sub>ρ²/(μk))<sup>0.2813</sup> | Tier 1 > 35, Tier 2 > 45 |
| FOM2 = k(ρc<sub>p</sub>/k)<sup>1/3</sup> | > 19 |
| Dynamic viscosity at 50 °C | < 0.015 Pa·s |
| Dielectric constant, 20 MHz to 40 GHz | ≤ 2.3 |
| Flash point | > 150 °C |
| Pour point | < −30 °C |
| Mineral oils | to be avoided |

## 2. Data

Eleven fluids from supplier data sheets (`thermalform/data/base_stocks.csv`): PAO 2/4/6, three Priolube synthetic esters, a natural ester, a mineral oil, a renewable isoparaffin, and two transformer-grade reference esters. Every value carries its basis (data sheet, estimate, or proxy). Estimates are flagged in the file: ester densities, some thermal conductivities and heat capacities, and all GHz permittivities. Prices are placeholders.

## 3. Property models

| Property | Model |
| --- | --- |
| Kinematic viscosity | Refutas blending index; Walther (ASTM D341) temperature dependence |
| Flash point | Hu–Burns blending index |
| Density, heat capacity | Ideal mixing; linear temperature slopes |
| Thermal conductivity | Filippov equation, n components |
| Thermal expansion β | Chemistry-class values |
| Dielectric constant | Lichtenecker logarithmic mixing; central and worst case at GHz |

Esters are polar. Their permittivity is published at low frequency (about 3 at 100 kHz) but not at 20–40 GHz, where it falls by an unknown amount. The model carries a central estimate and a worst case equal to the low-frequency value. Hydrocarbons are non-polar and frequency-flat (about 2.07–2.08, consistent with published PAO-type fluids at 20–40 GHz).

## 4. Chip-level heat transfer

A copper plate-fin heat sink (90 × 90 mm, 30 mm fins, 0.6 mm thick, 1.5 mm gap) on a 1 kW chip, fed 6 L/min of fluid at 40 °C.

- Laminar developing flow in the fin channels (Teertstra correlation), Reynolds number about 16.
- Fin efficiency from the classical straight-fin solution; caloric resistance from fluid heat capacity.
- Friction from the Muzychka–Yovanovich developing-flow model plus entrance and exit losses.

Because the flow is laminar, thermal conductivity and volumetric heat capacity dominate and viscosity matters little at the chip. The turbulent Mouromtseff figure of merit used in round 1 misranks fluids here and was retired.

## 5. Facility loop

A 1 MW module: chevron plate heat exchanger (VDI/Martin correlation, ε-NTU counterflow), 60 m of 150 mm pipe with minor losses, 65% pump efficiency, 10 K tank rise and New Jersey electricity at $0.1629/kWh. Fluid pumping comes to about 2% of IT power, roughly $20,000–34,000 per MW a year depending on the fluid. This corrected round 2, which had treated pumping as negligible.

## 6. Signal integrity

A microstrip designed for 50 Ω in air (Hammerstad–Jensen) on a Dk 3.5 laminate, with the fluid as cover dielectric through a filling factor. A PAO-class fluid (Dk 2.08) pulls the trace to 47.3 Ω; a fluid at Dk 3.0 pulls it to 45.3 Ω, the edge of a ±10% tolerance.

## 7. Design method

Blends of up to three stocks (each at least 5%) were searched with a multi-objective evolutionary optimizer (NSGA-II) under the OCP limits, in four scenarios: standard and premium flash class, central and worst-case ester permittivity. The optimizer belongs to ChemeNova's IntelliForm platform and is not in this repository; the resulting recipes are evaluated with the open models in `thermalform/evaluate_leads.py`.

## 8. Results

Reproduced by `python evaluate_leads.py` (1 kW chip, 40 °C supply, 6 L/min, 1 MW loop):

| Fluid | Recipe (mass) | Chip ΔT vs PAO | Max chip W at 70 °C | Pumping vs PAO, $/yr per MW | Worst-case Dk | Flash °C (model) |
| --- | --- | --- | --- | --- | --- | --- |
| PAO benchmark | PAO 4 88% + PAO 2 12% | 0 | 1,307 | 0 | 2.08 | 198 |
| Lead A | PAO 4 40% + Priolube 1976 35% + PAO 2 25% | −0.40 °C | 1,330 | −3,122 | 2.28 | 186 |
| Lead G | Renewable isoparaffin 62% + Priolube 1976 38% | −0.25 °C | 1,321 | −5,477 | 2.29 | 170 |
| Lead P | PAO 4 81% + Priolube EF 7010 19% | −0.28 °C | 1,323 | +118 | 2.22 | 216 |
| All-ester | EF 7010 52% + Priolube 1976 48% | −1.31 °C | 1,386 | −4,600 | 2.90 | 219 |

Lead G is 100% readily biodegradable. Lead P meets a 210 °C premium flash class. The all-ester blend is the best performer thermally but fails the dielectric limit unless ester GHz permittivity turns out low.

## 9. Uncertainty

Monte Carlo with 300 paired draws perturbed density (1%), conductivity (3%), heat capacity (4%) and viscosity (2–3%), with larger spreads on estimated values, and sampled ester GHz permittivity uniformly between 2.2 and its low-frequency value.

| Fluid | P(cooler than PAO) | P(less pumping) | P(Dk ≤ 2.3) |
| --- | --- | --- | --- |
| Lead A | 0.97 | 0.97 | 1.00 |
| Lead G | 0.73 | 0.97 | 1.00 |
| Lead P (round 2 "H") | 1.00 | 0.45 | 1.00 |
| All-ester (round 2 "S") | 1.00 | 0.93 | 0.05 |

## 10. Validation

- OCP FOM1 for PAO 4 at 50 °C: 45.4 against Evonik's published 46.
- Walther viscosity against a commercial PAO-class fluid (EC-120) fitted at 40 and 60 °C: +2.1% at 20 °C and +0.7% at 50 °C. Its published 30 °C point is inconsistent with its neighbors and was treated as an anomaly.
- Heat-capacity temperature slope within 6% of the same fluid's published curve.
- The browser selector (`docs/selector.html`) reproduces the Python models to within 0.01%.
- 26 unit tests (`thermalform/tests`).

## 11. Limits

Screening models, not certification. Flash points are open-cup model values and the PAO flash model is conservative; OCP's base spec uses closed cup. Oxidation stability, materials compatibility and breakdown voltage are not modeled and need lab work. Prices are placeholders.

## 12. Lab plan

1. Ester permittivity at 20–40 GHz (decides the ester ceiling).
2. Closed-cup flash and pour point for the three leads and PAO 2.
3. Thermal conductivity and heat capacity at 25–80 °C.
4. Oxidation (ASTM D2272) and elastomer compatibility for the ester-containing leads.

## References

- Open Compute Project, Base Specification for Immersion Fluids.
- Supplier data sheets: Chevron Phillips Synfluid PAO; Cargill Priolube and NatureCool; Ergon HyVolt; Motul EGEN DCC; M&I MIDEL.
- Evonik, PAO-based immersion fluid figure of merit (published FOM1 for PAO 4).
- Teertstra, Yovanovich and Culham, analytical forced convection modeling of plate-fin heat sinks.
- Muzychka and Yovanovich, pressure drop in laminar developing flow.
- VDI Heat Atlas, chevron plate heat exchanger correlations (Martin).
- Hammerstad and Jensen, accurate models for microstrip computer-aided design.
- US EIA, average retail electricity price, New Jersey, 2024.
