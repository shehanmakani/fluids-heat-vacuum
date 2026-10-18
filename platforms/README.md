# Running the work on ChemeNova's platforms

## IntelliForm v1.6 · ThermalForm vertical

`intelliform/thermalform_vertical.py` is a new IntelliForm module (branch `feature/thermalform-vertical` of [Cheme-Nova/IntelliForm](https://github.com/Cheme-Nova/IntelliForm)). It runs IntelliForm's NSGA-III Pareto engine, with its Das-Dennis reference directions and sparse seeding, using ThermalForm physics as the objectives. It then picks a blend with IntelliForm's own TOPSIS selector. IntelliForm's three axes map as follows:

| IntelliForm axis | ThermalForm meaning |
| --- | --- |
| cost | facility pumping cost, USD per MW-year |
| bio % | readily biodegradable share |
| performance | chip cooling score, 100 × R<sub>benchmark</sub> / R<sub>blend</sub> |

Constraints follow the OCP spec at worst-case ester permittivity. Mineral oil is excluded, and a blend may use at most 3 stocks, each ≥ 5%.

**Results** (`intelliform_pareto.json`, `run_intelliform.py`):

| Run | Pareto blends | Blend | Chip ΔT vs PAO | Pumping vs PAO, $/MW·yr | Worst Dk | Biodegradable |
| --- | --- | --- | --- | --- | --- | --- |
| Flash ≥ 150 °C, TOPSIS pick | 35 | renewable isoparaffin 68.7% + EF 7010 31.3% | +0.29 °C | −7,213 | 2.300 | 100% |
| Flash ≥ 150 °C, coolest 100%-biodegradable | | renewable isoparaffin 62.9% + Priolube 1976 37.1% | −0.22 °C | −5,534 | 2.288 | 100% |
| Flash ≥ 150 °C, coolest overall | | PAO 4 50.4% + Priolube 1976 37.2% + PAO 6 12.3% | −0.87 °C | +3,469 | 2.300 | 37% |
| Flash ≥ 210 °C, TOPSIS pick | 12 | PAO 4 54.9% + Priolube 1976 37.5% + renewable isoparaffin 7.7% | −0.76 °C | +686 | 2.300 | 45% |

- **IntelliForm independently recovered the round 3 leads.** Its coolest 100%-biodegradable blend is Lead G to within one percentage point. The Lead A family (PAO 4 + Priolube 1976 + PAO 2) also sits on its front.
- **Several picks sit exactly on the Dk 2.3 limit.** TOPSIS gives pumping cost a 40% weight, so it favors blends that ride the constraint. Ester permittivity at GHz is the unmeasured input, so the next step is a TOPSIS run with a safety margin (for example Dk ≤ 2.25).

**QSAR check** (`intelliform_qsar_twin.json`): IntelliForm's Morgan-fingerprint QSAR tier (trained on 398 ingredients) was run on representative structures for each base stock.

- It ranks esters above hydrocarbons correctly, at 79–89 against 74–76.
- It scores every hydrocarbon near 75, so it calls PAOs readily biodegradable. Only 5 of 9 stocks agree with data sheets.

Base oils are outside its training domain. Feeding base-oil biodegradation data into IntelliForm's feedback loop would fix that, and this run is the evidence that the data is needed.

**Digital twin:** a 1-tonne pilot batch of each lead fluid takes 2.6 h per cycle. That works out to 6.6 batches a day, 80 kWh per batch and 0.031 kg CO₂ per kg. All of these are simple blends that ChemRich USA can make.

## OrbitChem · ASTM E595 outgassing screen

OrbitChem's web app listed its standalone outgassing screen as "coming next sprint." Branch `feature/e595-outgassing-screen` of [Cheme-Nova/orbitchem-web](https://github.com/Cheme-Nova/orbitchem-web) builds it:

- `orbitchem/e595.ts`: the effusion-cell simulator, ported to TypeScript and matching the Python model.
- `orbitchem/outgassing-page.tsx`: the `/outgassing` page, using OrbitChem's own components and styling.

The page runs in the browser, so it works without the backend. It includes presets, a component editor, a per-component breakdown, the P<sub>125</sub>·√M design rule, and the NASA validation table. The app builds cleanly (`next build`), and the page reproduces the Python results (`outgassing-screen.png`).

**Not yet integrated:** the OrbitChem FastAPI backend (`Cheme-Nova/OrbitChem`) was not reachable from this work, and its hosted API was down. The radiation, thermal-cycling and contamination modules were therefore not run on these materials. OrbitChem's published ML notebooks (`cheminformatics-ml` 06–07) model space solar semiconductors and do not apply to lubricants or outgassing.
