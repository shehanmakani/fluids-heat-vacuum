# Fluids for Heat and Vacuum

**Live site:** https://shehanmakani.github.io/fluids-heat-vacuum/ · **OrbitChem by ChemeNova:** https://chemenova.com/orbitchem

Two research modules built by [Shehan Makani](https://github.com/shehanmakani) for ChemeNova's platforms: **ThermalForm**, an IntelliForm module for immersion cooling fluids, and a **vacuum-materials engine** for OrbitChem.

| Project | Question | Headline result |
| --- | --- | --- |
| **ThermalForm** | Can a PFAS-free fluid meet the OCP immersion-cooling spec and beat the PAO benchmark? | Three lead blends pass at worst-case permittivity. Lead A runs a 1 kW chip 0.40 °C cooler and saves ~$3,100/MW·yr in pumping; Lead G is fully biodegradable and saves ~$5,500. |
| **Space chemistry (OrbitChem)** | Where is the open gap in space chemistry, and can it be modeled from the molecule up? | 54% of 12,842 NASA outgassing tests fail. A metadata model predicts pass/fail on unseen products at AUC 0.81. A physics model of ASTM E595 predicts mass loss within 2.5× across a 100,000-fold volatility range with no fitting. |

## Reports

- [ThermalForm technical report](reports/thermalform-technical-report.md)
- [Space chemistry technical report](reports/space-chemistry-technical-report.md)

## Repository

```
thermalform/           property, heat-sink, loop, dielectric and OCP models; data; results; tests
  evaluate_leads.py    reproduces the lead-blend table from the published recipes
space/lubricants/      vacuum operating-window model for space lubricant base oils
space/outgassing/      NASA outgassing data, statistical model, ASTM E595 physics simulator
reports/               technical reports
docs/                  GitHub Pages site: showcase and the Fluid Selector tool
```

## Run

```bash
pip install -r requirements.txt

cd thermalform && python evaluate_leads.py && python -m pytest -q   # 26 tests
cd ../space/lubricants && python spacelub.py
cd ../outgassing && python profile.py && python model.py && python e595sim.py
```

## What is not here

The blend recipes were found with ThermalForm's multi-objective optimizer, which belongs to ChemeNova's IntelliForm platform and stays private. Everything needed to check the results is public: the data, the property and physics models, and the scripts that evaluate each recipe.

## Caveats

Screening models built from published data, not certification. Supplier prices in `base_stocks.csv` are placeholders. Estimated values are tagged in every data file.

## License

Code: MIT. Reports and figures: CC BY 4.0. NASA outgassing data: public domain (NASA Goddard Space Flight Center).
