# Space chemistry with OrbitChem: vacuum lubricants and outgassing

*Vacuum-materials engine for OrbitChem*

**Shehan Makani · ChemeNova LLC · October 2026**

## Abstract

A scan of eight space-chemistry niches, scored for market pull, timing, fit, capital, qualification burden, competition and story crowding, points to two linked opportunities: PFAS-free vacuum lubricants and low-outgassing materials. Both depend on the same molecular properties, vapor pressure and molar mass. Three analyses follow. A vacuum operating-window model shows the best non-fluorinated base oil beats PFPE at the hot end but falls 16 °C short at the cold end. A statistical study of NASA Goddard's 12,842-test outgassing database finds that 54% of materials fail ASTM E595 and that a metadata model predicts pass or fail on unseen products with AUC 0.81. A physics simulator of the E595 test predicts total mass loss within a factor of 2.5 across a more than 100,000-fold range of volatility with no fitted parameters, and yields a design rule usable for inverse molecular design.

## 1. Opportunity scan

Each niche was scored 1–5 on six criteria, with market pull, timing and fit weighted double (maximum 45), then re-scored with a seventh criterion for story crowding (5 = quiet).

| Niche | Score /45 | Crowding | Total /50 |
| --- | --- | --- | --- |
| Low-outgassing materials: prediction and fixing | 34 | 5 | 39 |
| PFAS-free space lubricants | 35 | 2 | 37 |
| Radiation-tolerant polymers | 26 | 3 | 29 |
| Thermal control and orbital data centers | 24 | 2 | 26 |
| Academic astrochemistry | 22 | 4 | 26 |
| Green propellants | 20 | 3 | 23 |
| In-space pharmaceutical manufacturing | 20 | 1 | 21 |
| In-situ resource use | 15 | 2 | 17 |

Scores are the author's judgment, informed by the sources listed.

## 2. PFAS-free vacuum lubricants

Space mechanisms rely on perfluoropolyether (PFPE) oils and greases, which are PFAS. ESA's 2025 ionic-liquid project still needed a PTFE thickener and reported that no fully PFAS-free space grease exists. Non-fluorinated candidates include multiply alkylated cyclopentane (MAC, Pennzane), ionic liquids and silahydrocarbons.

**Window model** (`space/lubricants/spacelub.py`):

- Cold limit: temperature where kinematic viscosity reaches 10,000 cSt, from a Walther fit through the two coldest published viscosities.
- Hot limit: temperature where free-molecular (Langmuir) evaporation from 1 cm² removes 5 mg over 15 years, from a Clausius–Clapeyron fit through published vapor pressures.
- Blends: Refutas viscosity mixing; Raoult's law for vapor pressure.

| Base oil | Cold limit | Hot limit | Width |
| --- | --- | --- | --- |
| PFPE (Brayco 815Z / Fomblin Z25) | −45.7 °C | 68.7 °C | 114 K |
| Penta-silahydrocarbon (NASA) | −29.9 °C | 81.4 °C | 111 K |
| MAC + 50% silahydrocarbon | −25.6 °C | 52.5 °C | 78 K |
| MAC (Pennzane 2001A) | −21.4 °C | 48.0 °C | 69 K |

The ordering holds when the cold limit varies from 5,000 to 20,000 cSt and the allowed loss varies 100-fold. The silahydrocarbon's vapor-pressure slope and molar mass are assumed equal to MAC's. **Inverse-design target:** a non-fluorinated oil with ≤ 10,000 cSt at −40 °C and ≤ 1×10⁻¹¹ Pa at 25 °C.

## 3. NASA outgassing database

NASA Goddard's database (`space/outgassing/data/nasa_outgassing_db.csv`, public domain) holds ASTM E595 results: total mass loss (TML) and collected volatile condensable material (CVCM) after 24 h at 125 °C in vacuum, with a 25 °C collector. Pass: TML ≤ 1.0% and CVCM ≤ 0.10%. After removing sentinel values, 12,842 tests from 1,476 manufacturers remain (`profile.py`).

- 46% pass. 23% fail both limits, 20% TML only, 11% CVCM only.
- 37% of failures lie within twice both limits.
- Devolatilized or vacuum-baked materials pass 69% of the time against 45% otherwise; silicones 68% against 41%. In 12 product families tested both ways, the pass rate rose from 1% to 59%.
- With the water-corrected criterion (TML − WVR ≤ 1.0%), the pass rate rises from 48% to 57%; for nylons and aramids from 6% to 44%.
- Of 414 families tested at least three times, 40% received both verdicts; a single test agrees with its family majority 86% of the time.
- Pass rates rose from 41% in the earliest fifth of tests to 53% in the latest.

Chemistry classes were assigned by keyword matching on product names (67% of rows assignable):

| Class | Tests | Pass rate | Median TML % | Median CVCM % |
| --- | --- | --- | --- | --- |
| PEEK / PEI / PPS / polysulfone | 71 | 96% | 0.31 | 0.00 |
| Fluoropolymers | 454 | 88% | 0.12 | 0.01 |
| PFPE oils and greases | 93 | 74% | 0.14 | 0.04 |
| Polyester / PET | 384 | 57% | 0.48 | 0.06 |
| Epoxy | 2,758 | 50% | 0.86 | 0.03 |
| Polyimide | 418 | 50% | 0.97 | 0.01 |
| Silicone | 2,026 | 44% | 0.48 | 0.13 |
| Polyurethane | 889 | 42% | 1.06 | 0.05 |
| Polyolefin | 191 | 36% | 0.58 | 0.13 |
| Acrylic | 543 | 34% | 1.33 | 0.05 |
| Rubber (EPDM, NBR, neoprene) | 191 | 26% | 1.35 | 0.30 |
| Hydrocarbon lubricants (MAC, PAO) | 75 | 20% | 1.42 | 0.34 |
| Polyamide / aramid | 281 | 8% | 2.43 | 0.04 |

## 4. Statistical pass/fail model

Logistic regression on one-hot class, usage, manufacturer, cure and process flags, plus character n-grams of the product name (`model.py`). Five-fold cross-validation grouped by product family (lot and batch tokens stripped) and, separately, by brand.

| Inputs | Family AUC | Family accuracy | Brand AUC | Brand accuracy |
| --- | --- | --- | --- | --- |
| Class + usage | 0.72 | 66% | 0.69 | 64% |
| + maker, cure, process flags | 0.76 | 69% | 0.71 | 65% |
| + name text | 0.81 | 74% | 0.74 | 67% |

Predicted probabilities are calibrated within about five points per decile. Testing only the top 30% raises the pass rate from 46% to 79%. The ceiling set by replicate disagreement is about 86%. The database records trade names, not compositions, which limits any metadata model.

## 5. Physics simulator of ASTM E595

The E595 specimen compartment is an effusion cell. `e595sim.py` integrates, for each component, effusion through the exit port (Raoult partial pressure × Hertz–Knudsen flux × effective area combining sample surface and a 6.3 mm port), capture of a fraction η = 0.85 of the effused vapor on the collector, and re-evaporation from the 25 °C deposit. Nothing is fitted to TML; η is set from the two diffusion-pump fluids.

| Fluid | P at 125 °C (Pa) | TML pred. | TML NASA | CVCM pred. | CVCM NASA |
| --- | --- | --- | --- | --- | --- |
| DC 704 | 0.19 | 100% | 99.9% | 85% | 11.3% |
| DC 705 | 0.028 | 100% | 97.2% | 85% | 81.0% |
| Polyphenyl ether (5P4E class)* | 0.0055 | 19% | 46% | 16% | 40% |
| MAC (7 lots) | 2.2×10⁻⁴ | 1.06% | 0.61% | 0.90% | 0.14% |
| PFPE Z25 / 815Z (2 lots) | 1.4×10⁻⁶ | 0.022% | 0.04% | 0.019% | 0.01% |

*Santovac 5 vapor pressure; NASA's entry is the related Convalex-10.

TML errors are within a factor of 2.5 throughout; plausible geometry changes move absolute values by up to about 2× without changing the order. CVCM matches three fluids and is overpredicted for MAC and DC 704. At DC 704's vapor pressure, a 25 °C deposit cannot re-evaporate more than about 0.2% of the sample in 24 h, so that entry likely reflects creep or a test artifact.

**Design rule.** For a single low-volatility fluid in the linear-loss regime:

- TML ≤ 1%: P<sub>125 °C</sub>·√M < 1.9×10⁻⁴ Pa·(kg/mol)<sup>½</sup>
- CVCM ≤ 0.10%: P<sub>125 °C</sub>·√M < 2.3×10⁻⁵ Pa·(kg/mol)<sup>½</sup> (binding)

The penta-silahydrocarbon is predicted at TML 0.010% / CVCM 0.009% and the inverse-design target oil at 0.034% / 0.029%, so the lubricant target is also outgassing-safe. An illustrative two-component case shows a 0.5% volatile tail failing an otherwise passing oil, and a 24 h, 150 °C thin-film vacuum strip restoring a pass.

## 6. Limits

Liquids only; cured polymers need diffusion coefficients per component. The exit-port size should be checked against the standard. Water is not modeled. Scan scores are judgment. Market sizes and program details are as published at the time of writing.

## 7. Next steps

Composition-aware outgassing model trained on 20–40 formulations with known recipes and E595 results; diffusion-limited extension for adhesives and potting compounds; OrbitChem inverse design against the base-oil target.

## References

- NASA Goddard Space Flight Center, Spacecraft Material Outgassing Data (data.nasa.gov).
- ASTM E595, Total Mass Loss and Collected Volatile Condensable Materials from Outgassing in a Vacuum Environment.
- NASA CR-2009-215637, MAC and PFPE bearing lifetimes and vapor pressures.
- NASA TM-2001-211196, Silahydrocarbons for space mechanisms.
- ESA executive summary, ionic-liquid vacuum lubricants (Materiales, ESTL), 2025.
- ESR Technology / ESTL, impact of PFAS restriction on spacecraft lubricants; ESMATS 2025 guidelines for new fluid lubricants.
- Eurospace, comments to ECHA on the universal PFAS restriction, May 2026.
- Duniway Stockroom and Neyco, diffusion pump fluid vapor pressures.
- Oneida Research Services, ASTM E595 testing, CMSE 2024.
