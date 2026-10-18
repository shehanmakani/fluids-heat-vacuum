// OrbitChem™ ASTM E595 effusion-cell simulator.
// The E595 specimen compartment is an effusion (Knudsen) cell: vapor leaves through the exit port
// toward a 25 °C collector. Per component: Raoult partial pressure × Hertz–Knudsen flux × effective
// port area; collector capture ETA; re-evaporation from the deposit at 25 °C. 24 h at 125 °C.
// Validated against NASA Goddard outgassing data (TML within 2.5× for five fluids, no fitting).

export interface Component {
  name: string
  pct: number          // mass %
  M: number            // molar mass, kg/mol
  vp25: number         // vapor pressure at 25 °C, Pa
  vp2?: [number, number] // second point [°C, Pa]
  dH?: number          // heat of vaporization, J/mol (used if vp2 absent)
}

export interface ComponentResult { name: string; pct: number; lostPct: number; depositPct: number; P125: number; rule: number }
export interface E595Result {
  tml: number; cvcm: number; nasaPass: boolean; components: ComponentResult[]
  ruleTml: number; ruleCvcm: number
}

const R = 8.314, T_S = 398.15, T_C = 298.15, HOURS = 24
export const GEOM = { m0: 300e-6, dPort: 6.3e-3, aSurf: 1.5e-4, aDep: 5e-4, eta: 0.85 }

function vpFun(c: Component): (T: number) => number {
  const t1 = 298.15
  const b = c.vp2 ? Math.log(c.vp2[1] / c.vp25) / (1 / t1 - 1 / (c.vp2[0] + 273.15)) : (c.dH ?? 1.2e5) / R
  return (T: number) => c.vp25 * Math.exp(b * (1 / t1 - 1 / T))
}

export function aEff(g = GEOM) {
  const ap = Math.PI * g.dPort ** 2 / 4
  return ap * g.aSurf / (ap + g.aSurf)
}

/** Largest P(125 °C)·√M (Pa·(kg/mol)^½) a single low-volatility fluid may have and meet the limit. */
export function passThreshold(limitPct: number, condensable: boolean, g = GEOM) {
  const allowed = limitPct / 100 * g.m0 / (condensable ? g.eta : 1)
  const J = allowed / (aEff(g) * HOURS * 3600)
  return J * Math.sqrt(2 * Math.PI * R * T_S)
}

export function simulate(comps: Component[], g = GEOM): E595Result {
  const n = comps.length
  const vps = comps.map(vpFun)
  const M = comps.map(c => c.M)
  const total = comps.reduce((a, c) => a + c.pct, 0) || 1
  const w0 = comps.map(c => c.pct / total * g.m0)
  const ae = aEff(g)
  const kS = comps.map((c, i) => vps[i](T_S) * Math.sqrt(M[i] / (2 * Math.PI * R * T_S)) * ae)
  const kC = comps.map((c, i) => vps[i](T_C) * Math.sqrt(M[i] / (2 * Math.PI * R * T_C)) * g.aDep)
  const liq = [...w0], dep = new Array(n).fill(0)
  const dt = 20
  for (let s = 0; s < HOURS * 3600 / dt; s++) {
    let ml = 0, md = 0
    for (let i = 0; i < n; i++) { ml += liq[i] / M[i]; md += dep[i] / M[i] }
    for (let i = 0; i < n; i++) {
      const xl = ml > 0 ? liq[i] / M[i] / ml : 0
      const xd = md > 0 ? dep[i] / M[i] / md : 0
      const out = Math.min(kS[i] * xl * dt, liq[i])
      const rev = Math.min(kC[i] * xd * dt, dep[i])
      liq[i] -= out; dep[i] += g.eta * out - rev
    }
  }
  const tml = 100 * (g.m0 - liq.reduce((a, b) => a + b, 0)) / g.m0
  const cvcm = 100 * dep.reduce((a, b) => a + b, 0) / g.m0
  return {
    tml, cvcm, nasaPass: tml <= 1.0 && cvcm <= 0.1,
    ruleTml: passThreshold(1.0, false, g), ruleCvcm: passThreshold(0.1, true, g),
    components: comps.map((c, i) => ({
      name: c.name, pct: c.pct, lostPct: 100 * (w0[i] - liq[i]) / g.m0, depositPct: 100 * dep[i] / g.m0,
      P125: vps[i](T_S), rule: vps[i](T_S) * Math.sqrt(M[i]),
    })),
  }
}

const TORR = 133.322
export interface LibraryFluid extends Omit<Component, 'pct'> { label: string; source: string; nasa?: { tml: number; cvcm: number; n: number } }
export const LIBRARY: LibraryFluid[] = [
  { name: 'PFPE Z25', label: 'PFPE (Fomblin Z25 / Brayco 815Z)', M: 9.5, vp25: 1.2e-11, vp2: [150, 1.1e-5], source: 'NASA CR-2009-215637', nasa: { tml: 0.04, cvcm: 0.01, n: 2 } },
  { name: 'MAC', label: 'MAC (Pennzane 2001A)', M: 0.91, vp25: 3.1e-10, vp2: [150, 2.3e-3], source: 'NASA CR-2009-215637', nasa: { tml: 0.61, cvcm: 0.14, n: 7 } },
  { name: 'Silahydrocarbon', label: 'Penta-silahydrocarbon (NASA)', M: 0.91, vp25: 3.0e-12, dH: 132700, source: 'NASA TM-2001-211196; slope assumed = MAC' },
  { name: 'Target oil', label: 'OrbitChem PFAS-free target oil', M: 0.91, vp25: 1.0e-11, dH: 132700, source: 'Design target: ≤1e-11 Pa at 25 °C' },
  { name: 'Polyphenyl ether', label: 'Polyphenyl ether (Santovac 5)', M: 0.4465, vp25: 4e-10 * TORR, vp2: [275, 0.5 * TORR], source: 'Duniway', nasa: { tml: 46.4, cvcm: 40.0, n: 1 } },
  { name: 'DC 705', label: 'DC 705 silicone pump fluid', M: 0.5469, vp25: 2.4e-9 * TORR, vp2: [243, 0.5 * TORR], source: 'Duniway', nasa: { tml: 97.2, cvcm: 81.0, n: 1 } },
  { name: 'DC 704', label: 'DC 704 silicone pump fluid', M: 0.4848, vp25: 2e-8 * TORR, vp2: [210, 0.5 * TORR], source: 'Duniway', nasa: { tml: 99.85, cvcm: 11.3, n: 1 } },
]
