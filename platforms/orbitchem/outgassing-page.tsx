'use client'
import { useMemo, useState } from 'react'
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Legend } from 'recharts'
import { QualCard, StatusBadge } from '@/components/qual-card'
import { simulate, LIBRARY, type Component } from '@/lib/e595'

const TT = { contentStyle: { background: '#12121A', border: '1px solid #1C2128', borderRadius: 6, fontSize: 11, color: '#F0F2F5' } }
const fmtE = (x: number) => x === 0 ? '0' : x.toExponential(1).replace('e', '×10^')
const pct = (x: number, d = 3) => `${x < 0.001 ? x.toExponential(1) : x.toFixed(d)}%`

const PRESETS: { label: string; comps: Component[] }[] = [
  { label: 'PFPE (flown today)', comps: [{ ...LIBRARY[0], pct: 100 }] },
  { label: 'MAC (Pennzane)', comps: [{ ...LIBRARY[1], pct: 100 }] },
  { label: 'PFAS-free target oil', comps: [{ ...LIBRARY[3], pct: 100 }] },
  { label: 'Base oil + 0.5% volatile tail', comps: [
    { name: 'Distilled base oil', pct: 99.5, M: 0.91, vp25: 3.1e-11, vp2: [150, 2.3e-4] },
    { name: 'Volatile homolog', pct: 0.5, M: 0.63, vp25: 3.1e-9, vp2: [150, 2.3e-2] } ] },
]

export default function OutgassingPage() {
  const [comps, setComps] = useState<Component[]>(PRESETS[3].comps)
  const res = useMemo(() => simulate(comps), [comps])
  const total = comps.reduce((a, c) => a + c.pct, 0)

  const upd = (i: number, k: keyof Component, v: string) =>
    setComps(cs => cs.map((c, j) => j === i ? { ...c, [k]: k === 'name' ? v : Number(v) } : c))
  const add = (name: string) => {
    const f = LIBRARY.find(l => l.name === name); if (!f) return
    setComps(cs => [...cs, { name: f.name, M: f.M, vp25: f.vp25, vp2: f.vp2, dH: f.dH, pct: 10 }])
  }

  const validation = LIBRARY.filter(l => l.nasa).map(l => {
    const r = simulate([{ ...l, pct: 100 }])
    return { ...l, pTml: r.tml, pCvcm: r.cvcm }
  })
  const chart = res.components.map(c => ({ name: c.name, 'Mass loss %': +c.lostPct.toPrecision(3), 'Condensed %': +c.depositPct.toPrecision(3) }))

  return (
    <div className="max-w-7xl mx-auto px-6 py-10">
      <div className="mb-8">
        <div className="section-label">// Outgassing Screen</div>
        <h1 className="section-title">ASTM E595 <em className="text-teal not-italic">Outgassing</em></h1>
        <p className="text-muted text-sm mt-2 max-w-2xl leading-relaxed">
          Physics simulation of the E595 test (24 h, 125 °C, vacuum, 25 °C collector) from each component’s vapor pressure and molar mass.
          The test cell is an effusion cell, so no fitting is needed: predicted mass loss lands within 2.5× of NASA Goddard measurements across a 100,000-fold range of volatility.
        </p>
      </div>

      <div className="grid xl:grid-cols-3 gap-6">
        <div className="space-y-4">
          <div className="card space-y-3">
            <div className="text-[10px] font-mono tracking-widest uppercase text-muted">Presets</div>
            <div className="flex flex-wrap gap-2">
              {PRESETS.map(p => <button key={p.label} className="btn-secondary" onClick={() => setComps(p.comps)}>{p.label}</button>)}
            </div>
          </div>
          <div className="card space-y-3">
            <div className="flex justify-between items-center">
              <div className="text-[10px] font-mono tracking-widest uppercase text-muted">Components</div>
              <div className={`text-[11px] font-mono ${Math.abs(total - 100) < 0.5 ? 'text-success' : 'text-warn'}`}>{total.toFixed(1)}%</div>
            </div>
            {comps.map((c, i) => (
              <div key={i} className="space-y-2 border-t border-white/5 pt-3">
                <div className="flex gap-2">
                  <input className="input" aria-label="Component name" value={c.name} onChange={e => upd(i, 'name', e.target.value)} />
                  <button className="btn-secondary" aria-label={`Remove ${c.name}`} onClick={() => setComps(cs => cs.filter((_, j) => j !== i))}>✕</button>
                </div>
                <div className="grid grid-cols-3 gap-2 text-[10px] font-mono text-muted">
                  <label>Mass %<input className="input mt-1" type="number" value={c.pct} onChange={e => upd(i, 'pct', e.target.value)} /></label>
                  <label>M, kg/mol<input className="input mt-1" type="number" step="0.01" value={c.M} onChange={e => upd(i, 'M', e.target.value)} /></label>
                  <label>P 25 °C, Pa<input className="input mt-1" type="number" value={c.vp25} onChange={e => upd(i, 'vp25', e.target.value)} /></label>
                </div>
              </div>
            ))}
            <select className="select" value="" onChange={e => add(e.target.value)} aria-label="Add a fluid from the library">
              <option value="">+ Add from library…</option>
              {LIBRARY.map(l => <option key={l.name} value={l.name}>{l.label}</option>)}
            </select>
          </div>
        </div>

        <div className="xl:col-span-2 space-y-4">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
            <QualCard label="TML" value={pct(res.tml)} sub="limit 1.0%" status={res.tml <= 1 ? 'pass' : 'fail'} />
            <QualCard label="CVCM" value={pct(res.cvcm)} sub="limit 0.10%" status={res.cvcm <= 0.1 ? 'pass' : 'fail'} />
            <QualCard label="NASA screen" value={res.nasaPass ? 'PASS' : 'FAIL'} status={res.nasaPass ? 'pass' : 'fail'} />
            <QualCard label="Design rule" value="P₁₂₅·√M" sub={`CVCM line ${fmtE(res.ruleCvcm)}`} status="info" />
          </div>

          <div className="card">
            <div className="text-[10px] font-mono tracking-widest uppercase text-muted mb-3">Where the loss comes from</div>
            <div style={{ height: 220 }}>
              <ResponsiveContainer>
                <BarChart data={chart} layout="vertical" margin={{ left: 20, right: 20 }}>
                  <XAxis type="number" stroke="#8A9099" fontSize={10} />
                  <YAxis type="category" dataKey="name" stroke="#8A9099" fontSize={10} width={130} />
                  <Tooltip {...TT} />
                  <Legend wrapperStyle={{ fontSize: 11 }} />
                  <Bar dataKey="Mass loss %" fill="#4A90D9" radius={[0, 4, 4, 0]} />
                  <Bar dataKey="Condensed %" fill="#00C9B1" radius={[0, 4, 4, 0]} />
                </BarChart>
              </ResponsiveContainer>
            </div>
            <table className="w-full text-[11px] font-mono mt-3">
              <thead className="text-muted"><tr><th className="text-left py-1">Component</th><th className="text-right">P at 125 °C</th><th className="text-right">P·√M</th><th className="text-right">vs CVCM line</th></tr></thead>
              <tbody>{res.components.map(c => (
                <tr key={c.name} className="border-t border-white/5">
                  <td className="py-1">{c.name}</td><td className="text-right">{fmtE(c.P125)} Pa</td><td className="text-right">{fmtE(c.rule)}</td>
                  <td className="text-right"><StatusBadge pass={c.rule <= res.ruleCvcm} label={`${(c.rule / res.ruleCvcm).toFixed(1)}×`} /></td>
                </tr>))}</tbody>
            </table>
          </div>

          <div className="card">
            <div className="text-[10px] font-mono tracking-widest uppercase text-muted mb-3">Validation against NASA Goddard data</div>
            <table className="w-full text-[11px] font-mono">
              <thead className="text-muted"><tr><th className="text-left py-1">Fluid</th><th className="text-right">TML pred.</th><th className="text-right">TML NASA</th><th className="text-right">CVCM pred.</th><th className="text-right">CVCM NASA</th></tr></thead>
              <tbody>{validation.map(v => (
                <tr key={v.name} className="border-t border-white/5">
                  <td className="py-1">{v.label}</td><td className="text-right">{pct(v.pTml, 2)}</td><td className="text-right">{pct(v.nasa!.tml, 2)}</td>
                  <td className="text-right">{pct(v.pCvcm, 2)}</td><td className="text-right">{pct(v.nasa!.cvcm, 2)}</td>
                </tr>))}</tbody>
            </table>
            <p className="text-[11px] text-muted mt-3 leading-relaxed">CVCM runs high for MAC and DC 704, so the screen errs on the safe side. Liquids only; cured polymers need diffusion data. Exit port 6.3 mm, collector capture 0.85. Screening estimate, not certification.</p>
          </div>
        </div>
      </div>
    </div>
  )
}
