import { useState, useMemo, useCallback, useRef } from "react";
import impuritiesData from "./impuritiesTable.json";

const COLORS = ["#22c55e","#f59e0b","#ef4444","#a78bfa","#ec4899","#06b6d4","#f97316","#84cc16"];
const GROUP_COLORS = ["#22c55e","#f59e0b","#a78bfa","#ec4899","#06b6d4","#f97316","#84cc16","#14b8a6","#e879f9","#fb923c","#64748b","#facc15"];

const SOLVENT_LABELS = {
  "THF-d8": "THF-d₈", "CD2Cl2": "CD₂Cl₂", "CDCl3": "CDCl₃", "toluene-d8": "Toluene-d₈",
  "C6D6": "C₆D₆", "C6D5Cl": "C₆D₅Cl", "(CD3)2CO": "Acetone-d₆", "(CD3)2SO": "DMSO-d₆",
  "CD3CN": "CD₃CN", "TFE-d3": "TFE-d₃", "CD3OD": "CD₃OD", "D2O": "D₂O",
};

const SOLVENT_CLASS = {
  "CDCl3": { cls: "I", desc: "Weak aprotic" }, "CD2Cl2": { cls: "I", desc: "Weak aprotic" },
  "C6D6": { cls: "I", desc: "Weak aprotic (aromatic)" }, "C6D5Cl": { cls: "I", desc: "Weak aprotic (aromatic)" },
  "toluene-d8": { cls: "I", desc: "Weak aprotic (aromatic)" },
  "(CD3)2SO": { cls: "II", desc: "Strong H-bond acceptor" }, "(CD3)2CO": { cls: "II", desc: "Strong H-bond acceptor" },
  "CD3CN": { cls: "II", desc: "Strong H-bond acceptor" }, "THF-d8": { cls: "II", desc: "Strong H-bond acceptor" },
  "CD3OD": { cls: "III", desc: "Protic / exchanging" }, "D2O": { cls: "III", desc: "Protic / exchanging" },
  "TFE-d3": { cls: "III", desc: "Protic / exchanging" },
};

const COMPOUND_ALIASES = {
  "dimethylformamide": ["dmf"], "dichloromethane": ["dcm", "ch2cl2"], "chloroform": ["chcl3", "hccl3"],
  "tetrahydrofuran": ["thf"], "methanol": ["meoh"], "diethyl ether": ["et2o", "ether"],
  "ethyl acetate": ["etoac", "ea"], "H grease": ["grease"], "silicone grease": ["silicone"],
  "n-hexane": ["hexane"], "n-pentane": ["pentane"], "triethylamine": ["tea", "et3n"],
  "pyridine": ["py"], "1,4-dioxane": ["dioxane"], "acetonitrile": ["mecn", "ch3cn"],
  "ethanol": ["etoh"], "2-propanol": ["ipa", "isopropanol"], "toluene": ["phme"],
  "dimethyl sulfoxide": ["dmso"],
};

const NMR_SOLV_NAMES = {
  "CDCl3": ["cdcl3", "chloroform-d"], "CD2Cl2": ["cd2cl2", "dcm"],
  "(CD3)2SO": ["(cd3)2so", "dmso", "dmso-d6"], "(CD3)2CO": ["(cd3)2co", "acetone", "acetone-d6"],
  "CD3CN": ["cd3cn", "mecn", "acetonitrile"], "CD3OD": ["cd3od", "meod", "methanol"],
  "D2O": ["d2o", "water"], "C6D6": ["c6d6", "benzene"], "C6D5Cl": ["c6d5cl", "chlorobenzene", "chlorobenzene-d5"],
  "THF-d8": ["thf-d8", "thf"], "toluene-d8": ["toluene-d8", "toluene"], "TFE-d3": ["tfe-d3"],
};

// Build SOLVENT_DB from impuritiesTable.json
const SOLVENT_DB = {};
for (const solv of impuritiesData.solvents) {
  const impurities = [];
  for (const [compName, groups] of Object.entries(impuritiesData.compounds)) {
    const ppmValues = [];
    for (const shifts of Object.values(groups)) {
      if (shifts[solv] != null) ppmValues.push(shifts[solv]);
    }
    if (ppmValues.length > 0) impurities.push({ name: compName, ppm: ppmValues });
  }
  SOLVENT_DB[solv] = {
    label: SOLVENT_LABELS[solv] || solv,
    solventPeaks: impuritiesData.solvent_signals[solv] || [],
    impurities,
  };
}

function parseMNOVA(text) {
  const lines = text.trim().split("\n").filter(l => l.trim());
  const peaks = [];
  for (const line of lines) {
    const parts = line.split("\t");
    if (parts.length < 5 || !/^\d+$/.test(parts[0].trim())) continue;
    const ppm = parseFloat(parts[1]), int_ = parseFloat(parts[2]), w = parseFloat(parts[3]), area = parseFloat(parts[4]);
    const type = (parts[5] || "").trim(), flags = (parts[6] || "").trim(), compound = (parts[7] || "").trim(), annotation = (parts[8] || "").trim();
    if (isNaN(ppm)) continue;
    peaks.push({ ppm, intensity: int_, width: isNaN(w) ? 1.5 : w, area: isNaN(area) ? 0 : area, type, flags, compound, annotation });
  }
  return peaks;
}

function parseChemDraw(text) {
  const lines = text.split("\n"); const nodes = []; let cur = null;
  for (const line of lines) {
    const m = line.match(/^(CH3|CH2|CH|C)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(.*)/);
    if (m) { if (cur) nodes.push(cur); cur = { type: m[1], shift: parseFloat(m[2]), desc: m[4].trim(), increments: [m[4].trim()] }; }
    else if (cur && line.match(/^\s+[-\d?.]+\s+/)) cur.increments.push(line.trim());
  }
  if (cur) nodes.push(cur);
  return nodes;
}

function lorentzian(x, x0, g, h) { return h * g * g / ((x - x0) * (x - x0) + g * g); }

function genSpectrum(peaks, min, max, n, gamma) {
  const pts = [], step = (max - min) / n;
  for (let i = 0; i <= n; i++) { const x = min + i * step; let y = 0; for (const p of peaks) y += lorentzian(x, p.ppm, gamma, p.intensity); pts.push([x, y]); }
  return pts;
}

function evalAt(x, peaks, gamma) { let y = 0; for (const p of peaks) y += lorentzian(x, p.ppm, gamma, p.intensity); return y; }

// Per-compound greedy match. Each compound matches independently.
// Within one compound, no actual peak reused for 2 predicted peaks.
// For major product: sort candidates by area (desc) first, then delta.
function matchSingleCompound(predicted, actual, tol, isProduct) {
  const candidates = [];
  for (let pi = 0; pi < predicted.length; pi++) {
    for (const a of actual) {
      const d = Math.abs(a.ppm - predicted[pi].shift);
      if (d <= tol) candidates.push({ pi, pred: predicted[pi], act: a, delta: d, area: a.area });
    }
  }
  // Major product: within tolerance, prefer highest area first, break ties by delta.
  // Non-major: prefer closest delta, area irrelevant.
  if (isProduct) {
    candidates.sort((a, b) => {
      if (Math.abs(b.area - a.area) > 100) return b.area - a.area;
      return a.delta - b.delta;
    });
  } else {
    candidates.sort((a, b) => a.delta - b.delta);
  }

  const usedPred = new Set(), usedAct = new Set();
  const results = [];
  for (const c of candidates) {
    if (usedPred.has(c.pi) || usedAct.has(c.act.ppm)) continue;
    usedPred.add(c.pi); usedAct.add(c.act.ppm);
    results.push({ ...c.pred, match: c.act, delta: c.delta, matched: true, _pi: c.pi });
  }
  for (let pi = 0; pi < predicted.length; pi++) {
    if (usedPred.has(pi)) continue;
    let best = null, bestD = Infinity;
    for (const a of actual) { const d = Math.abs(a.ppm - predicted[pi].shift); if (d < bestD) { bestD = d; best = a; } }
    results.push({ ...predicted[pi], match: best, delta: bestD, matched: false, _pi: pi });
  }
  results.sort((a, b) => a._pi - b._pi);
  return results;
}

function scoreMatch(m, maxArea, isProduct) {
  if (!m.matched || !m.match) return 0;
  const distScore = Math.max(0, 1 - m.delta / 6);
  const areaRatio = m.match.area / maxArea;
  if (isProduct) {
    // Major product: 50% proximity, 50% area (expect dominant peaks)
    const areaScore = Math.min(1, areaRatio * 2.5);
    return distScore * 0.5 + areaScore * 0.5;
  }
  // Contaminant: 85% proximity, 15% area (small peaks expected)
  return distScore * 0.85 + Math.min(1, areaRatio * 3) * 0.15;
}

function groupPeaks(peaks, gap = 2.0) {
  if (!peaks.length) return [];
  const s = [...peaks].sort((a, b) => a.ppm - b.ppm);
  const gs = []; let g = [s[0]];
  for (let i = 1; i < s.length; i++) { if (s[i].ppm - g[g.length - 1].ppm <= gap) g.push(s[i]); else { gs.push(g); g = [s[i]]; } }
  gs.push(g);
  return gs.map(gr => ({ center: gr.reduce((s, p) => s + p.ppm * p.intensity, 0) / gr.reduce((s, p) => s + p.intensity, 0), totalArea: gr.reduce((s, p) => s + p.area, 0), peaks: gr, min: gr[0].ppm, max: gr[gr.length - 1].ppm }));
}

const IMP_TOL = 2; // ppm window for matching impurity reference shifts to actual peaks

const TABS = ["input", "spectrum", "analysis"];

export default function NMRTool() {
  const [tab, setTab] = useState("input");
  const [mnovaText, setMnovaText] = useState("");
  const [compounds, setCompounds] = useState([{ name: "Product", text: "", color: COLORS[0], isProduct: true }]);
  const [tol, setTol] = useState(6);
  const [fullRange, setFullRange] = useState([0, 220]);
  const [viewRange, setViewRange] = useState([0, 220]);
  const [linewidth, setLinewidth] = useState(0.05);
  const [showSolvent, setShowSolvent] = useState(true);
  const [scaleToComp, setScaleToComp] = useState(true);
  const [solvent, setSolvent] = useState("CDCl3");
  const [enabledImp, setEnabledImp] = useState({});
  const [predSolvent, setPredSolvent] = useState("none");
  const [visibleComps, setVisibleComps] = useState({});

  const [zoomMode, setZoomMode] = useState(false);
  const [dragging, setDragging] = useState(false);
  const [dragStart, setDragStart] = useState(null);
  const [dragEnd, setDragEnd] = useState(null);
  const svgRef = useRef(null);
  const [zoomStack, setZoomStack] = useState([]);

  const solvData = SOLVENT_DB[solvent] || SOLVENT_DB["CDCl3"];
  const allPeaks = useMemo(() => parseMNOVA(mnovaText), [mnovaText]);
  const compPeaks = useMemo(() => allPeaks.filter(p => p.type !== "Solvent"), [allPeaks]);
  const solvPeaks = useMemo(() => allPeaks.filter(p => p.type === "Solvent"), [allPeaks]);
  const nmrSolvNames = useMemo(() => (NMR_SOLV_NAMES[solvent] || [solvent.toLowerCase()]), [solvent]);
  const nmrSolvPeaks = useMemo(() => solvPeaks.filter(p => {
    if (!p.compound) return true;
    const c = p.compound.toLowerCase();
    return nmrSolvNames.some(n => n === c || c.includes(n));
  }), [solvPeaks, nmrSolvNames]);
  const impSolvPeaks = useMemo(() => solvPeaks.filter(p => {
    if (!p.compound) return false;
    const c = p.compound.toLowerCase();
    return !nmrSolvNames.some(n => n === c || c.includes(n));
  }), [solvPeaks, nmrSolvNames]);
  const compData = useMemo(() => compounds.map(c => ({ ...c, predicted: parseChemDraw(c.text) })), [compounds]);
  const maxArea = useMemo(() => Math.max(...compPeaks.map(p => p.area), 1), [compPeaks]);

  // ── Assignment pipeline: major → impurities → minor ──────────────────────
  // Peaks are NOT exclusively claimed — the same peak can be matched by multiple
  // matchers. Priority order only determines what counts as "assigned" in the
  // unassigned display: a peak is assigned if ANY matcher claims it.

  // Step 1: match major products against all compound peaks
  const majorMatchRes = useMemo(() => {
    return compData
      .map((c, i) => ({ ...c, origIdx: i }))
      .filter(c => c.isProduct)
      .map(c => ({ ...c, matches: matchSingleCompound(c.predicted, compPeaks, tol, true) }));
  }, [compData, compPeaks, tol]);

  // Step 2: resolve impurities against all compound peaks (no exclusive claiming).
  // Each entry has { refPpm, matchedPeak, hasMatch, isAnnotated }.
  // Annotated entries (found in MNOVA compound/annotation field) override table
  // matching and use the exact annotated ppm — takes priority over everything.
  const resolvedImpurities = useMemo(() => {
    return solvData.impurities.filter(ki => enabledImp[ki.name]).map(ki => {
      const keys = [ki.name.toLowerCase(), ...(COMPOUND_ALIASES[ki.name] || [])];
      const annotated = allPeaks.filter(p => {
        const fields = [p.annotation, p.type, p.compound, p.flags].filter(Boolean).map(s => s.toLowerCase());
        return fields.some(f => keys.some(k => f.includes(k)));
      });
      if (annotated.length > 0) {
        return {
          name: ki.name, refPpm: ki.ppm, fromAnnotation: true,
          entries: annotated.map(p => ({ refPpm: p.ppm, matchedPeak: p, hasMatch: true, isAnnotated: true })),
        };
      }
      const entries = ki.ppm.map(refPpm => {
        let best = null, bestD = Infinity;
        for (const p of compPeaks) {
          const d = Math.abs(p.ppm - refPpm);
          if (d <= IMP_TOL && d < bestD) { best = p; bestD = d; }
        }
        return { refPpm, matchedPeak: best, hasMatch: best !== null, isAnnotated: false };
      });
      return { name: ki.name, refPpm: ki.ppm, fromAnnotation: false, entries };
    });
  }, [solvData, enabledImp, allPeaks, compPeaks]);

  // Step 3: match minor products against all compound peaks (no exclusive claiming)
  const minorMatchRes = useMemo(() => {
    return compData
      .map((c, i) => ({ ...c, origIdx: i }))
      .filter(c => !c.isProduct)
      .map(c => ({ ...c, matches: matchSingleCompound(c.predicted, compPeaks, tol, false) }));
  }, [compData, compPeaks, tol]);

  // Combined for display, original order preserved
  const matchRes = useMemo(() =>
    [...majorMatchRes, ...minorMatchRes].sort((a, b) => a.origIdx - b.origIdx),
    [majorMatchRes, minorMatchRes]);

  // A peak is "assigned" if claimed by any matcher (union, priority: major → imp → minor)
  const allAssigned = useMemo(() => {
    const a = new Set();
    majorMatchRes.forEach(c => c.matches.filter(m => m.matched && m.match).forEach(m => a.add(m.match.ppm)));
    resolvedImpurities.forEach(ki => ki.entries.forEach(e => { if (e.matchedPeak) a.add(e.matchedPeak.ppm); }));
    minorMatchRes.forEach(c => c.matches.filter(m => m.matched && m.match).forEach(m => a.add(m.match.ppm)));
    return a;
  }, [majorMatchRes, resolvedImpurities, minorMatchRes]);

  // ── Spectrum curves ───────────────────────────────────────────────────────
  const ppmSpan = viewRange[1] - viewRange[0];
  const nPts = Math.min(8000, Math.max(3000, Math.round(60000 / Math.max(ppmSpan, 1))));
  const compCurve = useMemo(() => compPeaks.length ? genSpectrum(compPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [compPeaks, viewRange, linewidth, nPts]);
  const nmrSolvCurve = useMemo(() => nmrSolvPeaks.length ? genSpectrum(nmrSolvPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [nmrSolvPeaks, viewRange, linewidth, nPts]);
  const impSolvCurve = useMemo(() => impSolvPeaks.length ? genSpectrum(impSolvPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [impSolvPeaks, viewRange, linewidth, nPts]);
  const compMaxY = useMemo(() => { const v = compCurve.map(p => p[1]); return v.length ? Math.max(...v, 1) : 1; }, [compCurve]);
  const solvMaxY = Math.max(...nmrSolvCurve.map(p => p[1]), ...impSolvCurve.map(p => p[1]), 0);
  const displayMaxY = scaleToComp ? compMaxY * 1.12 : Math.max(compMaxY, solvMaxY, 1) * 1.1;

  const unassigned = useMemo(() => compPeaks.filter(p => !allAssigned.has(p.ppm)), [compPeaks, allAssigned]);
  const groups = useMemo(() => groupPeaks(compPeaks), [compPeaks]);

  const solventPrediction = useMemo(() => {
    if (predSolvent === "none" || predSolvent === solvent) return null;
    const predCls = SOLVENT_CLASS[predSolvent];
    const actCls = SOLVENT_CLASS[solvent];
    const carbonyl = [], aromatic = [], aliphatic = [];
    for (const groups of Object.values(impuritiesData.compounds)) {
      for (const [groupName, shifts] of Object.entries(groups)) {
        const pS = shifts[predSolvent], aS = shifts[solvent];
        if (pS == null || aS == null) continue;
        const d = aS - pS;
        if (groupName.includes("CO") && pS > 150) carbonyl.push(d);
        else if (pS > 100) aromatic.push(d);
        else aliphatic.push(d);
      }
    }
    const summarize = arr => {
      if (!arr.length) return null;
      arr.sort((a, b) => a - b);
      return { min: arr[0], max: arr[arr.length - 1], avg: arr.reduce((s, v) => s + v, 0) / arr.length, n: arr.length };
    };
    const guidance = [];
    if (predCls && actCls) {
      const f = predCls.cls, t = actCls.cls;
      if (f === t) guidance.push("Same solvent class -- expect small shifts (<1 ppm for most carbons).");
      else if (f === "I" && t === "II") guidance.push("Class I->II: Polar carbons (C=O, C-O, C-N) shift most. Increase tolerance for carbons near H-bond donors.");
      else if (f === "I" && t === "III") guidance.push("Class I->III: Protic solvent -- carbonyl and polar carbons shift significantly. Consider +2-4 ppm tolerance.");
      else if (f === "II" && t === "I") guidance.push("Class II->I: Reverse polarity change. Polar carbons shift upfield. Increase tolerance for functional groups.");
      else if (f === "II" && t === "III") guidance.push("Class II->III: Protic vs acceptor -- moderate changes for polar carbons.");
      else if (f === "III" && t === "I") guidance.push("Class III->I: Large environment change. Polar functional groups shift most. Consider wider tolerance.");
      else if (f === "III" && t === "II") guidance.push("Class III->II: Moderate shift expected for polar carbons.");
      if (predCls.desc.includes("aromatic") || actCls.desc.includes("aromatic"))
        guidance.push("Aromatic solvent involved -- ring-current anisotropy may cause additional shifts, especially for sp2 carbons.");
    }
    return { predCls, actCls, carbonyl: summarize(carbonyl), aromatic: summarize(aromatic), aliphatic: summarize(aliphatic), guidance };
  }, [predSolvent, solvent]);

  const addComp = () => setCompounds(c => [...c, { name: `Compound ${c.length + 1}`, text: "", color: COLORS[c.length % COLORS.length], isProduct: false }]);
  const rmComp = i => setCompounds(c => c.filter((_, j) => j !== i));
  const updComp = (i, k, v) => setCompounds(c => c.map((x, j) => j === i ? { ...x, [k]: v } : x));
  const isCompVisible = (name) => visibleComps[name] !== false;
  const toggleCompVis = (name) => setVisibleComps(v => ({ ...v, [name]: v[name] === false ? true : false }));

  const numC = compData.filter(c => c.predicted.length > 0).length;
  const mzH = Math.max(50, numC * 24 + 20);
  const svgW = 900, svgH = 360 + mzH;
  const maxLabelWidth = useMemo(() => {
    const names = compData.filter(c => c.predicted.length > 0).map(c => c.name);
    if (!names.length) return 20;
    const longest = Math.max(...names.map(n => n.length));
    return Math.min(63, longest * 6.5 + 8); // ~6.5px per char at fontSize 10, plus padding
  }, [compData]);

  const pad = { t: mzH + 10, b: 50, l: maxLabelWidth, r: 20 };
  const plotW = svgW - pad.l - pad.r, plotH = svgH - pad.t - pad.b;
  const xSc = useCallback(ppm => pad.l + plotW * (1 - (ppm - viewRange[0]) / (viewRange[1] - viewRange[0])), [viewRange, plotW, pad.l]);
  const ySc = useCallback(v => pad.t + plotH * (1 - Math.min(v, displayMaxY) / displayMaxY), [displayMaxY, plotH, pad.t]);
  const xInv = useCallback(px => viewRange[1] - (px - pad.l) / plotW * (viewRange[1] - viewRange[0]), [viewRange, plotW, pad.l]);
  const makePath = pts => pts.map((p, i) => `${i === 0 ? "M" : "L"}${xSc(p[0]).toFixed(1)},${ySc(p[1]).toFixed(1)}`).join(" ");

  const svgToPpm = useCallback((e) => { if (!svgRef.current) return 0; const r = svgRef.current.getBoundingClientRect(); return xInv((e.clientX - r.left) / r.width * svgW); }, [xInv, svgW]);
  const handleMouseDown = useCallback((e) => { if (!zoomMode) return; const p = svgToPpm(e); setDragging(true); setDragStart(p); setDragEnd(p); }, [zoomMode, svgToPpm]);
  const handleMouseMove = useCallback((e) => { if (dragging) setDragEnd(svgToPpm(e)); }, [dragging, svgToPpm]);
  const handleMouseUp = useCallback(() => {
    if (!dragging || dragStart === null || dragEnd === null) { setDragging(false); return; }
    const lo = Math.min(dragStart, dragEnd), hi = Math.max(dragStart, dragEnd);
    if (hi - lo > 0.2) { setZoomStack(s => [...s, [...viewRange]]); setViewRange([Math.max(lo, fullRange[0]), Math.min(hi, fullRange[1])]); }
    setDragging(false); setDragStart(null); setDragEnd(null);
  }, [dragging, dragStart, dragEnd, viewRange, fullRange]);
  const resetZoom = () => { setViewRange([...fullRange]); setZoomStack([]); };
  const zoomBack = () => { if (zoomStack.length) { setViewRange(zoomStack[zoomStack.length - 1]); setZoomStack(s => s.slice(0, -1)); } };
  const isZoomed = viewRange[0] !== fullRange[0] || viewRange[1] !== fullRange[1];
  const dragLo = dragStart !== null && dragEnd !== null ? Math.min(dragStart, dragEnd) : 0;
  const dragHi = dragStart !== null && dragEnd !== null ? Math.max(dragStart, dragEnd) : 0;

  const predSolvLabel = predSolvent === "none" ? null : SOLVENT_DB[predSolvent]?.label || predSolvent;
  const actSolvLabel = solvData.label;

  const st = { bg: "#0f172a", card: "#1e293b", bdr: "#334155", tx: "#e2e8f0", mt: "#94a3b8", btn: "#3b82f6" };
  const Btn = ({ active, children, onClick, style: s2 }) => (
    <button onClick={onClick} style={{ padding: "8px 16px", borderRadius: 8, border: "none", cursor: "pointer", background: active ? st.btn : st.card, color: active ? "#fff" : st.mt, fontWeight: 600, fontSize: 13, ...(s2 || {}) }}>{children}</button>
  );
  const Sel = ({ value, onChange, children, style: s2 }) => (
    <select value={value} onChange={onChange} style={{ background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 8px", fontSize: 13, ...(s2 || {}) }}>{children}</select>
  );

  return (
    <div style={{ fontFamily: "system-ui,sans-serif", background: st.bg, color: st.tx, minHeight: "100vh", padding: 16 }}>
      <h1 style={{ fontSize: 22, fontWeight: 700, margin: 0 }}>¹³C NMR Spectrum Analyzer</h1>
      <p style={{ color: st.mt, fontSize: 13, marginBottom: 16 }}>Paste MNOVA peak tables & ChemDraw predictions and auto-compare</p>
      <div style={{ display: "flex", gap: 8, marginBottom: 16, flexWrap: "wrap", alignItems: "center" }}>
        {TABS.map(t => <Btn key={t} active={tab === t} onClick={() => setTab(t)}>{t === "input" ? "📋 Input" : t === "spectrum" ? "📊 Spectrum" : "🔬 Analysis"}</Btn>)}
        {allPeaks.length > 0 && <span style={{ color: st.mt, fontSize: 12, marginLeft: 8 }}>{allPeaks.length} peaks · {compData.filter(c => c.predicted.length > 0).length} compounds</span>}
      </div>

      {/* ===== INPUT TAB ===== */}
      {tab === "input" && (
        <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
          <div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
            <h3 style={{ fontSize: 15, fontWeight: 700, marginBottom: 8 }}>MNOVA Peak Table</h3>
            <p style={{ fontSize: 12, color: st.mt, marginBottom: 8 }}>Paste tab-separated peak table (index, ppm, intensity, width, area, type, flags, impurity/compound, annotation)</p>
            <textarea value={mnovaText} onChange={e => setMnovaText(e.target.value)}
              placeholder={"1\t188.81\t1099.8\t1.44\t9171.18\tCompound\tNone\t\t\n2\t77.16\t16109.1\t2.27\t205295.30\tSolvent\tNone\tCDCl3\t"}
              style={{ width: "100%", height: 150, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 8, padding: 10, fontSize: 12, fontFamily: "monospace", resize: "vertical", boxSizing: "border-box" }} />
            {allPeaks.length > 0 && <p style={{ fontSize: 12, color: "#4ade80", marginTop: 6 }}>✓ {allPeaks.length} peaks ({solvPeaks.length} solvent, {compPeaks.length} compound)</p>}
          </div>

          <div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 12 }}>
              <h3 style={{ fontSize: 15, fontWeight: 700, margin: 0 }}>ChemDraw Predictions</h3>
              <button onClick={addComp} style={{ background: st.btn, color: "#fff", border: "none", borderRadius: 8, padding: "6px 14px", cursor: "pointer", fontSize: 13, fontWeight: 600 }}>+ Add Compound</button>
            </div>
            {compounds.map((c, i) => (
              <div key={i} style={{ background: "#0f172a", borderRadius: 8, padding: 12, marginBottom: 10, borderLeft: `4px solid ${c.color}` }}>
                <div style={{ display: "flex", gap: 8, alignItems: "center", marginBottom: 8, flexWrap: "wrap" }}>
                  <input value={c.name} onChange={e => updComp(i, "name", e.target.value)}
                    style={{ background: st.card, color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 10px", fontSize: 13, fontWeight: 600, flex: 1, minWidth: 100 }} />
                  <select value={c.color} onChange={e => updComp(i, "color", e.target.value)}
                    style={{ background: st.card, color: c.color, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 8px", fontSize: 13, fontWeight: 700 }}>
                    {COLORS.map(cl => <option key={cl} value={cl} style={{ color: cl }}>{cl}</option>)}
                  </select>
                  <label style={{ fontSize: 12, display: "flex", alignItems: "center", gap: 4 }}>
                    <input type="checkbox" checked={c.isProduct} onChange={e => updComp(i, "isProduct", e.target.checked)} /> Major product
                  </label>
                  {compounds.length > 1 && <button onClick={() => rmComp(i)} style={{ background: "#7f1d1d", color: "#fca5a5", border: "none", borderRadius: 6, padding: "4px 10px", cursor: "pointer", fontSize: 12 }}>✕</button>}
                </div>
                <textarea value={c.text} onChange={e => updComp(i, "text", e.target.value)} placeholder="Paste ChemDraw C-13 NMR Prediction here..."
                  style={{ width: "100%", height: 100, background: st.card, color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: 8, fontSize: 11, fontFamily: "monospace", resize: "vertical", boxSizing: "border-box" }} />
                {compData[i]?.predicted.length > 0 && <p style={{ fontSize: 12, color: "#4ade80", marginTop: 4 }}>✓ {compData[i].predicted.length} predicted carbons</p>}
              </div>
            ))}
          </div>

          <div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
            <h3 style={{ fontSize: 15, fontWeight: 700, marginBottom: 10 }}>Solvent & Impurities</h3>
            <div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center", marginBottom: 12, fontSize: 13 }}>
              <label>NMR solvent: <Sel value={solvent} onChange={e => { setSolvent(e.target.value); setEnabledImp({}); }}>
                {Object.entries(SOLVENT_DB).map(([k, v]) => <option key={k} value={k}>{v.label}</option>)}
              </Sel></label>
              <label>Prediction solvent: <Sel value={predSolvent} onChange={e => setPredSolvent(e.target.value)}>
                <option value="none">Unknown / N/A</option>
                {Object.entries(SOLVENT_DB).map(([k, v]) => <option key={k} value={k}>{v.label}</option>)}
              </Sel></label>
            </div>
            {solventPrediction && (
              <div style={{ background: "#1a1a2e", border: "1px solid #4f46e5", borderRadius: 8, padding: 12, marginBottom: 12 }}>
                <div style={{ fontSize: 13, fontWeight: 700, marginBottom: 8, color: "#818cf8" }}>
                  Solvent Shift Predictions: {SOLVENT_LABELS[predSolvent]} {"→"} {SOLVENT_LABELS[solvent]}
                </div>
                <div style={{ display: "flex", gap: 4, marginBottom: 8, fontSize: 12, flexWrap: "wrap", alignItems: "center" }}>
                  <span style={{ padding: "2px 8px", borderRadius: 12, background: "#312e81", color: "#c7d2fe" }}>
                    {solventPrediction.predCls?.desc} (Class {solventPrediction.predCls?.cls})
                  </span>
                  <span style={{ color: st.mt }}>{"→"}</span>
                  <span style={{ padding: "2px 8px", borderRadius: 12, background: "#312e81", color: "#c7d2fe" }}>
                    {solventPrediction.actCls?.desc} (Class {solventPrediction.actCls?.cls})
                  </span>
                </div>
                <div style={{ fontSize: 12, color: st.tx, marginBottom: 6 }}>
                  <span style={{ fontWeight: 600 }}>Expected 13C shift (actual - predicted):</span>
                </div>
                <div style={{ display: "grid", gridTemplateColumns: "repeat(3, 1fr)", gap: 8, marginBottom: 8 }}>
                  {[
                    { label: "Carbonyl C", data: solventPrediction.carbonyl, color: "#ef4444" },
                    { label: "Aromatic/sp2 C", data: solventPrediction.aromatic, color: "#f59e0b" },
                    { label: "Aliphatic C", data: solventPrediction.aliphatic, color: "#22c55e" },
                  ].map(({ label, data, color }) => (
                    <div key={label} style={{ background: "#0f172a", borderRadius: 6, padding: 8, borderLeft: `3px solid ${color}` }}>
                      <div style={{ fontSize: 11, color: st.mt, marginBottom: 4 }}>{label}</div>
                      {data ? (
                        <>
                          <div style={{ fontSize: 14, fontWeight: 700, color }}>
                            {data.avg >= 0 ? "+" : ""}{data.avg.toFixed(1)} ppm
                          </div>
                          <div style={{ fontSize: 10, color: st.mt }}>
                            range: {data.min >= 0 ? "+" : ""}{data.min.toFixed(1)} to {data.max >= 0 ? "+" : ""}{data.max.toFixed(1)}
                          </div>
                        </>
                      ) : (
                        <div style={{ fontSize: 11, color: st.mt }}>No data</div>
                      )}
                    </div>
                  ))}
                </div>
                {solventPrediction.guidance.map((g, i) => (
                  <p key={i} style={{ fontSize: 11, color: "#fef3c7", margin: "4px 0", lineHeight: 1.4 }}>* {g}</p>
                ))}
                <p style={{ fontSize: 11, color: st.mt, marginTop: 6 }}>
                  Computed from {Object.keys(impuritiesData.compounds).length} reference compounds. Consider adjusting tolerance to account for solvent effects.
                </p>
              </div>
            )}
            <p style={{ fontSize: 12, color: st.mt, marginBottom: 8 }}>Known impurities in {solvData.label} ({solvData.impurities.length} compounds):</p>
            <div style={{ display: "flex", flexWrap: "wrap", gap: 6, maxHeight: 240, overflowY: "auto", padding: 4 }}>
              {solvData.impurities.map((ki, i) => {
                const en = !!enabledImp[ki.name];
                const keys = [ki.name.toLowerCase(), ...(COMPOUND_ALIASES[ki.name] || [])];
                const annotated = allPeaks.some(p => {
                  const fields = [p.annotation, p.type, p.compound, p.flags].filter(Boolean).map(s => s.toLowerCase());
                  return fields.some(f => keys.some(k => f.includes(k)));
                });
                const displayName = ki.name.charAt(0).toUpperCase() + ki.name.slice(1);
                return (
                  <button key={i} onClick={() => setEnabledImp(prev => ({ ...prev, [ki.name]: !en }))}
                    style={{ padding: "5px 12px", borderRadius: 8, border: en ? "2px solid #d946ef" : `1px solid ${st.bdr}`, background: en ? "#2e1065" : "#0f172a", color: en ? "#e9d5ff" : st.mt, cursor: "pointer", fontSize: 12, fontWeight: en ? 600 : 400 }}>
                    {displayName} {ki.ppm.length > 0 && <span style={{ opacity: 0.6 }}>({ki.ppm.slice(0, 3).map(p => p.toFixed(1)).join(", ")}{ki.ppm.length > 3 ? ", ..." : ""})</span>}
                    {annotated && <span style={{ color: "#4ade80", marginLeft: 4 }}>*</span>}
                  </button>
                );
              })}
            </div>
            <p style={{ fontSize: 11, color: st.mt, marginTop: 6 }}>
              * = found in MNOVA annotations (annotated ppm used directly).
              Non-annotated impurities are matched within ±{IMP_TOL} ppm of the table reference shift.
              Shifts sourced from experimental ¹³C data table.
            </p>
          </div>

          <div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
            <h3 style={{ fontSize: 15, fontWeight: 700, marginBottom: 10 }}>Settings</h3>
            <div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center", fontSize: 13 }}>
              <label>Tolerance: <input type="number" value={tol} onChange={e => setTol(parseFloat(e.target.value) || 6)} min={0.5} max={20} step={0.5}
                style={{ width: 55, marginLeft: 4, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} /> ppm</label>
              <label>Linewidth: <input type="number" value={linewidth} onChange={e => setLinewidth(parseFloat(e.target.value) || 0.05)} min={0.001} max={2} step={0.001}
                style={{ width: 65, marginLeft: 4, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} /> ppm</label>
              <label>Range: <input type="number" value={fullRange[0]} onChange={e => { const v = parseFloat(e.target.value) || 0; setFullRange([v, fullRange[1]]); setViewRange([v, fullRange[1]]); setZoomStack([]); }}
                style={{ width: 50, marginLeft: 4, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} />–
                <input type="number" value={fullRange[1]} onChange={e => { const v = parseFloat(e.target.value) || 220; setFullRange([fullRange[0], v]); setViewRange([fullRange[0], v]); setZoomStack([]); }}
                  style={{ width: 50, marginLeft: 2, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} /></label>
              <label><input type="checkbox" checked={showSolvent} onChange={e => setShowSolvent(e.target.checked)} /> Solvent</label>
              <label><input type="checkbox" checked={scaleToComp} onChange={e => setScaleToComp(e.target.checked)} /> Scale to compound</label>
            </div>
          </div>

          {allPeaks.length > 0 && <div style={{ textAlign: "center" }}>
            <button onClick={() => setTab("spectrum")} style={{ background: st.btn, color: "#fff", border: "none", borderRadius: 10, padding: "12px 32px", cursor: "pointer", fontSize: 15, fontWeight: 700 }}>View Spectrum →</button>
          </div>}
        </div>
      )}

      {/* ===== SPECTRUM TAB ===== */}
      {tab === "spectrum" && (
        <div>
          {allPeaks.length === 0 ? (
            <div style={{ background: st.card, borderRadius: 12, padding: 40, textAlign: "center" }}><p style={{ color: st.mt }}>No peaks. Go to Input tab.</p></div>
          ) : (
            <div style={{ background: st.card, borderRadius: 12, padding: 16, overflowX: "auto" }}>
              <div style={{ display: "flex", gap: 8, marginBottom: 8, alignItems: "center", flexWrap: "wrap" }}>
                <Btn active={zoomMode} onClick={() => setZoomMode(!zoomMode)} style={{ background: zoomMode ? "#7c3aed" : st.card }}>
                  {zoomMode ? "🔍 Drag to Zoom" : "🔍 Zoom"}
                </Btn>
                {isZoomed && <Btn onClick={zoomBack}>← Back</Btn>}
                {isZoomed && <Btn onClick={resetZoom}>⊞ Full</Btn>}
                {isZoomed && <span style={{ fontSize: 12, color: st.mt }}>{viewRange[0].toFixed(1)}–{viewRange[1].toFixed(1)} ppm</span>}
              </div>

              {/* Compound visibility toggles */}
              {matchRes.filter(c => c.predicted.length > 0).length > 0 && (
                <div style={{ display: "flex", gap: 6, marginBottom: 10, flexWrap: "wrap", alignItems: "center" }}>
                  <span style={{ fontSize: 12, color: st.mt }}>Show overlays:</span>
                  {matchRes.filter(c => c.predicted.length > 0).map((c, i) => {
                    const vis = isCompVisible(c.name);
                    return (
                      <button key={i} onClick={() => toggleCompVis(c.name)}
                        style={{ padding: "4px 12px", borderRadius: 8, border: vis ? `2px solid ${c.color}` : `1px solid ${st.bdr}`,
                          background: vis ? st.card : "#0f172a", color: vis ? c.color : st.mt,
                          cursor: "pointer", fontSize: 12, fontWeight: vis ? 700 : 400, opacity: vis ? 1 : 0.5 }}>
                        {vis ? "●" : "○"} {c.name}
                      </button>
                    );
                  })}
                </div>
              )}

              <svg ref={svgRef} viewBox={`0 0 ${svgW} ${svgH}`}
                style={{ width: "100%", height: "auto", minWidth: 600, cursor: zoomMode ? "crosshair" : "default", userSelect: "none" }}
                onMouseDown={handleMouseDown} onMouseMove={handleMouseMove} onMouseUp={handleMouseUp}
                onMouseLeave={() => { if (dragging) { setDragging(false); setDragStart(null); setDragEnd(null); } }}>

                <rect x={pad.l} y={pad.t} width={plotW} height={plotH} fill="#0f172a" rx={4} />

                {dragging && dragStart !== null && dragEnd !== null && (
                  <rect x={Math.min(xSc(dragLo), xSc(dragHi))} y={pad.t} width={Math.abs(xSc(dragHi) - xSc(dragLo))} height={plotH} fill="#7c3aed" opacity={0.15} stroke="#7c3aed" strokeWidth={1} />
                )}

                {(() => {
                  const span = viewRange[1] - viewRange[0];
                  const step = span > 150 ? 20 : span > 60 ? 10 : span > 30 ? 5 : span > 10 ? 2 : span > 4 ? 1 : span > 1 ? 0.5 : 0.1;
                  const start = Math.ceil(viewRange[0] / step) * step;
                  const lines = [];
                  for (let v = start; v <= viewRange[1]; v += step) { const x = xSc(v); if (x >= pad.l && x <= pad.l + plotW) lines.push(v); }
                  return lines.map(v => (
                    <g key={v}>
                      <line x1={xSc(v)} x2={xSc(v)} y1={pad.t} y2={pad.t + plotH} stroke={st.bdr} strokeWidth={0.5} strokeDasharray="3,6" />
                      <text x={xSc(v)} y={svgH - 24} textAnchor="middle" fill={st.mt} fontSize={11}>{step < 1 ? v.toFixed(1) : v}</text>
                    </g>
                  ));
                })()}
                <text x={svgW / 2} y={svgH - 6} textAnchor="middle" fill={st.mt} fontSize={10}>Chemical Shift (ppm)</text>

                {showSolvent && nmrSolvPeaks.map((p, i) => (
                  <rect key={`sv${i}`} x={xSc(p.ppm) - 3} y={pad.t} width={6} height={plotH} fill="#475569" opacity={0.1} />
                ))}
                {showSolvent && impSolvPeaks.map((p, i) => (
                  <g key={`isv${i}`}>
                    <rect x={xSc(p.ppm) - 3} y={pad.t} width={6} height={plotH} fill="#f59e0b" opacity={0.08} />
                    <text x={xSc(p.ppm)} y={pad.t + plotH - 6} textAnchor="middle" fill="#f59e0b" fontSize={7} opacity={0.5}>{p.compound}</text>
                  </g>
                ))}
                {(() => {
                  const items = [];
                  resolvedImpurities.forEach((ki, ki_i) => ki.entries.forEach((e, e_i) => {
                    const x = xSc(e.refPpm);
                    items.push({ x, ki_i, e_i, name: ki.name, hasMatch: e.hasMatch, isAnnotated: e.isAnnotated, matchedPeak: e.matchedPeak, refPpm: e.refPpm });
                  }));
                  items.sort((a, b) => a.x - b.x);
                  const placed = [];
                  for (const lb of items) {
                    let yOff = 0, show = false;
                    while (yOff <= 80) {
                      const y = pad.t + 14 + yOff;
                      if (!placed.some(p => {
                        const minDist = (lb.name.length + p.name.length) * 2.2 + 8;
                        return Math.abs(p.x - lb.x) < minDist && Math.abs(p.y - y) < 10;
                      })) { show = true; break; }
                      yOff += 11;
                    }
                    placed.push({ ...lb, y: pad.t + 14 + yOff, show });
                  }
                  return placed.map(lb => {
                    const lineOpacity = lb.isAnnotated ? 0.6 : lb.hasMatch ? 0.45 : 0.35;
                    const lineColor = lb.hasMatch || lb.isAnnotated ? "#d946ef" : "#94a3b8";
                    const labelColor = lb.hasMatch ? "#d946ef" : "#94a3b8";
                    const labelOpacity = lb.hasMatch ? 0.75 : 0.55;
                    return (
                      <g key={`ki${lb.ki_i}${lb.e_i}`}>
                        {/* Reference dashed line — bright if matched/annotated, dim if no match */}
                        <line x1={lb.x} x2={lb.x} y1={pad.t} y2={pad.t + plotH}
                          stroke={lineColor} strokeWidth={1} strokeDasharray="5,5" opacity={lineOpacity} />
                        {/* Connecting line to actual matched peak on the compound curve */}
                        {lb.hasMatch && !lb.isAnnotated && lb.matchedPeak && (
                          <line x1={lb.x} x2={xSc(lb.matchedPeak.ppm)}
                            y1={pad.t} y2={ySc(Math.min(evalAt(lb.matchedPeak.ppm, compPeaks, linewidth), displayMaxY))}
                            stroke="#a855f7" strokeWidth={1.5} opacity={0.45} strokeDasharray="3,3" />
                        )}
                        {lb.show && (
                          <text x={lb.x} y={lb.y} textAnchor="middle" fill={labelColor} fontSize={8} opacity={labelOpacity}>
                            {lb.name}{lb.isAnnotated && <tspan fill="#4ade80"> ●</tspan>}
                          </text>
                        )}
                      </g>
                    );
                  });
                })()}

                {showSolvent && nmrSolvCurve.length > 0 && <path d={makePath(nmrSolvCurve)} fill="none" stroke="#64748b" strokeWidth={1} opacity={0.5} />}
                {showSolvent && impSolvCurve.length > 0 && <path d={makePath(impSolvCurve)} fill="none" stroke="#f59e0b" strokeWidth={1.2} opacity={0.6} />}
                {compCurve.length > 0 && <path d={makePath(compCurve)} fill="none" stroke="#3b82f6" strokeWidth={1.5} />}

                {/* Marker zone */}
                <rect x={pad.l} y={4} width={plotW} height={mzH - 8} fill="#1e293b" rx={6} opacity={0.5} />

                {matchRes.filter(c => c.predicted.length > 0).map((comp, ci) => {
                  const vis = isCompVisible(comp.name);
                  const yB = 18 + ci * 24;
                  return (
                    <g key={`cg${ci}`} opacity={vis ? 1 : 0.15}>
                      <text x={pad.l - 4} y={yB + 4} textAnchor="end" fill={comp.color} fontSize={10} fontWeight="600">{comp.name}</text>
                      <line x1={pad.l} x2={pad.l + plotW} y1={yB} y2={yB} stroke={comp.color} strokeWidth={0.4} opacity={0.2} />
                      {comp.matches.map((m, mi) => {
                        const x = xSc(m.shift);
                        if (x < pad.l - 10 || x > pad.l + plotW + 10) return null;
                        return (
                          <g key={`m${ci}${mi}`}>
                            <circle cx={x} cy={yB} r={4} fill={m.matched ? comp.color : "transparent"} stroke={comp.color} strokeWidth={1.5} opacity={0.9} />
                            <text x={x} y={yB - 7} textAnchor="middle" fill={comp.color} fontSize={7.5} opacity={0.8}>{m.shift.toFixed(0)}</text>
                            {vis && <line x1={x} x2={x} y1={yB + 5} y2={pad.t} stroke={comp.color} strokeWidth={0.6} opacity={0.2} strokeDasharray="2,3" />}
                            {vis && m.matched && m.match && (
                              <line x1={x} x2={xSc(m.match.ppm)} y1={pad.t} y2={ySc(Math.min(evalAt(m.match.ppm, compPeaks, linewidth), displayMaxY))}
                                stroke={comp.color} strokeWidth={1.8} opacity={0.55} strokeDasharray="4,3" />
                            )}
                          </g>
                        );
                      })}
                    </g>
                  );
                })}

                {(() => {
                  const labels = compPeaks
                    .filter(p => p.ppm >= viewRange[0] && p.ppm <= viewRange[1])
                    .filter(p => p.intensity > compMaxY * 0.08)
                    .map(p => ({ ppm: p.ppm, x: xSc(p.ppm), baseY: ySc(Math.min(evalAt(p.ppm, compPeaks, linewidth), displayMaxY)) - 5, text: p.ppm.toFixed(1) }))
                    .sort((a, b) => a.x - b.x);
                  const placed = [];
                  for (const lb of labels) {
                    let yOff = 0, ok = false;
                    while (yOff <= 35) {
                      const y = lb.baseY - yOff;
                      if (y < pad.t - 4) break;
                      if (!placed.some(p => Math.abs(p.x - lb.x) < 28 && Math.abs(p.y - y) < 10)) { ok = true; break; }
                      yOff += 11;
                    }
                    if (ok) placed.push({ ...lb, y: lb.baseY - yOff });
                  }
                  return placed.map((lb, i) => (
                    <text key={`pl${i}`} x={lb.x} y={lb.y} textAnchor="middle" fill="#93c5fd" fontSize={8} opacity={0.65}>{lb.text}</text>
                  ));
                })()}
              </svg>

              <div style={{ display: "flex", gap: 14, flexWrap: "wrap", marginTop: 10, fontSize: 12 }}>
                <span><span style={{ color: "#3b82f6" }}>━━</span> Compound</span>
                {showSolvent && <span><span style={{ color: "#64748b" }}>━━</span> NMR Solvent</span>}
                {showSolvent && impSolvPeaks.length > 0 && <span><span style={{ color: "#f59e0b" }}>━━</span> Solvent (impurity)</span>}
                <span><span style={{ color: "#d946ef" }}>┄┄</span> Impurities</span>
              </div>
            </div>
          )}
        </div>
      )}

      {/* ===== ANALYSIS TAB ===== */}
      {tab === "analysis" && (
        <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
          {matchRes.filter(c => c.predicted.length > 0).map((comp, ci) => {
            const uniq = comp.matches;
            const mc = uniq.filter(m => m.matched).length;
            const pct = uniq.length > 0 ? Math.round(mc / uniq.length * 100) : 0;
            const avgScore = mc > 0 ? uniq.filter(m => m.matched).reduce((s, m) => s + scoreMatch(m, maxArea, comp.isProduct), 0) / mc : 0;
            return (
              <div key={ci} style={{ background: st.card, borderRadius: 12, padding: 16, borderLeft: `4px solid ${comp.color}` }}>
                <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 4, flexWrap: "wrap", gap: 8 }}>
                  <h3 style={{ fontSize: 16, fontWeight: 700, color: comp.color, margin: 0 }}>
                    {comp.name} {comp.isProduct && <span style={{ fontSize: 11, color: st.mt }}>(major)</span>}
                  </h3>
                  <div style={{ display: "flex", gap: 8 }}>
                    <span style={{ padding: "4px 12px", borderRadius: 20, fontSize: 12, fontWeight: 700, background: pct >= 70 ? "#166534" : pct >= 40 ? "#854d0e" : "#7f1d1d", color: "#fff" }}>
                      {mc}/{uniq.length} ({pct}%)
                    </span>
                    {mc > 0 && <span style={{ padding: "4px 12px", borderRadius: 20, fontSize: 12, fontWeight: 600, background: avgScore > 0.7 ? "#166534" : avgScore > 0.4 ? "#854d0e" : "#7f1d1d", color: "#fff" }}>
                      Score: {(avgScore * 100).toFixed(0)}%
                    </span>}
                  </div>
                </div>
                {comp.isProduct && <p style={{ fontSize: 11, color: st.mt, marginBottom: 6 }}>Score weights area since major product peaks should dominate the spectrum</p>}
                <table style={{ width: "100%", fontSize: 12, borderCollapse: "collapse" }}>
                  <thead><tr style={{ color: st.mt, textAlign: "left", borderBottom: `1px solid ${st.bdr}` }}>
                    <th style={{ padding: "6px 8px" }}>Type</th><th style={{ padding: "6px 8px" }}>Predicted</th><th style={{ padding: "6px 8px" }}>Matched</th><th style={{ padding: "6px 8px" }}>Δ</th><th style={{ padding: "6px 8px" }}>Area</th><th style={{ padding: "6px 8px" }}>Score</th><th style={{ padding: "6px 8px" }}>Status</th>
                  </tr></thead>
                  <tbody>{uniq.map((m, j) => {
                    const sc = m.matched ? scoreMatch(m, maxArea, comp.isProduct) : 0;
                    return (
                      <tr key={j} style={{ borderBottom: "1px solid #1e293b" }}>
                        <td style={{ padding: "6px 8px", fontWeight: 600 }}>{m.type}</td>
                        <td style={{ padding: "6px 8px" }}>{m.shift.toFixed(1)}</td>
                        <td style={{ padding: "6px 8px" }}>{m.match ? m.match.ppm.toFixed(2) : "—"}</td>
                        <td style={{ padding: "6px 8px" }}>{m.delta.toFixed(1)}</td>
                        <td style={{ padding: "6px 8px" }}>{m.match && m.matched ? m.match.area.toFixed(0) : "—"}</td>
                        <td style={{ padding: "6px 8px" }}>
                          {m.matched && <div style={{ display: "flex", alignItems: "center", gap: 4 }}>
                            <div style={{ width: 50, height: 8, background: "#0f172a", borderRadius: 4, overflow: "hidden" }}>
                              <div style={{ width: `${sc * 100}%`, height: "100%", background: sc > 0.7 ? "#22c55e" : sc > 0.4 ? "#f59e0b" : "#ef4444", borderRadius: 4 }} />
                            </div>
                            <span style={{ fontSize: 10, color: st.mt }}>{(sc * 100).toFixed(0)}</span>
                          </div>}
                        </td>
                        <td style={{ padding: "6px 8px" }}>
                          <span style={{ padding: "2px 10px", borderRadius: 12, fontSize: 11, fontWeight: 600, background: m.matched ? (m.delta < 3 ? "#166534" : "#854d0e") : "#7f1d1d", color: "#fff" }}>
                            {m.matched ? (m.delta < 3 ? "✓ Good" : "~ Approx") : "✕ No match"}
                          </span>
                        </td>
                      </tr>
                    );
                  })}</tbody>
                </table>
                <p style={{ fontSize: 11, color: st.mt, marginTop: 6 }}>
                  Peaks are matched 1:1, so overlapping carbons at the same ppm will show as unmatched. Manually review these unmatched carbons.
                </p>
              </div>
            );
          })}

          <div style={{ background: st.card, borderRadius: 12, padding: 16, borderLeft: "4px solid #64748b" }}>
            <h3 style={{ fontSize: 16, fontWeight: 700, color: "#94a3b8", marginBottom: 10 }}>Unassigned Peaks ({unassigned.length})</h3>
            {unassigned.length === 0 ? <p style={{ color: "#4ade80", fontSize: 13 }}>All peaks accounted for!</p> : (
              <div style={{ display: "flex", flexWrap: "wrap", gap: 8 }}>{unassigned.map((p, i) => (
                <div key={i} style={{ background: "#0f172a", borderRadius: 8, padding: "6px 12px", fontSize: 12 }}>
                  <span style={{ fontWeight: 700, color: "#fbbf24" }}>{p.ppm.toFixed(2)}</span>
                  <span style={{ color: st.mt, marginLeft: 4 }}>I:{p.intensity.toFixed(0)}</span>
                  <span style={{ color: st.mt, marginLeft: 4 }}>A:{p.area.toFixed(0)}</span>
                </div>
              ))}</div>
            )}
            {compPeaks.length > 0 && (() => {
              const W = 700, H = 120, P = { l: 10, r: 10, t: 8, b: 24 };
              const pw = W - P.l - P.r, ph = H - P.t - P.b;
              const allP = [...compPeaks, ...solvPeaks];
              const r = [Math.min(...allP.map(p => p.ppm)) - 5, Math.max(...allP.map(p => p.ppm)) + 5];
              const n = 1500, gm = 0.1;
              const asgn = compPeaks.filter(p => allAssigned.has(p.ppm));
              const unasgn = compPeaks.filter(p => !allAssigned.has(p.ppm));
              const aC = asgn.length ? genSpectrum(asgn, r[0], r[1], n, gm) : [];
              const uC = unasgn.length ? genSpectrum(unasgn, r[0], r[1], n, gm) : [];
              const sC = solvPeaks.length ? genSpectrum(solvPeaks, r[0], r[1], n, gm) : [];
              const maxY = Math.max(...aC.map(p => p[1]), ...uC.map(p => p[1]), ...sC.map(p => p[1]), 1) * 1.1;
              const xS = v => P.l + (r[1] - v) / (r[1] - r[0]) * pw;
              const yS = v => P.t + ph - (v / maxY) * ph;
              const mp = c => { let d='',w=false; for(const p of c){if(p[1]>maxY*.002){d+=(w?'L':'M')+xS(p[0]).toFixed(1)+','+yS(p[1]).toFixed(1)+' ';w=true}else w=false} return d };
              const span = r[1] - r[0], step = span > 150 ? 20 : span > 60 ? 10 : span > 30 ? 5 : 2;
              const ticks = []; for (let v = Math.ceil(r[0] / step) * step; v <= r[1]; v += step) ticks.push(v);
              return (
                <svg viewBox={`0 0 ${W} ${H}`} style={{ width: "100%", height: "auto", marginTop: 12 }}>
                  <rect x={P.l} y={P.t} width={pw} height={ph} fill="#0f172a" rx={4} />
                  {ticks.map(v => <g key={v}><line x1={xS(v)} x2={xS(v)} y1={P.t} y2={P.t + ph} stroke="#334155" strokeWidth={0.5} strokeDasharray="2,4" /><text x={xS(v)} y={H - 6} textAnchor="middle" fill="#64748b" fontSize={8}>{v}</text></g>)}
                  {sC.length > 0 && <path d={mp(sC)} fill="none" stroke="#64748b" strokeWidth={1} opacity={0.5} />}
                  {aC.length > 0 && <path d={mp(aC)} fill="none" stroke="#3b82f6" strokeWidth={1.5} />}
                  {uC.length > 0 && <path d={mp(uC)} fill="none" stroke="#ef4444" strokeWidth={1.5} />}
                  <g transform={`translate(${P.l + 8},${P.t + 14})`}>
                    <line x1={0} x2={16} y1={0} y2={0} stroke="#3b82f6" strokeWidth={2} /><text x={20} y={4} fill="#94a3b8" fontSize={9}>Assigned</text>
                    <line x1={76} x2={92} y1={0} y2={0} stroke="#ef4444" strokeWidth={2} /><text x={96} y={4} fill="#94a3b8" fontSize={9}>Unassigned</text>
                    {solvPeaks.length > 0 && <><line x1={168} x2={184} y1={0} y2={0} stroke="#64748b" strokeWidth={2} /><text x={188} y={4} fill="#94a3b8" fontSize={9}>Solvent</text></>}
                  </g>
                </svg>
              );
            })()}
          </div>

          <div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
            <h3 style={{ fontSize: 16, fontWeight: 700, marginBottom: 10 }}>Peak Groups</h3>
            <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(180px, 1fr))", gap: 8 }}>
              {groups.map((g, i) => (
                <div key={i} style={{ background: "#0f172a", borderRadius: 8, padding: 10, borderLeft: `3px solid ${GROUP_COLORS[i % GROUP_COLORS.length]}` }}>
                  <div style={{ fontWeight: 700, fontSize: 16, color: GROUP_COLORS[i % GROUP_COLORS.length] }}>~{g.center.toFixed(1)}</div>
                  <div style={{ fontSize: 11, color: st.mt }}>{g.peaks.length} peak{g.peaks.length > 1 ? "s" : ""} · {g.min.toFixed(1)}–{g.max.toFixed(1)}</div>
                  <div style={{ fontSize: 11, color: st.mt }}>Area: {g.totalArea.toFixed(0)}</div>
                </div>
              ))}
            </div>
            {groups.length > 0 && (() => {
              const W = 700, H = 120, P = { l: 10, r: 10, t: 8, b: 24 };
              const pw = W - P.l - P.r, ph = H - P.t - P.b;
              const r = [Math.min(...compPeaks.map(p => p.ppm)) - 5, Math.max(...compPeaks.map(p => p.ppm)) + 5];
              const n = 1500, gm = 0.1;
              const curves = groups.map(gr => gr.peaks.length ? genSpectrum(gr.peaks, r[0], r[1], n, gm) : []);
              const maxY = Math.max(...curves.flatMap(c => c.map(p => p[1])), 1) * 1.1;
              const xS = v => P.l + (r[1] - v) / (r[1] - r[0]) * pw;
              const yS = v => P.t + ph - (v / maxY) * ph;
              const mp = c => { let d='',w=false; for(const p of c){if(p[1]>maxY*.002){d+=(w?'L':'M')+xS(p[0]).toFixed(1)+','+yS(p[1]).toFixed(1)+' ';w=true}else w=false} return d };
              const span = r[1] - r[0], step = span > 150 ? 20 : span > 60 ? 10 : span > 30 ? 5 : 2;
              const ticks = []; for (let v = Math.ceil(r[0] / step) * step; v <= r[1]; v += step) ticks.push(v);
              return (
                <svg viewBox={`0 0 ${W} ${H}`} style={{ width: "100%", height: "auto", marginTop: 12 }}>
                  <rect x={P.l} y={P.t} width={pw} height={ph} fill="#0f172a" rx={4} />
                  {ticks.map(v => <g key={v}><line x1={xS(v)} x2={xS(v)} y1={P.t} y2={P.t + ph} stroke="#334155" strokeWidth={0.5} strokeDasharray="2,4" /><text x={xS(v)} y={H - 6} textAnchor="middle" fill="#64748b" fontSize={8}>{v}</text></g>)}
                  {curves.map((c, i) => c.length > 0 && <path key={i} d={mp(c)} fill="none" stroke={GROUP_COLORS[i % GROUP_COLORS.length]} strokeWidth={1.5} />)}
                </svg>
              );
            })()}
          </div>
        </div>
      )}
    </div>
  );
}