import { useState, useMemo, useCallback, useRef } from "react";
import protonImpurities from "./1HimpuritiesTable.json";

const COLORS = ["#22c55e", "#f59e0b", "#ef4444", "#a78bfa", "#ec4899", "#06b6d4", "#f97316", "#84cc16"];
const GROUP_COLORS = ["#22c55e", "#f59e0b", "#a78bfa", "#ec4899", "#06b6d4", "#f97316", "#84cc16", "#14b8a6", "#e879f9", "#fb923c", "#64748b", "#facc15"];

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
	"dimethyl sulfoxide": ["dmso"], "water": ["h2o", "hdo"],
};

const NMR_SOLV_NAMES = {
	"CDCl3": ["cdcl3", "chloroform-d"], "CD2Cl2": ["cd2cl2", "dcm"],
	"(CD3)2SO": ["(cd3)2so", "dmso", "dmso-d6"], "(CD3)2CO": ["(cd3)2co", "acetone", "acetone-d6"],
	"CD3CN": ["cd3cn", "mecn", "acetonitrile"], "CD3OD": ["cd3od", "meod", "methanol"],
	"D2O": ["d2o", "water"], "C6D6": ["c6d6", "benzene"], "C6D5Cl": ["c6d5cl", "chlorobenzene"],
	"THF-d8": ["thf-d8", "thf"], "toluene-d8": ["toluene-d8", "toluene"], "TFE-d3": ["tfe-d3"],
};

// ── Flag parsing ─────────────────────────────────────────────────────────────
function parseFlags(flagStr) {
	if (!flagStr || flagStr === "None") return new Set();
	return new Set(flagStr.split(/\s*\+\s*/).map(f => f.trim().toLowerCase()).filter(Boolean));
}

const FLAG_C13_SAT = "c13_sat";
const FLAG_WEAK = "weak";
const FLAG_HIDDEN = "hidden";
const FLAG_2ARY_SOLVENT = "2ary_solvent";

// ── Build 1H impurity DB ─────────────────────────────────────────────────────
const SOLVENT_DB = {};
for (const solv of protonImpurities.solvents) {
	const impurities = [];
	for (const [compName, groups] of Object.entries(protonImpurities.compounds)) {
		const ppmValues = [];
		for (const shifts of Object.values(groups)) {
			if (shifts[solv] != null) {
				const vals = Array.isArray(shifts[solv]) ? shifts[solv] : [shifts[solv]];
				ppmValues.push(...vals.filter(v => typeof v === 'number'));
			}
		}
		if (ppmValues.length > 0) impurities.push({ name: compName, ppm: ppmValues });
	}
	SOLVENT_DB[solv] = {
		label: SOLVENT_LABELS[solv] || solv,
		solventPeaks: protonImpurities.solvent_signals[solv] || [],
		impurities,
	};
}

// ── Parsers ──────────────────────────────────────────────────────────────────

function parseMNOVA_1H(text) {
	const lines = text.trim().split("\n").filter(l => l.trim());
	const peaks = [];
	for (const line of lines) {
		const parts = line.split("\t");
		if (parts.length < 5 || !/^\d+$/.test(parts[0].trim())) continue;
		const ppm = parseFloat(parts[1]);
		const int_ = parseFloat(parts[2]);
		const w = parseFloat(parts[3]);
		const area = parseFloat(parts[4]);
		const type = (parts[5] || "").trim();
		const flagStr = (parts[6] || "").trim();
		const compound = (parts[7] || "").trim();
		const annotation = (parts[8] || "").trim();
		if (isNaN(ppm)) continue;
		const flags = parseFlags(flagStr);
		peaks.push({
			ppm, intensity: int_, width: isNaN(w) ? 1.5 : w, area: isNaN(area) ? 0 : area,
			type, flagStr, flags, compound, annotation,
			isC13Sat: flags.has(FLAG_C13_SAT),
			isWeak: flags.has(FLAG_WEAK),
			isHidden: flags.has(FLAG_HIDDEN),
			is2arySolvent: flags.has(FLAG_2ARY_SOLVENT),
		});
	}
	return peaks;
}

function parseChemDraw1H(text) {
	const lines = text.split("\n");
	const nodes = [];
	let cur = null;
	for (const line of lines) {
		const m = line.match(/^(CH3|CH2|CH|NH2|NH|OH)\s+([-\d.;]+)\s+([-\d.]+)\s+(.*)/);
		if (m) {
			if (cur) nodes.push(cur);
			const shifts = m[2].split(";").map(s => parseFloat(s.trim())).filter(s => !isNaN(s));
			const hCount = m[1] === "CH3" ? 3 : (m[1] === "CH2" || m[1] === "NH2") ? 2 : 1;
			cur = { type: m[1], shifts, hCount, desc: m[4].trim(), baseShift: parseFloat(m[3]), increments: [m[4].trim()] };
		} else if (cur && line.match(/^\s+[-\d?.]+\s+/)) {
			cur.increments.push(line.trim());
		}
	}
	if (cur) nodes.push(cur);
	return nodes;
}

// ── Spectrum ─────────────────────────────────────────────────────────────────

function lorentzian(x, x0, g, h) { return h * g * g / ((x - x0) * (x - x0) + g * g); }

function genSpectrum(peaks, min, max, n, gamma) {
	const pts = [], step = (max - min) / n;
	for (let i = 0; i <= n; i++) { const x = min + i * step; let y = 0; for (const p of peaks) y += lorentzian(x, p.ppm, gamma, p.intensity); pts.push([x, y]); }
	return pts;
}

function evalAt(x, peaks, gamma) { let y = 0; for (const p of peaks) y += lorentzian(x, p.ppm, gamma, p.intensity); return y; }

// ── Matching ─────────────────────────────────────────────────────────────────

function matchProtonCompound(predicted, actual, tol, isProduct) {
	const flatPred = [];
	for (let ni = 0; ni < predicted.length; ni++) {
		for (let si = 0; si < predicted[ni].shifts.length; si++) {
			flatPred.push({ nodeIdx: ni, shiftIdx: si, shift: predicted[ni].shifts[si], node: predicted[ni] });
		}
	}
	const candidates = [];
	for (let fi = 0; fi < flatPred.length; fi++) {
		for (const a of actual) {
			const d = Math.abs(a.ppm - flatPred[fi].shift);
			if (d <= tol) candidates.push({ fi, pred: flatPred[fi], act: a, delta: d, area: a.area });
		}
	}
	if (isProduct) {
		candidates.sort((a, b) => { if (Math.abs(b.area - a.area) > 100) return b.area - a.area; return a.delta - b.delta; });
	} else {
		candidates.sort((a, b) => a.delta - b.delta);
	}
	const usedPred = new Set(), usedAct = new Set();
	const results = [];
	for (const c of candidates) {
		if (usedPred.has(c.fi) || usedAct.has(c.act.ppm)) continue;
		usedPred.add(c.fi); usedAct.add(c.act.ppm);
		results.push({ ...c.pred, match: c.act, delta: c.delta, matched: true, _fi: c.fi });
	}
	for (let fi = 0; fi < flatPred.length; fi++) {
		if (usedPred.has(fi)) continue;
		let best = null, bestD = Infinity;
		for (const a of actual) { const d = Math.abs(a.ppm - flatPred[fi].shift); if (d < bestD) { bestD = d; best = a; } }
		results.push({ ...flatPred[fi], match: best, delta: bestD, matched: false, _fi: fi });
	}
	results.sort((a, b) => a._fi - b._fi);

	const grouped = predicted.map((node, ni) => {
		const nr = results.filter(r => r.nodeIdx === ni);
		const allMatched = nr.every(r => r.matched);
		const anyMatched = nr.some(r => r.matched);
		const avgDelta = nr.length > 0 ? nr.reduce((s, r) => s + r.delta, 0) / nr.length : Infinity;
		const matchedArea = nr.filter(r => r.matched && r.match).reduce((s, r) => s + r.match.area, 0);
		return { node, nodeIdx: ni, results: nr, allMatched, anyMatched, avgDelta, matchedArea };
	});

	return { flat: results, grouped };
}

function scoreProtonMatch(nm, maxArea, isProduct) {
	if (!nm.anyMatched) return 0;
	const mr = nm.results.filter(r => r.matched);
	const avgDist = mr.reduce((s, r) => s + r.delta, 0) / mr.length;
	const distScore = Math.max(0, 1 - avgDist / 1.0);
	const frac = mr.length / nm.results.length;
	const areaRatio = nm.matchedArea / maxArea;
	if (isProduct) return distScore * 0.4 + Math.min(1, areaRatio * 2.5) * 0.3 + frac * 0.3;
	return distScore * 0.6 + Math.min(1, areaRatio * 3) * 0.15 + frac * 0.25;
}

function groupPeaks(peaks, gap = 0.08) {
	if (!peaks.length) return [];
	const s = [...peaks].sort((a, b) => a.ppm - b.ppm);
	const gs = []; let g = [s[0]];
	for (let i = 1; i < s.length; i++) { if (s[i].ppm - g[g.length - 1].ppm <= gap) g.push(s[i]); else { gs.push(g); g = [s[i]]; } }
	gs.push(g);
	return gs.map(gr => ({ center: gr.reduce((s, p) => s + p.ppm * p.intensity, 0) / gr.reduce((s, p) => s + p.intensity, 0), totalArea: gr.reduce((s, p) => s + p.area, 0), peaks: gr, min: gr[0].ppm, max: gr[gr.length - 1].ppm }));
}

const IMP_TOL = 0.15;
const TABS = ["input", "spectrum", "analysis"];

export default function ProtonNMR() {
	const [tab, setTab] = useState("input");
	const [mnovaText, setMnovaText] = useState("");
	const [compounds, setCompounds] = useState([{ name: "Product", text: "", color: COLORS[0], isProduct: true }]);
	const [tol, setTol] = useState(1.0);
	const [fullRange, setFullRange] = useState([-0.5, 14]);
	const [viewRange, setViewRange] = useState([-0.5, 14]);
	const [linewidth, setLinewidth] = useState(0.003);
	const [showSolvent, setShowSolvent] = useState(true);
	const [showC13Sats, setShowC13Sats] = useState(false);
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

	// ── Parse & categorize peaks ───────────────────────────────────────────────
	const allPeaks = useMemo(() => parseMNOVA_1H(mnovaText), [mnovaText]);
	const c13SatPeaks = useMemo(() => allPeaks.filter(p => p.isC13Sat), [allPeaks]);
	const solvPeaks = useMemo(() => allPeaks.filter(p => p.type === "Solvent" && !p.is2arySolvent), [allPeaks]);
	const secondarySolvPeaks = useMemo(() => allPeaks.filter(p => p.is2arySolvent || (p.type === "Solvent" && p.is2arySolvent)), [allPeaks]);
	const impPeaks = useMemo(() => allPeaks.filter(p => p.type === "Impurity"), [allPeaks]);
	const artifactPeaks = useMemo(() => allPeaks.filter(p => p.type === "Artifact"), [allPeaks]);

	// Compound peaks: include hidden (drawn but not labeled), exclude C13 sats
	const compPeaks = useMemo(() => allPeaks.filter(p => p.type === "Compound" && !p.isC13Sat), [allPeaks]);
	// For matching: use all compound peaks (hidden included — they're real peaks)
	const matchablePeaks = compPeaks;

	const nmrSolvNames = useMemo(() => (NMR_SOLV_NAMES[solvent] || [solvent.toLowerCase()]), [solvent]);
	const nmrSolvPeaks = useMemo(() => solvPeaks.filter(p => {
		if (!p.compound) return true;
		const c = p.compound.toLowerCase();
		return nmrSolvNames.some(n => n === c || c.includes(n));
	}), [solvPeaks, nmrSolvNames]);

	const compData = useMemo(() => compounds.map(c => ({ ...c, predicted: parseChemDraw1H(c.text) })), [compounds]);
	const maxArea = useMemo(() => Math.max(...matchablePeaks.map(p => p.area), 1), [matchablePeaks]);

	const solventPrediction = useMemo(() => {
		if (predSolvent === "none" || predSolvent === solvent) return null;
		const predCls = SOLVENT_CLASS[predSolvent];
		const actCls = SOLVENT_CLASS[solvent];
		const exchangeable = [], aromatic = [], alphaHetero = [], aliphatic = [];

		function classify(groupName, shift) {
			const g = groupName.toLowerCase();
			if (g === "oh" || g === "nh" || g === "nh2" || g === "h2") return "exchangeable";
			if (g.startsWith("ar") || g.startsWith("ch(")) return "aromatic"; // ArH, ArCH3, CH(2,6) ring positions
			if (g.startsWith("och") || g.startsWith("nch") || g.includes("co")) return "alphaHetero"; // OCH3, NCH3, CH3CO
			if (shift > 6) return "aromatic";
			if (shift >= 3) return "alphaHetero";
			return "aliphatic";
		}

		for (const groups of Object.values(protonImpurities.compounds)) {
			for (const [groupName, info] of Object.entries(groups)) {
				const pS = info[predSolvent], aS = info[solvent];
				if (pS == null || aS == null) continue;
				const pVal = Array.isArray(pS) ? pS[0] : pS;
				const aVal = Array.isArray(aS) ? aS[0] : aS;
				if (typeof pVal !== "number" || typeof aVal !== "number") continue;
				const d = aVal - pVal;
				const cls = classify(groupName, pVal);
				if (cls === "exchangeable") exchangeable.push(d);
				else if (cls === "aromatic") aromatic.push(d);
				else if (cls === "alphaHetero") alphaHetero.push(d);
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
			const predAro = predCls.desc.includes("aromatic");
			const actAro = actCls.desc.includes("aromatic");
			if (f === t) {
				if (predAro !== actAro) {
					guidance.push("Same class but aromatic ↔ non-aromatic -- ring-current anisotropy shifts nearby protons. Consider increasing tolerance.");
				} else {
					guidance.push("Same solvent sub-type -- most CH shifts are small. OH/NH may still vary between solvents.");
				}
			} else if (f === "I" && t === "III") {
				guidance.push("Class I → III: OH/NH shift strongly downfield or vanish via H/D exchange. Other protons shift slightly downfield. Widen tolerance near polar groups.");
			} else if (f === "III" && t === "I") {
				guidance.push("Class III → I: OH/NH shift strongly upfield or reappear (were absent in protic solvent). Other protons shift slightly upfield.");
			} else if (f === "II" && t === "III") {
				guidance.push("Class II → III: exchangeable protons shift or broaden, may vanish via H/D exchange. Non-exchangeable CH shifts are small.");
			} else if (f === "III" && t === "II") {
				guidance.push("Class III → II: exchangeable protons shift upfield and may sharpen. Non-exchangeable CH shifts are small.");
			} else if (f === "I" && t === "II") {
				guidance.push("Class I → II: OH/NH shift downfield. Aliphatic and aromatic CH largely unaffected.");
			} else if (f === "II" && t === "I") {
				guidance.push("Class II → I: OH/NH shift upfield. Aliphatic and aromatic CH largely unaffected.");
			}
			if (predAro || actAro) {
				if (f !== t)
					guidance.push("Aromatic solvent involved -- ring-current anisotropy causes additional shifts, especially for aromatic protons.");
			}
			if (solvent === "D2O" || predSolvent === "D2O" || solvent === "CD3OD" || predSolvent === "CD3OD") {
				guidance.push("Protic deuterated solvent -- OH, NH, NH₂ protons undergo H/D exchange and will most likely not be observed.");
			}
		}
		return { predCls, actCls, exchangeable: summarize(exchangeable), aromatic: summarize(aromatic), alphaHetero: summarize(alphaHetero), aliphatic: summarize(aliphatic), guidance };
	}, [predSolvent, solvent]);

	// ── Matching pipeline ──────────────────────────────────────────────────────
	const majorMatchRes = useMemo(() =>
		compData.map((c, i) => ({ ...c, origIdx: i })).filter(c => c.isProduct)
			.map(c => ({ ...c, matchResult: matchProtonCompound(c.predicted, matchablePeaks, tol, true) })),
		[compData, matchablePeaks, tol]);

	const resolvedImpurities = useMemo(() => {
		return solvData.impurities.filter(ki => enabledImp[ki.name]).map(ki => {
			const keys = [ki.name.toLowerCase(), ...(COMPOUND_ALIASES[ki.name] || [])];
			const annotated = allPeaks.filter(p => {
				const fields = [p.annotation, p.type, p.compound, p.flagStr].filter(Boolean).map(s => s.toLowerCase());
				return fields.some(f => keys.some(k => f.includes(k)));
			});
			if (annotated.length > 0) {
				return {
					name: ki.name, refPpm: ki.ppm, fromAnnotation: true,
					entries: annotated.map(p => ({ refPpm: p.ppm, matchedPeak: p, hasMatch: true, isAnnotated: true }))
				};
			}
			const entries = ki.ppm.map(refPpm => {
				let best = null, bestD = Infinity;
				for (const p of matchablePeaks) { const d = Math.abs(p.ppm - refPpm); if (d <= IMP_TOL && d < bestD) { best = p; bestD = d; } }
				return { refPpm, matchedPeak: best, hasMatch: best !== null, isAnnotated: false };
			});
			return { name: ki.name, refPpm: ki.ppm, fromAnnotation: false, entries };
		});
	}, [solvData, enabledImp, allPeaks, matchablePeaks]);

	const minorMatchRes = useMemo(() =>
		compData.map((c, i) => ({ ...c, origIdx: i })).filter(c => !c.isProduct)
			.map(c => ({ ...c, matchResult: matchProtonCompound(c.predicted, matchablePeaks, tol, false) })),
		[compData, matchablePeaks, tol]);

	const matchRes = useMemo(() =>
		[...majorMatchRes, ...minorMatchRes].sort((a, b) => a.origIdx - b.origIdx),
		[majorMatchRes, minorMatchRes]);

	const allAssigned = useMemo(() => {
		const a = new Set();
		majorMatchRes.forEach(c => c.matchResult.flat.filter(m => m.matched && m.match).forEach(m => a.add(m.match.ppm)));
		resolvedImpurities.forEach(ki => ki.entries.forEach(e => { if (e.matchedPeak) a.add(e.matchedPeak.ppm); }));
		minorMatchRes.forEach(c => c.matchResult.flat.filter(m => m.matched && m.match).forEach(m => a.add(m.match.ppm)));
		return a;
	}, [majorMatchRes, resolvedImpurities, minorMatchRes]);

	// ── Spectrum curves ────────────────────────────────────────────────────────
	const ppmSpan = viewRange[1] - viewRange[0];
	const nPts = Math.min(8000, Math.max(3000, Math.round(60000 / Math.max(ppmSpan, 1))));
	const compCurve = useMemo(() => compPeaks.length ? genSpectrum(compPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [compPeaks, viewRange, linewidth, nPts]);
	const c13SatCurve = useMemo(() => showC13Sats && c13SatPeaks.length ? genSpectrum(c13SatPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [showC13Sats, c13SatPeaks, viewRange, linewidth, nPts]);
	const solvCurve = useMemo(() => solvPeaks.length ? genSpectrum(solvPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [solvPeaks, viewRange, linewidth, nPts]);
	const impCurve = useMemo(() => impPeaks.length ? genSpectrum(impPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [impPeaks, viewRange, linewidth, nPts]);
	const secSolvCurve = useMemo(() => secondarySolvPeaks.length ? genSpectrum(secondarySolvPeaks, viewRange[0], viewRange[1], nPts, linewidth) : [], [secondarySolvPeaks, viewRange, linewidth, nPts]);

	const compMaxY = useMemo(() => { const v = compCurve.map(p => p[1]); return v.length ? Math.max(...v, 1) : 1; }, [compCurve]);
	const solvMaxY = Math.max(...solvCurve.map(p => p[1]), ...impCurve.map(p => p[1]), 0);
	const displayMaxY = scaleToComp ? compMaxY * 1.12 : Math.max(compMaxY, solvMaxY, 1) * 1.1;

	const unassigned = useMemo(() => compPeaks.filter(p => !allAssigned.has(p.ppm)), [compPeaks, allAssigned]);
	const groups = useMemo(() => groupPeaks(compPeaks, 0.08), [compPeaks]);

	// ── Helpers ────────────────────────────────────────────────────────────────
	const addComp = () => setCompounds(c => [...c, { name: `Compound ${c.length + 1}`, text: "", color: COLORS[c.length % COLORS.length], isProduct: false }]);
	const rmComp = i => setCompounds(c => c.filter((_, j) => j !== i));
	const updComp = (i, k, v) => setCompounds(c => c.map((x, j) => j === i ? { ...x, [k]: v } : x));
	const isCompVisible = (name) => visibleComps[name] !== false;
	const toggleCompVis = (name) => setVisibleComps(v => ({ ...v, [name]: v[name] === false ? true : false }));

	// ── Peak stats ─────────────────────────────────────────────────────────────
	const peakStats = useMemo(() => ({
		total: allPeaks.length, compound: compPeaks.length, solvent: solvPeaks.length,
		impurity: impPeaks.length, artifact: artifactPeaks.length, c13sat: c13SatPeaks.length,
		weak: allPeaks.filter(p => p.isWeak && !p.isC13Sat).length,
		hidden: allPeaks.filter(p => p.isHidden).length,
		secSolv: secondarySolvPeaks.length,
	}), [allPeaks, compPeaks, solvPeaks, impPeaks, artifactPeaks, c13SatPeaks, secondarySolvPeaks]);

	// ── SVG layout ─────────────────────────────────────────────────────────────
	const numC = compData.filter(c => c.predicted.length > 0).length;
	const mzH = Math.max(50, numC * 24 + 20);
	const svgW = 900, svgH = 360 + mzH;
	const maxLabelWidth = useMemo(() => {
		const names = compData.filter(c => c.predicted.length > 0).map(c => c.name);
		if (!names.length) return 20;
		return Math.min(63, Math.max(...names.map(n => n.length)) * 6.5 + 8);
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
		if (hi - lo > 0.05) { setZoomStack(s => [...s, [...viewRange]]); setViewRange([Math.max(lo, fullRange[0]), Math.min(hi, fullRange[1])]); }
		setDragging(false); setDragStart(null); setDragEnd(null);
	}, [dragging, dragStart, dragEnd, viewRange, fullRange]);
	const resetZoom = () => { setViewRange([...fullRange]); setZoomStack([]); };
	const zoomBack = () => { if (zoomStack.length) { setViewRange(zoomStack[zoomStack.length - 1]); setZoomStack(s => s.slice(0, -1)); } };
	const isZoomed = viewRange[0] !== fullRange[0] || viewRange[1] !== fullRange[1];
	const dragLo = dragStart !== null && dragEnd !== null ? Math.min(dragStart, dragEnd) : 0;
	const dragHi = dragStart !== null && dragEnd !== null ? Math.max(dragStart, dragEnd) : 0;

	const st = { bg: "#0f172a", card: "#1e293b", bdr: "#334155", tx: "#e2e8f0", mt: "#94a3b8", btn: "#3b82f6" };
	const Btn = ({ active, children, onClick, style: s2 }) => (
		<button onClick={onClick} style={{ padding: "8px 16px", borderRadius: 8, border: "none", cursor: "pointer", background: active ? st.btn : st.card, color: active ? "#fff" : st.mt, fontWeight: 600, fontSize: 13, ...(s2 || {}) }}>{children}</button>
	);
	const Sel = ({ value, onChange, children, style: s2 }) => (
		<select value={value} onChange={onChange} style={{ background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 8px", fontSize: 13, ...(s2 || {}) }}>{children}</select>
	);

	return (
		<div style={{ fontFamily: "system-ui,sans-serif", background: st.bg, color: st.tx, minHeight: "100vh", padding: 16 }}>
			<h1 style={{ fontSize: 22, fontWeight: 700, margin: 0 }}>¹H NMR Spectrum Analyzer</h1>
			<p style={{ color: st.mt, fontSize: 13, marginBottom: 16 }}>Paste MNOVA peak tables & ChemDraw predictions and auto-compare</p>

			<div style={{ display: "flex", gap: 8, marginBottom: 16, flexWrap: "wrap", alignItems: "center" }}>
				{TABS.map(t => <Btn key={t} active={tab === t} onClick={() => setTab(t)}>{t === "input" ? "📋 Input" : t === "spectrum" ? "📊 Spectrum" : "🔬 Analysis"}</Btn>)}
				{allPeaks.length > 0 && (
					<span style={{ color: st.mt, fontSize: 12, marginLeft: 8 }}>
						{peakStats.total} peaks · {peakStats.compound} cmpd · {peakStats.solvent} solv · {peakStats.impurity} imp
						{peakStats.c13sat > 0 && ` · ${peakStats.c13sat} ¹³C sat`}
						{peakStats.hidden > 0 && ` · ${peakStats.hidden} hidden`}
						{peakStats.weak > 0 && ` · ${peakStats.weak} weak`}
					</span>
				)}
			</div>

			{/* ===== INPUT TAB ===== */}
			{tab === "input" && (
				<div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
					<div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
						<h3 style={{ fontSize: 15, fontWeight: 700, marginBottom: 8 }}>MNOVA ¹H Peak Table</h3>
						<p style={{ fontSize: 12, color: st.mt, marginBottom: 8 }}>Paste tab-separated peak table (index, ppm, intensity, width, area, type, flags, impurity/compound, annotation)</p>
						<textarea value={mnovaText} onChange={e => setMnovaText(e.target.value)}
							placeholder={"1\t7.26\t1000.0\t1.50\t5000.00\tSolvent\tNone\tCDCl3\t\n2\t7.36\t76.2\t1.03\t1373.62\tCompound\tHidden\t\t\n3\t1.67\t47.5\t2.02\t1527.12\tCompound\tWeak + C13_Sat\t\t"}
							style={{ width: "100%", height: 150, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 8, padding: 10, fontSize: 12, fontFamily: "monospace", resize: "vertical", boxSizing: "border-box" }} />
						{allPeaks.length > 0 && (
							<div style={{ fontSize: 12, color: "#4ade80", marginTop: 6 }}>
								✓ {peakStats.total} peaks parsed
								<span style={{ color: st.mt, marginLeft: 8 }}>
									({peakStats.compound} compound, {peakStats.solvent} solvent, {peakStats.impurity} impurity, {peakStats.artifact} artifact)
								</span>
								{peakStats.c13sat > 0 && <span style={{ color: "#f97316", marginLeft: 8 }}>⚡ {peakStats.c13sat} ¹³C satellites</span>}
								{peakStats.hidden > 0 && <span style={{ color: "#64748b", marginLeft: 8 }}>👁 {peakStats.hidden} hidden</span>}
							</div>
						)}
					</div>

					<div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
						<div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 12 }}>
							<h3 style={{ fontSize: 15, fontWeight: 700, margin: 0 }}>ChemDraw ¹H Predictions</h3>
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
								<textarea value={c.text} onChange={e => updComp(i, "text", e.target.value)} placeholder={"Paste ChemDraw ¹H NMR Prediction here...\nCH 7.95 7.62 benzylidenimin\nCH2 3.38;3.135 1.37 methylene"}
									style={{ width: "100%", height: 100, background: st.card, color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: 8, fontSize: 11, fontFamily: "monospace", resize: "vertical", boxSizing: "border-box" }} />
								{compData[i]?.predicted.length > 0 && (
									<p style={{ fontSize: 12, color: "#4ade80", marginTop: 4 }}>
										✓ {compData[i].predicted.length} nodes ({compData[i].predicted.reduce((s, n) => s + n.hCount, 0)}H total)
										— {compData[i].predicted.map(n => `${n.type} ${n.shifts.map(s => s.toFixed(2)).join(";")}`).join(", ")}
									</p>
								)}
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
									<span style={{ fontWeight: 600 }}>Expected ¹H shift (actual − predicted):</span>
								</div>
								<div style={{ display: "grid", gridTemplateColumns: "repeat(4, 1fr)", gap: 8, marginBottom: 8 }}>
									{[
										{ label: "Exchangeable (OH/NH)", data: solventPrediction.exchangeable, color: "#ef4444" },
										{ label: "Aromatic (>6 ppm)", data: solventPrediction.aromatic, color: "#f59e0b" },
										{ label: "α-Heteroatom (3-6 ppm)", data: solventPrediction.alphaHetero, color: "#a78bfa" },
										{ label: "Aliphatic (<3 ppm)", data: solventPrediction.aliphatic, color: "#22c55e" },
									].map(({ label, data, color }) => (
										<div key={label} style={{ background: "#0f172a", borderRadius: 6, padding: 8, borderLeft: `3px solid ${color}` }}>
											<div style={{ fontSize: 11, color: st.mt, marginBottom: 4 }}>{label}</div>
											{data ? (
												<>
													<div style={{ fontSize: 14, fontWeight: 700, color }}>
														{data.avg >= 0 ? "+" : ""}{data.avg.toFixed(2)} ppm
													</div>
													<div style={{ fontSize: 10, color: st.mt }}>
														range: {data.min >= 0 ? "+" : ""}{data.min.toFixed(2)} to {data.max >= 0 ? "+" : ""}{data.max.toFixed(2)}
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
									Computed from {Object.keys(protonImpurities.compounds).length} reference compounds. Consider adjusting tolerance to account for solvent effects.
								</p>
							</div>
						)}
						<p style={{ fontSize: 12, color: st.mt, marginBottom: 8 }}>Known ¹H impurities in {solvData.label} ({solvData.impurities.length} compounds):</p>
						<div style={{ display: "flex", flexWrap: "wrap", gap: 6, maxHeight: 240, overflowY: "auto", padding: 4 }}>
							{solvData.impurities.map((ki, i) => {
								const en = !!enabledImp[ki.name];
								const keys = [ki.name.toLowerCase(), ...(COMPOUND_ALIASES[ki.name] || [])];
								const annotated = allPeaks.some(p => {
									const fields = [p.annotation, p.type, p.compound, p.flagStr].filter(Boolean).map(s => s.toLowerCase());
									return fields.some(f => keys.some(k => f.includes(k)));
								});
								return (
									<button key={i} onClick={() => setEnabledImp(prev => ({ ...prev, [ki.name]: !en }))}
										style={{ padding: "5px 12px", borderRadius: 8, border: en ? "2px solid #d946ef" : `1px solid ${st.bdr}`, background: en ? "#2e1065" : "#0f172a", color: en ? "#e9d5ff" : st.mt, cursor: "pointer", fontSize: 12, fontWeight: en ? 600 : 400 }}>
										{ki.name.charAt(0).toUpperCase() + ki.name.slice(1)} {ki.ppm.length > 0 && <span style={{ opacity: 0.6 }}>({ki.ppm.slice(0, 3).map(p => p.toFixed(2)).join(", ")}{ki.ppm.length > 3 ? ", …" : ""})</span>}
										{annotated && <span style={{ color: "#4ade80", marginLeft: 4 }}>*</span>}
									</button>
								);
							})}
						</div>
						<p style={{ fontSize: 11, color: st.mt, marginTop: 6 }}>
							* = found in MNOVA annotations (annotated ppm used directly).
							Non-annotated impurities are matched within ±{IMP_TOL} ppm of the table reference shift.
							Shifts sourced from experimental ¹H data table.
						</p>
					</div>

					<div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
						<h3 style={{ fontSize: 15, fontWeight: 700, marginBottom: 10 }}>Settings</h3>
						<div style={{ display: "flex", gap: 16, flexWrap: "wrap", alignItems: "center", fontSize: 13 }}>
							<label>Tolerance: <input type="number" value={tol} onChange={e => setTol(parseFloat(e.target.value) || 0.3)} min={0.01} max={2} step={0.05}
								style={{ width: 55, marginLeft: 4, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} /> ppm</label>
							<label>Linewidth: <input type="number" value={linewidth} onChange={e => setLinewidth(parseFloat(e.target.value) || 0.003)} min={0.0005} max={0.5} step={0.0005}
								style={{ width: 65, marginLeft: 4, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} /> ppm</label>
							<label>Range: <input type="number" value={fullRange[0]} onChange={e => { const v = parseFloat(e.target.value) || -0.5; setFullRange([v, fullRange[1]]); setViewRange([v, fullRange[1]]); setZoomStack([]); }}
								style={{ width: 50, marginLeft: 4, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} />–
								<input type="number" value={fullRange[1]} onChange={e => { const v = parseFloat(e.target.value) || 14; setFullRange([fullRange[0], v]); setViewRange([fullRange[0], v]); setZoomStack([]); }}
									style={{ width: 50, marginLeft: 2, background: "#0f172a", color: st.tx, border: `1px solid ${st.bdr}`, borderRadius: 6, padding: "4px 6px", fontSize: 13 }} /></label>
							<label><input type="checkbox" checked={showSolvent} onChange={e => setShowSolvent(e.target.checked)} /> Solvent</label>
							<label><input type="checkbox" checked={showC13Sats} onChange={e => setShowC13Sats(e.target.checked)} /> ¹³C satellites</label>
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
								{isZoomed && <span style={{ fontSize: 12, color: st.mt }}>{viewRange[0].toFixed(2)}–{viewRange[1].toFixed(2)} ppm</span>}
							</div>

							{matchRes.filter(c => c.predicted.length > 0).length > 0 && (
								<div style={{ display: "flex", gap: 6, marginBottom: 10, flexWrap: "wrap", alignItems: "center" }}>
									<span style={{ fontSize: 12, color: st.mt }}>Show overlays:</span>
									{matchRes.filter(c => c.predicted.length > 0).map((c, i) => {
										const vis = isCompVisible(c.name);
										return (
											<button key={i} onClick={() => toggleCompVis(c.name)}
												style={{
													padding: "4px 12px", borderRadius: 8, border: vis ? `2px solid ${c.color}` : `1px solid ${st.bdr}`,
													background: vis ? st.card : "#0f172a", color: vis ? c.color : st.mt,
													cursor: "pointer", fontSize: 12, fontWeight: vis ? 700 : 400, opacity: vis ? 1 : 0.5
												}}>
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

								{/* Grid lines */}
								{(() => {
									const span = viewRange[1] - viewRange[0];
									const step = span > 10 ? 2 : span > 5 ? 1 : span > 2 ? 0.5 : span > 0.5 ? 0.1 : 0.05;
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

								{/* Solvent shading */}
								{showSolvent && nmrSolvPeaks.map((p, i) => (
									<rect key={`sv${i}`} x={xSc(p.ppm) - 3} y={pad.t} width={6} height={plotH} fill="#475569" opacity={0.1} />
								))}

								{/* Impurity reference lines */}
								{(() => {
									const items = [];
									resolvedImpurities.forEach((ki, ki_i) => ki.entries.forEach((e, e_i) => {
										items.push({ x: xSc(e.refPpm), ki_i, e_i, name: ki.name, hasMatch: e.hasMatch, isAnnotated: e.isAnnotated, matchedPeak: e.matchedPeak, refPpm: e.refPpm });
									}));
									items.sort((a, b) => a.x - b.x);
									const placed = [];
									for (const lb of items) {
										let yOff = 0, show = false;
										while (yOff <= 80) {
											const y = pad.t + 14 + yOff;
											if (!placed.some(p => { const minD = (lb.name.length + p.name.length) * 2.2 + 8; return Math.abs(p.x - lb.x) < minD && Math.abs(p.y - y) < 10; })) { show = true; break; }
											yOff += 11;
										}
										placed.push({ ...lb, y: pad.t + 14 + yOff, show });
									}
									return placed.map(lb => (
										<g key={`ki${lb.ki_i}${lb.e_i}`}>
											<line x1={lb.x} x2={lb.x} y1={pad.t} y2={pad.t + plotH}
												stroke={lb.hasMatch || lb.isAnnotated ? "#d946ef" : "#94a3b8"} strokeWidth={1} strokeDasharray="5,5" opacity={lb.isAnnotated ? 0.6 : lb.hasMatch ? 0.45 : 0.35} />
											{lb.hasMatch && !lb.isAnnotated && lb.matchedPeak && (
												<line x1={lb.x} x2={xSc(lb.matchedPeak.ppm)}
													y1={pad.t} y2={ySc(Math.min(evalAt(lb.matchedPeak.ppm, compPeaks, linewidth), displayMaxY))}
													stroke="#a855f7" strokeWidth={1.5} opacity={0.45} strokeDasharray="3,3" />
											)}
											{lb.show && (
												<text x={lb.x} y={lb.y} textAnchor="middle" fill={lb.hasMatch ? "#d946ef" : "#94a3b8"} fontSize={8} opacity={lb.hasMatch ? 0.75 : 0.55}>
													{lb.name}{lb.isAnnotated && <tspan fill="#4ade80"> ●</tspan>}
												</text>
											)}
										</g>
									));
								})()}

								{/* Curves */}
								{showSolvent && solvCurve.length > 0 && <path d={makePath(solvCurve)} fill="none" stroke="#64748b" strokeWidth={1} opacity={0.5} />}
								{showSolvent && secSolvCurve.length > 0 && <path d={makePath(secSolvCurve)} fill="none" stroke="#f59e0b" strokeWidth={1} opacity={0.5} />}
								{impCurve.length > 0 && <path d={makePath(impCurve)} fill="none" stroke="#f97316" strokeWidth={1} opacity={0.6} />}
								{showC13Sats && c13SatCurve.length > 0 && <path d={makePath(c13SatCurve)} fill="none" stroke="#f43f5e" strokeWidth={0.8} opacity={0.4} strokeDasharray="3,3" />}
								{compCurve.length > 0 && <path d={makePath(compCurve)} fill="none" stroke="#3b82f6" strokeWidth={1.5} />}

								{/* Compound prediction marker zone */}
								<rect x={pad.l} y={4} width={plotW} height={mzH - 8} fill="#1e293b" rx={6} opacity={0.5} />
								{matchRes.filter(c => c.predicted.length > 0).map((comp, ci) => {
									const vis = isCompVisible(comp.name);
									const yB = 18 + ci * 24;
									return (
										<g key={`cg${ci}`} opacity={vis ? 1 : 0.15}>
											<text x={pad.l - 4} y={yB + 4} textAnchor="end" fill={comp.color} fontSize={10} fontWeight="600">{comp.name}</text>
											<line x1={pad.l} x2={pad.l + plotW} y1={yB} y2={yB} stroke={comp.color} strokeWidth={0.4} opacity={0.2} />
											{comp.matchResult.flat.map((m, mi) => {
												const x = xSc(m.shift);
												if (x < pad.l - 10 || x > pad.l + plotW + 10) return null;
												return (
													<g key={`m${ci}${mi}`}>
														<circle cx={x} cy={yB} r={4} fill={m.matched ? comp.color : "transparent"} stroke={comp.color} strokeWidth={1.5} opacity={0.9} />
														<text x={x} y={yB - 7} textAnchor="middle" fill={comp.color} fontSize={7.5} opacity={0.8}>{m.shift.toFixed(2)}</text>
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

								{/* Peak ppm labels — SKIP hidden peaks */}
								{(() => {
									const labels = compPeaks
										.filter(p => p.ppm >= viewRange[0] && p.ppm <= viewRange[1])
										.filter(p => !p.isHidden) // hidden peaks: drawn but NOT labeled
										.filter(p => p.intensity > compMaxY * 0.08)
										.map(p => ({ ppm: p.ppm, x: xSc(p.ppm), baseY: ySc(Math.min(evalAt(p.ppm, compPeaks, linewidth), displayMaxY)) - 5, text: p.ppm.toFixed(2), isWeak: p.isWeak }))
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
										<text key={`pl${i}`} x={lb.x} y={lb.y} textAnchor="middle"
											fill={lb.isWeak ? "#fbbf24" : "#93c5fd"} fontSize={8} opacity={lb.isWeak ? 0.5 : 0.65}>
											{lb.text}
										</text>
									));
								})()}
							</svg>

							<div style={{ display: "flex", gap: 14, flexWrap: "wrap", marginTop: 10, fontSize: 12 }}>
								<span><span style={{ color: "#3b82f6" }}>━━</span> Compound</span>
								{showSolvent && <span><span style={{ color: "#64748b" }}>━━</span> Solvent</span>}
								{impPeaks.length > 0 && <span><span style={{ color: "#f97316" }}>━━</span> Impurity</span>}
								{showC13Sats && c13SatPeaks.length > 0 && <span><span style={{ color: "#f43f5e" }}>┄┄</span> ¹³C satellites</span>}
								{secondarySolvPeaks.length > 0 && <span><span style={{ color: "#f59e0b" }}>━━</span> 2° solvent</span>}
								<span><span style={{ color: "#d946ef" }}>┄┄</span> Ref impurities</span>
							</div>
						</div>
					)}
				</div>
			)}

			{/* ===== ANALYSIS TAB ===== */}
			{tab === "analysis" && (
				<div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
					{/* Compound match tables (grouped by node) */}
					{matchRes.filter(c => c.predicted.length > 0).map((comp, ci) => {
						const gm = comp.matchResult.grouped;
						const totalShifts = gm.reduce((s, g) => s + g.results.length, 0);
						const matchedShifts = gm.reduce((s, g) => s + g.results.filter(r => r.matched).length, 0);
						const pct = totalShifts > 0 ? Math.round(matchedShifts / totalShifts * 100) : 0;
						const avgScore = gm.filter(g => g.anyMatched).length > 0
							? gm.filter(g => g.anyMatched).reduce((s, g) => s + scoreProtonMatch(g, maxArea, comp.isProduct), 0) / gm.filter(g => g.anyMatched).length : 0;
						return (
							<div key={ci} style={{ background: st.card, borderRadius: 12, padding: 16, borderLeft: `4px solid ${comp.color}` }}>
								<div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 4, flexWrap: "wrap", gap: 8 }}>
									<h3 style={{ fontSize: 16, fontWeight: 700, color: comp.color, margin: 0 }}>
										{comp.name} {comp.isProduct && <span style={{ fontSize: 11, color: st.mt }}>(major)</span>}
									</h3>
									<div style={{ display: "flex", gap: 8 }}>
										<span style={{ padding: "4px 12px", borderRadius: 20, fontSize: 12, fontWeight: 700, background: pct >= 70 ? "#166534" : pct >= 40 ? "#854d0e" : "#7f1d1d", color: "#fff" }}>
											{matchedShifts}/{totalShifts} shifts ({pct}%)
										</span>
										<span style={{ padding: "4px 12px", borderRadius: 20, fontSize: 12, fontWeight: 600, background: avgScore > 0.7 ? "#166534" : avgScore > 0.4 ? "#854d0e" : "#7f1d1d", color: "#fff" }}>
											Score: {(avgScore * 100).toFixed(0)}%
										</span>
									</div>
								</div>
								<table style={{ width: "100%", fontSize: 12, borderCollapse: "collapse", marginTop: 8 }}>
									<thead><tr style={{ color: st.mt, textAlign: "left", borderBottom: `1px solid ${st.bdr}` }}>
										<th style={{ padding: "6px 8px" }}>Type</th>
										<th style={{ padding: "6px 8px" }}>H</th>
										<th style={{ padding: "6px 8px" }}>Predicted</th>
										<th style={{ padding: "6px 8px" }}>Matched</th>
										<th style={{ padding: "6px 8px" }}>Δ avg</th>
										<th style={{ padding: "6px 8px" }}>Score</th>
										<th style={{ padding: "6px 8px" }}>Status</th>
									</tr></thead>
									<tbody>{gm.map((nm, j) => {
										const sc = scoreProtonMatch(nm, maxArea, comp.isProduct);
										const predStr = nm.node.shifts.map(s => s.toFixed(2)).join("; ");
										const matchStr = nm.results.filter(r => r.matched && r.match).map(r => r.match.ppm.toFixed(2)).join("; ") || "—";
										const hiddenMatches = nm.results.filter(r => r.matched && r.match && r.match.isHidden);
										return (
											<tr key={j} style={{ borderBottom: "1px solid #1e293b" }}>
												<td style={{ padding: "6px 8px", fontWeight: 600 }}>{nm.node.type}</td>
												<td style={{ padding: "6px 8px" }}>{nm.node.hCount}</td>
												<td style={{ padding: "6px 8px" }}>{predStr}</td>
												<td style={{ padding: "6px 8px" }}>
													{matchStr}
													{hiddenMatches.length > 0 && <span style={{ color: "#64748b", fontSize: 10, marginLeft: 4 }}>(hidden)</span>}
												</td>
												<td style={{ padding: "6px 8px" }}>{nm.avgDelta < 100 ? nm.avgDelta.toFixed(3) : "—"}</td>
												<td style={{ padding: "6px 8px" }}>
													{nm.anyMatched && <div style={{ display: "flex", alignItems: "center", gap: 4 }}>
														<div style={{ width: 50, height: 8, background: "#0f172a", borderRadius: 4, overflow: "hidden" }}>
															<div style={{ width: `${sc * 100}%`, height: "100%", background: sc > 0.7 ? "#22c55e" : sc > 0.4 ? "#f59e0b" : "#ef4444", borderRadius: 4 }} />
														</div>
														<span style={{ fontSize: 10, color: st.mt }}>{(sc * 100).toFixed(0)}</span>
													</div>}
												</td>
												<td style={{ padding: "6px 8px" }}>
													<span style={{
														padding: "2px 10px", borderRadius: 12, fontSize: 11, fontWeight: 600,
														background: nm.allMatched ? (nm.avgDelta < 0.15 ? "#166534" : "#854d0e") : nm.anyMatched ? "#854d0e" : "#7f1d1d", color: "#fff"
													}}>
														{nm.allMatched ? (nm.avgDelta < 0.15 ? "✓ Good" : "~ Approx") : nm.anyMatched ? "◐ Partial" : "✕ No match"}
													</span>
												</td>
											</tr>
										);
									})}</tbody>
								</table>
								<p style={{ fontSize: 11, color: st.mt, marginTop: 6 }}>
									Diastereotopic CH₂ protons show two predicted shifts. Hidden peaks are matched but noted — they may have unreliable intensities.
								</p>
							</div>
						);
					})}

					{/* Flag summary */}
					{(peakStats.c13sat > 0 || peakStats.hidden > 0 || peakStats.weak > 0) && (
						<div style={{ background: st.card, borderRadius: 12, padding: 16, borderLeft: "4px solid #f97316" }}>
							<h3 style={{ fontSize: 16, fontWeight: 700, marginBottom: 10, color: "#f97316" }}>Flag Summary</h3>
							<div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(200px, 1fr))", gap: 8 }}>
								{peakStats.c13sat > 0 && (
									<div style={{ background: "#0f172a", borderRadius: 8, padding: 10, borderLeft: "3px solid #f43f5e" }}>
										<div style={{ fontWeight: 700, color: "#f43f5e" }}>{peakStats.c13sat} ¹³C Satellites</div>
										<div style={{ fontSize: 11, color: st.mt }}>Excluded from matching. Small doublets flanking strong peaks (~0.55% intensity).</div>
									</div>
								)}
								{peakStats.hidden > 0 && (
									<div style={{ background: "#0f172a", borderRadius: 8, padding: 10, borderLeft: "3px solid #64748b" }}>
										<div style={{ fontWeight: 700, color: "#64748b" }}>{peakStats.hidden} Hidden</div>
										<div style={{ fontSize: 11, color: st.mt }}>Drawn but not labeled. Peaks obscured by overlap — included in matching but integration unreliable.</div>
									</div>
								)}
								{peakStats.weak > 0 && (
									<div style={{ background: "#0f172a", borderRadius: 8, padding: 10, borderLeft: "3px solid #fbbf24" }}>
										<div style={{ fontWeight: 700, color: "#fbbf24" }}>{peakStats.weak} Weak</div>
										<div style={{ fontSize: 11, color: st.mt }}>Low-confidence peaks. Labeled in yellow. Included in matching but may be noise.</div>
									</div>
								)}
							</div>
						</div>
					)}

					{/* Unassigned */}
					<div style={{ background: st.card, borderRadius: 12, padding: 16, borderLeft: "4px solid #64748b" }}>
						<h3 style={{ fontSize: 16, fontWeight: 700, color: "#94a3b8", marginBottom: 10 }}>Unassigned Compound Peaks ({unassigned.length})</h3>
						{unassigned.length === 0 ? <p style={{ color: "#4ade80", fontSize: 13 }}>All compound peaks accounted for!</p> : (
							<div style={{ display: "flex", flexWrap: "wrap", gap: 8 }}>{unassigned.map((p, i) => (
								<div key={i} style={{
									background: "#0f172a", borderRadius: 8, padding: "6px 12px", fontSize: 12,
									border: p.isHidden ? "1px solid #475569" : p.isWeak ? "1px solid #854d0e" : "1px solid transparent"
								}}>
									<span style={{ fontWeight: 700, color: p.isHidden ? "#64748b" : p.isWeak ? "#fbbf24" : "#fbbf24" }}>{p.ppm.toFixed(2)}</span>
									<span style={{ color: st.mt, marginLeft: 4 }}>I:{p.intensity.toFixed(0)}</span>
									<span style={{ color: st.mt, marginLeft: 4 }}>A:{p.area.toFixed(0)}</span>
									{p.isHidden && <span style={{ color: "#64748b", marginLeft: 4, fontSize: 10 }}>hidden</span>}
									{p.isWeak && <span style={{ color: "#fbbf24", marginLeft: 4, fontSize: 10 }}>weak</span>}
								</div>
							))}</div>
						)}

						{/* Mini spectrum: assigned vs unassigned */}
						{compPeaks.length > 0 && (() => {
							const W = 700, H = 120, P = { l: 10, r: 10, t: 8, b: 24 };
							const pw = W - P.l - P.r, ph = H - P.t - P.b;
							const allP = [...compPeaks, ...solvPeaks];
							const r = [Math.min(...allP.map(p => p.ppm)) - 0.5, Math.max(...allP.map(p => p.ppm)) + 0.5];
							const n = 1500, gm = linewidth;
							const asgn = compPeaks.filter(p => allAssigned.has(p.ppm));
							const unasgn = compPeaks.filter(p => !allAssigned.has(p.ppm));
							const aC = asgn.length ? genSpectrum(asgn, r[0], r[1], n, gm) : [];
							const uC = unasgn.length ? genSpectrum(unasgn, r[0], r[1], n, gm) : [];
							const sC = solvPeaks.length ? genSpectrum(solvPeaks, r[0], r[1], n, gm) : [];
							const maxY = Math.max(...aC.map(p => p[1]), ...uC.map(p => p[1]), ...sC.map(p => p[1]), 1) * 1.1;
							const xS = v => P.l + (r[1] - v) / (r[1] - r[0]) * pw;
							const yS = v => P.t + ph - (v / maxY) * ph;
							const mp = c => { let d = '', w = false; for (const p of c) { if (p[1] > maxY * .002) { d += (w ? 'L' : 'M') + xS(p[0]).toFixed(1) + ',' + yS(p[1]).toFixed(1) + ' '; w = true } else w = false } return d };
							const span = r[1] - r[0], step = span > 10 ? 2 : span > 5 ? 1 : 0.5;
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
									</g>
								</svg>
							);
						})()}
					</div>

					{/* Peak groups */}
					<div style={{ background: st.card, borderRadius: 12, padding: 16 }}>
						<h3 style={{ fontSize: 16, fontWeight: 700, marginBottom: 10 }}>Peak Groups (gap ≤ 0.08 ppm)</h3>
						<div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(180px, 1fr))", gap: 8 }}>
							{groups.map((g, i) => (
								<div key={i} style={{ background: "#0f172a", borderRadius: 8, padding: 10, borderLeft: `3px solid ${GROUP_COLORS[i % GROUP_COLORS.length]}` }}>
									<div style={{ fontWeight: 700, fontSize: 16, color: GROUP_COLORS[i % GROUP_COLORS.length] }}>~{g.center.toFixed(2)}</div>
									<div style={{ fontSize: 11, color: st.mt }}>{g.peaks.length} peak{g.peaks.length > 1 ? "s" : ""} · {g.min.toFixed(2)}–{g.max.toFixed(2)}</div>
									<div style={{ fontSize: 11, color: st.mt }}>Area: {g.totalArea.toFixed(0)}</div>
								</div>
							))}
						</div>
					</div>
				</div>
			)}
		</div>
	);
}