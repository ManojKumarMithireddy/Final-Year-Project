import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Dna, FlaskConical, RotateCcw } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import toast from 'react-hot-toast';
import api from '../lib/api';
import GroverStepNavigator from './GroverStepNavigator';

// ── Pipeline step card ────────────────────────────────────────────────────────
function PipelineStep({ n, label, type, done, optional, locked, last, children }) {
  const isLocked = !!locked && !done;
  const T = {
    input:  { pill: 'bg-amber-500/10 border-amber-500/30 text-amber-400',  txt: '▶ INPUT'  },
    action: { pill: 'bg-sky-500/10  border-sky-500/30   text-sky-400',     txt: '⚡ ACTION' },
    output: { pill: 'bg-cyan-500/10 border-cyan-500/30  text-cyan-400',    txt: '◀ OUTPUT' },
  };
  const cfg = T[type] ?? T.input;
  const leftBorder = done     ? 'rgba(16,185,129,0.55)'
                   : isLocked ? 'rgba(51,65,85,0.35)'
                   :            'rgba(245,158,11,0.55)';
  return (
    <div className="flex gap-4 items-start">
      {/* Step indicator + vertical connector */}
      <div className="flex flex-col items-center w-10 shrink-0 self-stretch">
        <motion.div
          initial={false}
          animate={{
            boxShadow: done     ? '0 0 0 6px rgba(16,185,129,0.12)'
                     : isLocked ? 'none'
                     :            '0 0 0 6px rgba(245,158,11,0.12)',
          }}
          transition={{ duration: 0.4 }}
          className={`w-10 h-10 rounded-full border-2 flex items-center justify-center text-sm font-bold relative bg-slate-950 shrink-0 transition-colors duration-500 ${
            done     ? 'border-emerald-500 text-emerald-400'
          : isLocked ? 'border-slate-700   text-slate-600'
          :            'border-amber-500   text-amber-400'
          }`}
        >
          {done ? '✓' : n}
        </motion.div>
        {!last && (
          <div className={`w-0.5 flex-1 min-h-[24px] mt-1 transition-colors duration-500 ${
            done ? 'bg-gradient-to-b from-emerald-500/40 to-emerald-500/10'
                 : 'bg-gradient-to-b from-slate-600/30 to-transparent'
          }`} />
        )}
      </div>

      {/* Content card */}
      <motion.div
        initial={{ opacity: 0, y: 8 }}
        animate={{ opacity: isLocked ? 0.35 : 1, y: 0 }}
        transition={{ duration: 0.3 }}
        className={`flex-1 min-w-0 rounded-2xl border p-5 mb-6 ${
          done     ? 'border-emerald-900/50 bg-emerald-950/10'
        : isLocked ? 'border-slate-800/40  bg-slate-950/10 pointer-events-none'
        :            'border-slate-700/50  bg-slate-900/40'
        }`}
        style={{ borderLeftWidth: '3px', borderLeftColor: leftBorder }}
      >
        {/* Card header */}
        <div className="flex items-center gap-2 mb-4 flex-wrap">
          <span className={`inline-flex items-center px-2 py-0.5 rounded-full border text-[9px] font-black uppercase tracking-widest shrink-0 ${cfg.pill}`}>
            {cfg.txt}
          </span>
          <span className={`text-sm font-semibold ${
            done ? 'text-emerald-300/90' : isLocked ? 'text-slate-600' : 'text-slate-100'
          }`}>{label}</span>
          {optional && <span className="text-[10px] text-slate-600 uppercase tracking-widest ml-1">· optional</span>}
          <div className="ml-auto flex items-center gap-2">
            {isLocked  && <span className="text-[10px] text-slate-600 italic">↑ complete previous step first</span>}
            {!isLocked && !done && type !== 'output' && <span className="text-[10px] text-slate-500 uppercase tracking-widest">👆 requires action</span>}
            {done      && <span className="text-[10px] text-emerald-500/60 uppercase tracking-widest">✓ complete</span>}
          </div>
        </div>
        {children}
      </motion.div>
    </div>
  );
}

// ── Input / Output section divider ───────────────────────────────────────────
function SectionLabel({ type, label }) {
  const isInput = type === 'input';
  return (
    <div className={`flex items-center gap-3 ${isInput ? 'text-amber-400/70' : 'text-cyan-400/70'}`}>
      <span className={`inline-flex items-center px-2 py-0.5 rounded-full border text-[9px] font-black uppercase tracking-widest shrink-0 ${
        isInput
          ? 'bg-amber-500/10 border-amber-500/30 text-amber-400'
          : 'bg-cyan-500/10 border-cyan-500/30 text-cyan-400'
      }`}>{isInput ? '▶ INPUT' : '◀ OUTPUT'}</span>
      <span className="text-xs font-bold uppercase tracking-widest">{label}</span>
      <span className={`flex-1 h-px ${isInput ? 'bg-amber-500/15' : 'bg-cyan-500/15'}`} />
    </div>
  );
}

// ── Encoding helpers ──────────────────────────────────────────────────────────
const DEC = { '00': 'A', '01': 'C', 10: 'G', 11: 'T' };
const bitsToDna = (bits) => {
  let s = '';
  for (let i = 0; i < bits.length; i += 2) s += DEC[bits.substring(i, i + 2)] ?? '?';
  return s;
};
// ── Mini detection badge (for comparison row) ─────────────────────────────────
function MiniDetectionBadge({ result, confidence, label, active, onClick }) {
  const found = result === 'FOUND';
  return (
    <button onClick={onClick} className={`w-full text-left rounded-2xl border p-5 transition-all ${
      active ? 'ring-2 ring-offset-2 ring-offset-slate-900' : 'opacity-80 hover:opacity-100'
    } ${found
      ? `bg-red-900/30 border-red-500/50 ${active ? 'ring-red-500' : ''}`
      : `bg-emerald-900/20 border-emerald-500/40 ${active ? 'ring-emerald-500' : ''}`
    }`}>
      <div className="text-xs font-medium text-slate-400 mb-2 uppercase tracking-wider">{label}</div>
      <div className={`text-xl font-bold mb-1 ${found ? 'text-red-400' : 'text-emerald-400'}`}>
        {found ? '🧬 DETECTED' : '✅ NOT FOUND'}
      </div>
      <div className={`text-xs ${found ? 'text-red-300/80' : 'text-emerald-300/80'}`}>
        {found ? 'c.5266dupC mutation present' : 'No matching mutation found'}
      </div>
      <div className="mt-2 text-xs text-slate-400">
        Confidence: <span className="font-mono text-white">{(confidence * 100).toFixed(1)}%</span>
        {confidence >= 0.995 && (
          <span className="ml-1.5 text-emerald-400/80" title="Grover's algorithm achieves near-perfect amplification when the optimal iteration count lands the target amplitude at sin²(π/2) ≈ 1">
            ✦ optimal
          </span>
        )}
      </div>
      {active && <div className="mt-2 text-xs text-slate-500">↓ step detail below</div>}
    </button>
  );
}

// ── Compact node table ────────────────────────────────────────────────────────
function NodeTable({ nodes, targetBits }) {
  const targetRowRef = useRef(null);
  const scrollContainerRef = useRef(null);

  // Auto-scroll to the target (marker) row after results arrive
  useEffect(() => {
    if (!targetRowRef.current || !scrollContainerRef.current) return;
    const timer = setTimeout(() => {
      const container = scrollContainerRef.current;
      const row = targetRowRef.current;
      if (!container || !row) return;
      const containerTop = container.getBoundingClientRect().top;
      const rowTop = row.getBoundingClientRect().top;
      const offset = rowTop - containerTop - container.clientHeight / 2 + row.clientHeight / 2;
      container.scrollBy({ top: offset, behavior: 'smooth' });
    }, 400);
    return () => clearTimeout(timer);
  }, [nodes]);

  const targetNode = nodes.find((nd) => nd.is_target);

  return (
    <div className="border border-slate-800 rounded-xl overflow-hidden">
      <div className="px-4 py-2.5 bg-slate-900/80 border-b border-slate-800 flex items-center justify-between">
        <span className="text-xs font-medium text-slate-400 uppercase tracking-widest">Patient DNA Node Table</span>
        <div className="flex items-center gap-3">
          {targetNode && (
            <span className="flex items-center gap-1 text-xs text-amber-400 font-semibold">
              <span className="relative flex h-2 w-2">
                <span className="animate-ping absolute inline-flex h-full w-full rounded-full bg-amber-400 opacity-75" />
                <span className="relative inline-flex rounded-full h-2 w-2 bg-amber-400" />
              </span>
              Marker at node #{targetNode.position}
            </span>
          )}
          <span className="text-xs text-slate-500">{nodes.length} nodes shown</span>
        </div>
      </div>
      <div ref={scrollContainerRef} className="max-h-52 overflow-y-auto">
        <table className="w-full text-left text-xs">
          <thead className="sticky top-0 bg-slate-900 border-b border-slate-800 z-10">
            <tr>
              <th className="px-4 py-2 text-slate-500 font-medium">#</th>
              <th className="px-4 py-2 text-slate-500 font-medium">DNA</th>
              <th className="px-4 py-2 text-slate-500 font-medium">Bits</th>
              <th className="px-4 py-2 text-slate-500 font-medium">Match</th>
            </tr>
          </thead>
          <tbody>
            {nodes.map((nd) =>
              nd.is_target ? (
                <motion.tr
                  key={nd.position}
                  ref={targetRowRef}
                  initial={{ opacity: 0, backgroundColor: 'rgba(251,191,36,0)' }}
                  animate={{ opacity: 1, backgroundColor: ['rgba(251,191,36,0.05)', 'rgba(251,191,36,0.30)', 'rgba(251,191,36,0.18)'] }}
                  transition={{ duration: 1.2, ease: 'easeInOut', delay: 0.5 }}
                  className="border-b border-amber-500/40"
                  style={{ boxShadow: 'inset 3px 0 0 #f59e0b' }}
                >
                  <td className="px-4 py-1.5 font-mono text-amber-300/70 font-bold">{nd.position}</td>
                  <td className="px-4 py-1.5 font-mono tracking-widest text-amber-300 font-black text-sm">
                    {nd.dna}
                  </td>
                  <td className="px-4 py-1.5 font-mono text-amber-400 font-bold tracking-widest">
                    {nd.bits}
                  </td>
                  <td className="px-4 py-1.5">
                    <span className="inline-flex items-center gap-1 bg-amber-500/20 border border-amber-500/50 text-amber-300 text-[10px] font-bold px-2 py-0.5 rounded-full uppercase tracking-wider">
                      🎯 MARKER
                    </span>
                  </td>
                </motion.tr>
              ) : (
                <tr
                  key={nd.position}
                  className="border-b border-slate-800/30 hover:bg-slate-800/20"
                >
                  <td className="px-4 py-1.5 font-mono text-slate-500">{nd.position}</td>
                  <td className="px-4 py-1.5 font-mono tracking-widest text-emerald-400">{nd.dna}</td>
                  <td className="px-4 py-1.5 font-mono text-slate-400">{nd.bits}</td>
                  <td className="px-4 py-1.5" />
                </tr>
              )
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
}

// ── Main component ────────────────────────────────────────────────────────────
export default function GroverPOC() {
  const [nCodons, setNCodons]           = useState(2);
  const [isRunning, setIsRunning]       = useState(false);

  // Simulator comparison state
  const [carrierResult, setCarrierResult] = useState(null);
  const [healthyResult, setHealthyResult] = useState(null);
  const [activePatient, setActivePatient] = useState('carrier');
  const [carrierStep, setCarrierStep]     = useState(0);
  const [healthyStep, setHealthyStep]     = useState(0);

  // Live marker info (refreshable)
  const [markerInfo, setMarkerInfo]           = useState(null);
  const [markerLoading, setMarkerLoading]     = useState(false);
  const [markerFetchedAt, setMarkerFetchedAt] = useState(null);
  const [markerFlash, setMarkerFlash]         = useState(false);
  const [markerOffset, setMarkerOffset]       = useState(0);  // index into nearbyWindows (0 = centre)
  const [nearbyWindows, setNearbyWindows]     = useState([]);  // pre-fetched from backend

  // Nonce ref: incremented on every codon change so in-flight fetches for old
  // codon counts cannot overwrite state after the user switches.
  const fetchNonce = useRef(0);

  const fetchMarker = useCallback(async (codons, showToast = false, offset = 0) => {
    const myNonce = fetchNonce.current;
    setMarkerLoading(true);
    try {
      const r = await api.get(`/search/quantum-poc/bio-marker?n_codons=${codons}&force=${showToast}&offset=${offset}`);
      if (fetchNonce.current !== myNonce) return;
      setMarkerInfo(r.data);
      // store the pre-computed windows array if the backend returned it
      if (r.data.nearby_windows?.length) setNearbyWindows(r.data.nearby_windows);
      setMarkerFetchedAt(new Date());
      setMarkerFlash(true);
      setTimeout(() => setMarkerFlash(false), 800);
      if (showToast) toast.success('Marker window shifted — adjacent node around c.5266dupC.');
    } catch (err) {
      if (fetchNonce.current !== myNonce) return;
      // Always surface errors — user explicitly clicked to fetch
      toast.error(err.response?.data?.detail ?? 'Failed to reach NCBI. Try again shortly.');
    } finally {
      if (fetchNonce.current === myNonce) setMarkerLoading(false);
    }
  }, []);

  // Reset offset, nearby windows, custom marker, and stale marker data when codon count changes.
  // Auto-fetch is intentionally removed — the user must click "Fetch from NCBI" explicitly.
  // Incrementing the nonce cancels any in-flight fetch for the previous codon count.
  useEffect(() => {
    fetchNonce.current += 1;
    setMarkerOffset(0);
    setNearbyWindows([]);
    setMarkerInfo(null);
    setMarkerLoading(false);
  }, [nCodons]);

  // Derived: what the marker display should show right now
  // nearbyWindows is sorted offset -5..+5; index 5 = centre (offset 0)
  const CENTRE_IDX = 5;
  const displayWindow = nearbyWindows.length
    ? nearbyWindows[Math.min(Math.max(CENTRE_IDX + markerOffset, 0), nearbyWindows.length - 1)]
    : null;
  const displayDna  = displayWindow?.dna  ?? markerInfo?.marker_dna  ?? null;
  const displayBits = displayWindow?.bits ?? markerInfo?.marker_bits ?? null;

  // Specific windows: exclude windows whose DNA appears in the healthy reference sequence
  // (those produce false positives — Grover finds them in both patients).
  const specificWindows = nearbyWindows.filter(w => !w.in_reference);

  // Auto-correct markerOffset when nearbyWindows loads and the selected window is non-specific.
  useEffect(() => {
    if (!nearbyWindows.length) return;
    const current = nearbyWindows.find(w => w.offset === markerOffset);
    if (current?.in_reference) {
      const specific = nearbyWindows.filter(w => !w.in_reference);
      if (specific.length) {
        setMarkerOffset(
          specific.reduce((best, w) =>
            Math.abs(w.offset) < Math.abs(best.offset) ? w : best
          ).offset
        );
      }
    }
  }, [nearbyWindows]); // eslint-disable-line react-hooks/exhaustive-deps

  const nQubits   = nCodons * 6;
  const maxStates = Math.pow(2, nQubits);

  const handleRun = async () => {
    setIsRunning(true);
    setCarrierResult(null);
    setHealthyResult(null);
    setCarrierStep(0);
    setHealthyStep(0);

    try {
      // Fire both patient scenarios in parallel; share run_id so history groups them
      const ts = new Date().toISOString();
      const runId = crypto.randomUUID();
      const [r1, r2] = await Promise.all([
        api.post('/search/quantum-poc/bio-local', { n_codons: nCodons, has_mutation: true,  marker_offset: markerOffset, client_timestamp: ts, run_id: runId }),
        api.post('/search/quantum-poc/bio-local', { n_codons: nCodons, has_mutation: false, marker_offset: markerOffset, client_timestamp: ts, run_id: runId }),
      ]);
      setCarrierResult(r1.data);
      setHealthyResult(r2.data);
      setActivePatient('carrier');
    } catch (err) {
      toast.error(err.response?.data?.detail ?? 'Quantum execution failed.');
    }
    setIsRunning(false);
  };
  const selectedResult = activePatient === 'carrier' ? carrierResult : healthyResult;
  const selectedStep   = activePatient === 'carrier' ? carrierStep   : healthyStep;
  const setSelectedStep = activePatient === 'carrier' ? setCarrierStep : setHealthyStep;
  const comparisonReady = carrierResult && healthyResult;

  return (
    <div className="w-full max-w-4xl mx-auto space-y-4">

      {/* ── Header ──────────────────────────────────────────────────────────── */}
      <div className="bg-slate-900/60 rounded-3xl px-8 py-6 border border-slate-800 shadow-2xl backdrop-blur-md flex items-center gap-4">
        <div className="w-12 h-12 rounded-2xl bg-red-500/20 flex items-center justify-center ring-1 ring-red-500/15 shrink-0">
          <Dna className="w-6 h-6 text-red-400" />
        </div>
        <div>
          <h2 className="text-2xl font-bold text-white">BRCA1 Constrained Quantum Search</h2>
          <p className="text-slate-400 text-sm mt-0.5">
            Grover's algorithm — search space restricted to real NCBI patient DNA nodes only.
          </p>
        </div>
      </div>

      {/* ── Vertical Pipeline ───────────────────────────────────────────────── */}
      <div className="bg-slate-900/40 rounded-3xl p-6 border border-slate-800 shadow-xl backdrop-blur-md">

        {/* Step 1 — Choose Codon Size */}
        <PipelineStep n={1} label="Choose Codon Size" type="input" done={true}>
          <label className="block text-sm font-medium text-slate-300 mb-3">
            Codons per node
            <span className="ml-2 text-slate-500 font-normal text-xs">(sets qubit count + data volume)</span>
          </label>
          <div className="grid grid-cols-3 gap-2">
            {[1, 2, 3].map((c) => (
              <button key={c} onClick={() => setNCodons(c)}
                className={`py-3 rounded-xl text-sm border transition-all ${
                  nCodons === c
                    ? 'bg-amber-500/20 border-amber-500/60 text-amber-300 font-bold'
                    : 'bg-slate-900 border-slate-800 text-slate-400 hover:border-slate-600'
                }`}>
                <div className="font-mono font-bold text-lg">{c}</div>
                <div className="text-xs opacity-70">{c * 6} qubits</div>
                <div className="text-xs opacity-50">2<sup>{c * 6}</sup> states</div>
              </button>
            ))}
          </div>
          <div className="mt-3 text-xs text-slate-500 bg-slate-950/40 border border-slate-800 rounded-lg px-3 py-2">
            {nCodons} codon{nCodons > 1 ? 's' : ''} = {nCodons * 3} nt/node ·{' '}
            {nQubits} qubits · up to {maxStates.toLocaleString()} states ·{' '}
            <span className="text-amber-400">constrained to actual DNA nodes</span>
          </div>
        </PipelineStep>

        {/* Step 2 — Fetch NCBI Marker */}
        <PipelineStep n={2} label="Fetch NCBI Marker" type="action" done={!!markerInfo}>
          <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4 space-y-3 text-sm">
            <div className="flex items-center justify-between">
              <div className="text-xs font-medium text-slate-400 uppercase tracking-widest">NCBI Data Sources</div>
            </div>
            <div className="flex items-center justify-between gap-3">
              <div className="flex items-center gap-2">
                <span className="text-blue-400 shrink-0">🧬</span>
                <div>
                  <span className="text-slate-300 font-medium">Patient: </span>
                  <span className="font-mono text-blue-300">NM_007294.4</span>
                  <span className="text-slate-500 ml-1 text-xs">(BRCA1 mRNA from NCBI)</span>
                </div>
              </div>
              {!markerInfo ? (
                <button
                  onClick={() => fetchMarker(nCodons, false, 0)}
                  disabled={markerLoading}
                  title="Fetch BRCA1 marker window from NCBI"
                  className="flex items-center gap-1.5 text-xs bg-amber-500/20 border border-amber-500/50 text-amber-300 hover:bg-amber-500/30 px-2.5 py-1.5 rounded-lg transition-all disabled:opacity-40 font-medium shrink-0"
                >
                  <RotateCcw className={`w-3 h-3 ${markerLoading ? 'animate-spin' : ''}`} />
                  {markerLoading ? 'Fetching…' : 'Fetch from NCBI'}
                </button>
              ) : (
                <button
                  onClick={() => {
                    const offsets = specificWindows.map(w => w.offset);
                    if (!offsets.length) return;
                    const idx = offsets.indexOf(markerOffset);
                    const next = offsets[(idx + 1) % offsets.length];
                    setMarkerOffset(next);
                    setMarkerFlash(true);
                    setTimeout(() => setMarkerFlash(false), 800);
                  }}
                  disabled={markerLoading}
                  title="Cycle to next node window around c.5266dupC hotspot"
                  className="flex items-center gap-1 text-xs text-slate-500 hover:text-amber-400 transition-colors disabled:opacity-40 shrink-0"
                >
                  <RotateCcw className={`w-3 h-3 ${markerLoading ? 'animate-spin' : ''}`} />
                  {markerLoading ? 'Fetching…' : `Refresh marker${markerOffset !== 0 ? ` (${markerOffset > 0 ? '+' : ''}${markerOffset})` : ''}`}
                </button>
              )}
            </div>
            <div className="flex items-start gap-2">
              <span className="text-red-400 shrink-0">🎯</span>
              <div className="space-y-0.5">
                <div>
                  <span className="text-slate-300 font-medium">Marker: </span>
                  <span className="font-mono text-red-300">{markerInfo?.marker_variant ?? 'c.5266dupC'}</span>
                  <span className="text-slate-500 ml-1 text-xs">pathogenic — exon 20</span>
                </div>
                <div className={`flex items-center gap-2 flex-wrap transition-all duration-300 ${markerFlash ? 'ring-1 ring-amber-400/60 bg-amber-400/10 rounded px-1' : ''}`}>
                  {markerLoading && !markerInfo ? (
                    <span className="text-xs text-slate-500 italic animate-pulse">fetching marker…</span>
                  ) : !markerInfo ? (
                    <span className="text-xs text-amber-400/60 italic">↑ click "Fetch from NCBI" to load marker</span>
                  ) : (
                    <>
                      <span className="inline-flex items-center px-1.5 py-0.5 rounded-full border text-[9px] font-black uppercase tracking-widest bg-cyan-500/10 border-cyan-500/30 text-cyan-400/80 shrink-0">◀ output</span>
                      <span className="font-mono text-red-400 tracking-widest text-xs bg-red-950/40 px-2 py-0.5 rounded">
                        {displayDna ?? '—'}
                      </span>
                      <span className="font-mono text-slate-500 text-xs">{displayBits ?? '—'}</span>
                      <span className="text-slate-500 text-xs">{markerInfo?.n_qubits ?? nCodons * 6} qubits</span>
                      {markerOffset !== 0 && (
                        <span className="text-xs text-amber-400/60 font-mono">
                          node{markerOffset > 0 ? `+${markerOffset}` : markerOffset}
                        </span>
                      )}
                      {markerLoading && (
                        <span className="text-xs text-amber-400/70 italic animate-pulse">updating…</span>
                      )}
                    </>
                  )}
                </div>
              </div>
            </div>
            <div className="pt-1 border-t border-slate-800 text-xs text-slate-500 space-y-1">
              <div className="flex items-center gap-2">
                <span className="w-2 h-2 rounded-full bg-red-500 shrink-0" />
                Patient A: NM_007294.4 + c.5266dupC insertion (carrier)
              </div>
              <div className="flex items-center gap-2">
                <span className="w-2 h-2 rounded-full bg-emerald-500 shrink-0" />
                Patient B: NM_007294.4 reference only (healthy)
              </div>
            </div>
          </div>
        </PipelineStep>

        {/* Step 3 — Select Marker Window */}
        <PipelineStep n={3} label="Select Marker Window" type="input" optional done={markerOffset !== 0} locked={nearbyWindows.length === 0}>
          <p className="text-xs text-slate-400 mb-2">
            Choose a DNA window around c.5266dupC as the Grover target
            <span className="text-slate-600 ml-1">— optional, defaults to hotspot centre</span>
          </p>
          <select
            value={String(markerOffset)}
            onChange={(e) => setMarkerOffset(parseInt(e.target.value, 10))}
            disabled={specificWindows.length === 0}
            className="w-full bg-slate-950/60 border border-slate-700 rounded-xl px-4 py-2.5 text-sm text-slate-200 outline-none focus:border-amber-500/60 transition-colors disabled:opacity-40 disabled:cursor-not-allowed"
          >
            {specificWindows.map((w) => (
              <option key={w.offset} value={String(w.offset)}>
                {w.offset === 0 ? '★ hotspot (recommended)' : w.offset > 0 ? `node +${w.offset}` : `node ${w.offset}`}
                {' · '}{w.dna}{' · '}{w.bits}
              </option>
            ))}
          </select>
          {nearbyWindows.length === 0 ? (
            <p className="mt-1.5 text-[11px] text-slate-600 italic">Complete step 2 (Fetch from NCBI) to load options.</p>
          ) : specificWindows.length === 0 ? (
            <p className="mt-1.5 text-[11px] text-amber-600/80 italic">⚠ All windows appear in the healthy reference — try increasing codons (step 1).</p>
          ) : markerOffset !== 0 ? (
            <div className="mt-1.5 flex items-center gap-2 text-xs text-emerald-400">
              <span>✓</span>
              <span className="font-mono text-red-400 bg-red-950/40 px-1.5 py-0.5 rounded">
                {nearbyWindows.find((w) => w.offset === markerOffset)?.dna ?? displayDna}
              </span>
              <span className="text-slate-500">selected as Grover target</span>
              <button
                onClick={() => setMarkerOffset(0)}
                className="ml-auto text-slate-500 hover:text-red-400 transition-colors text-[11px]"
              >✕ reset to auto</button>
            </div>
          ) : (
            <p className="mt-1.5 text-[11px] text-slate-500">★ hotspot window selected (default).</p>
          )}
        </PipelineStep>

        {/* Step 4 — Run Quantum Search */}
        <PipelineStep n={4} label="Run Quantum Search" type="action" done={comparisonReady} locked={!markerInfo}>
          <button onClick={handleRun}
            disabled={isRunning}
            className="w-full bg-gradient-to-r from-red-700 to-rose-600 hover:from-red-600 hover:to-rose-500 text-white py-4 text-base rounded-xl font-bold transition-all disabled:opacity-50 shadow-lg shadow-red-900/30 flex items-center justify-center gap-3">
            {isRunning
              ? <><span className="animate-spin inline-block">⚙</span> Running search…</>
              : <><FlaskConical className="w-5 h-5" /> Compare Both Patients
                  <span className="text-sm font-normal opacity-70">({nQubits} qubits each)</span>
                </>
            }
          </button>
        </PipelineStep>

        {/* Step 5 — Results */}
        <AnimatePresence>
          {comparisonReady && (
            <motion.div key="results-step" initial={{ opacity: 0, y: 12 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.4 }}>
              <PipelineStep n={5} label="Quantum Search Results" type="output" done={true} last>

                {/* Shared marker info */}
                <div className="bg-red-950/20 border border-red-800/30 rounded-xl p-4 flex flex-wrap gap-4 items-center text-sm mb-4">
                  <div className="text-xs text-red-400 uppercase tracking-wider shrink-0">Disease Marker</div>
                  <div className="font-mono text-red-300 tracking-widest">{carrierResult.marker_dna}</div>
                  <div className="font-mono text-slate-500 text-xs">{carrierResult.marker_bits}</div>
                  <div className="text-slate-500 text-xs ml-auto">{carrierResult.marker_variant} · {carrierResult.marker_region}</div>
                </div>

                {/* Specificity warning */}
                {carrierResult.marker_in_reference && (
                  <div className="bg-amber-950/40 border border-amber-500/30 rounded-2xl p-4 text-sm space-y-2 mb-4">
                    <div className="font-semibold text-amber-400 flex items-center gap-2">
                      ⚠ Marker Not Specific Enough — {nCodons}-codon ({nCodons * 3}-nt) marker appears in both sequences
                    </div>
                    <p className="text-amber-200/70 text-xs leading-relaxed">
                      With only {nCodons * 3} nucleotides, there are only 4<sup>{nCodons * 3}</sup> = {Math.pow(4, nCodons * 3).toLocaleString()} possible sequences.
                      The BRCA1 reference (~7 kb) has enough nodes that this specific {nCodons * 3}-nt pattern coincidentally appears in the healthy sequence too — so Grover detects it in both patients.
                      This is a real challenge in molecular diagnostics: <strong className="text-amber-300">markers must be long enough to be unique</strong>.
                    </p>
                    <p className="text-amber-300/80 text-xs">
                      {markerOffset !== 0
                        ? '→ Try a different marker window, or increase the codon count.'
                        : <><strong>→ Increase to 2 codons (12 qubits)</strong> or <strong>3 codons (18 qubits)</strong> for a marker specific enough to discriminate.</>
                      }
                    </p>
                  </div>
                )}

                {/* 100% confidence explanation */}
                {(carrierResult.confidence >= 0.995 || healthyResult.confidence >= 0.995) && (
                  <div className="bg-slate-800/40 border border-slate-700/40 rounded-xl px-4 py-3 text-xs text-slate-400 flex items-start gap-2 mb-4">
                    <span className="text-emerald-400 mt-0.5 shrink-0">ℹ</span>
                    <span>
                      <strong className="text-slate-300">Why 100% confidence?</strong>{' '}
                      Grover's algorithm uses exactly <em>k = ⌊π/4·√N⌋</em> iterations, where N is the
                      number of unique DNA nodes. This places the target amplitude at
                      sin²((2k+1)·arcsin(1/√N)) ≈ 1 — so with {carrierResult.n_qubits}-qubit
                      search over ~{carrierResult.total_nodes} nodes, nearly all 1024 shots land on
                      the target state. This is correct and expected behaviour, not an error.
                    </span>
                  </div>
                )}

                {/* Per-patient cards */}
                <SectionLabel type="output" label="Quantum search results" />
                <div className="text-xs text-slate-500 mt-1 mb-4">
                  Same marker searched in two patients — click a card to view step detail
                </div>
                <div className="grid grid-cols-1 sm:grid-cols-2 gap-5">

                  {/* Patient A */}
                  <div className="space-y-3">
                    <MiniDetectionBadge
                      result={carrierResult.detection_result}
                      confidence={carrierResult.confidence}
                      label="Patient A — Carrier (c.5266dupC present)"
                      active={activePatient === 'carrier'}
                      onClick={() => setActivePatient('carrier')}
                    />
                    <div className="grid grid-cols-3 gap-2 text-xs">
                      <div className="bg-slate-950/60 border border-slate-800 rounded-lg p-2">
                        <div className="text-slate-500 uppercase tracking-wider mb-0.5">Nodes</div>
                        <div className="font-mono text-amber-400 font-bold">{carrierResult.total_nodes}</div>
                        <div className="text-slate-500">{carrierResult.n_qubits} qubits</div>
                      </div>
                      <div className="bg-slate-950/60 border border-slate-800 rounded-lg p-2">
                        <div className="text-slate-500 uppercase tracking-wider mb-0.5">Iter</div>
                        <div className="font-mono text-blue-300 font-bold">{carrierResult.iterations}</div>
                        <div className="text-slate-500">Grover steps</div>
                      </div>
                      <div className="bg-slate-950/60 border border-slate-800 rounded-lg p-2">
                        <div className="text-slate-500 uppercase tracking-wider mb-0.5">Marker</div>
                        <div className={`font-bold ${carrierResult.target_found ? 'text-red-400' : 'text-emerald-400'}`}>
                          {carrierResult.target_found ? '⚠ YES' : '✓ NO'}
                        </div>
                        <div className="text-slate-500">{carrierResult.target_found ? 'present' : 'absent'}</div>
                      </div>
                    </div>
                    {carrierResult.nodes_preview && (
                      <NodeTable nodes={carrierResult.nodes_preview} targetBits={carrierResult.marker_bits} />
                    )}
                  </div>

                  {/* Patient B */}
                  <div className="space-y-3">
                    <MiniDetectionBadge
                      result={healthyResult.detection_result}
                      confidence={healthyResult.confidence}
                      label="Patient B — Healthy Control"
                      active={activePatient === 'healthy'}
                      onClick={() => setActivePatient('healthy')}
                    />
                    <div className="grid grid-cols-3 gap-2 text-xs">
                      <div className="bg-slate-950/60 border border-slate-800 rounded-lg p-2">
                        <div className="text-slate-500 uppercase tracking-wider mb-0.5">Nodes</div>
                        <div className="font-mono text-amber-400 font-bold">{healthyResult.total_nodes}</div>
                        <div className="text-slate-500">{healthyResult.n_qubits} qubits</div>
                      </div>
                      <div className="bg-slate-950/60 border border-slate-800 rounded-lg p-2">
                        <div className="text-slate-500 uppercase tracking-wider mb-0.5">Iter</div>
                        <div className="font-mono text-blue-300 font-bold">{healthyResult.iterations}</div>
                        <div className="text-slate-500">Grover steps</div>
                      </div>
                      <div className="bg-slate-950/60 border border-slate-800 rounded-lg p-2">
                        <div className="text-slate-500 uppercase tracking-wider mb-0.5">Marker</div>
                        <div className={`font-bold ${healthyResult.target_found ? 'text-red-400' : 'text-emerald-400'}`}>
                          {healthyResult.target_found ? '⚠ YES' : '✓ NO'}
                        </div>
                        <div className="text-slate-500">{healthyResult.target_found ? 'present' : 'absent'}</div>
                      </div>
                    </div>
                    {healthyResult.nodes_preview && (
                      <NodeTable nodes={healthyResult.nodes_preview} targetBits={healthyResult.marker_bits} />
                    )}
                  </div>
                </div>

                {/* Step navigator */}
                {selectedResult.step_circuits && selectedResult.step_circuits.length > 0 && (
                  <div className="mt-5">
                    <GroverStepNavigator
                      key={activePatient}
                      stepCircuits={selectedResult.step_circuits}
                      activeStep={selectedStep}
                      onStepChange={setSelectedStep}
                      targetBits={selectedResult.marker_bits}
                    />
                  </div>
                )}

              </PipelineStep>
            </motion.div>
          )}
        </AnimatePresence>

      </div>
    </div>
  );
}
