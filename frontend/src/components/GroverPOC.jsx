import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Dna, FlaskConical, RotateCcw } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import toast from 'react-hot-toast';
import api from '../lib/api';
import GroverStepNavigator from './GroverStepNavigator';

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

  const fetchMarker = useCallback(async (codons, showToast = false, offset = 0) => {
    setMarkerLoading(true);
    try {
      const r = await api.get(`/search/quantum-poc/bio-marker?n_codons=${codons}&force=${showToast}&offset=${offset}`);
      setMarkerInfo(r.data);
      // store the pre-computed windows array if the backend returned it
      if (r.data.nearby_windows?.length) setNearbyWindows(r.data.nearby_windows);
      setMarkerFetchedAt(new Date());
      setMarkerFlash(true);
      setTimeout(() => setMarkerFlash(false), 800);
      if (showToast) toast.success('Marker window shifted — adjacent node around c.5266dupC.');
    } catch (err) {
      if (showToast) toast.error(err.response?.data?.detail ?? 'Failed to reach NCBI. Try again shortly.');
    }
    setMarkerLoading(false);
  }, []);

  // Reset offset, nearby windows, and custom marker when codon count changes
  useEffect(() => {
    setMarkerOffset(0);
    setNearbyWindows([]);
    setCustomMarkerDna('');
    fetchMarker(nCodons, false, 0);
  }, [nCodons, fetchMarker]);

  // Derived: what the marker display should show right now
  // nearbyWindows is sorted offset -5..+5; index 5 = centre (offset 0)
  const CENTRE_IDX = 5;
  const displayWindow = nearbyWindows.length
    ? nearbyWindows[Math.min(Math.max(CENTRE_IDX + markerOffset, 0), nearbyWindows.length - 1)]
    : null;
  const displayDna  = displayWindow?.dna  ?? markerInfo?.marker_dna  ?? null;
  const displayBits = displayWindow?.bits ?? markerInfo?.marker_bits ?? null;

  // Custom marker state
  const [useCustomMarker, setUseCustomMarker] = useState(false);
  const [customMarkerDna, setCustomMarkerDna] = useState('');

  const nQubits   = nCodons * 6;
  const maxStates = Math.pow(2, nQubits);

  // Derived: custom marker validation
  const customMarkerLen   = nCodons * 3;
  const customMarkerValid = useCustomMarker
    && customMarkerDna.length === customMarkerLen
    && /^[AGCT]+$/i.test(customMarkerDna);

  const handleRun = async () => {
    if (useCustomMarker && customMarkerDna.length > 0 && !customMarkerValid) {
      toast.error(`Custom marker must be exactly ${customMarkerLen} AGCT characters.`); return;
    }
    setIsRunning(true);
    setCarrierResult(null);
    setHealthyResult(null);
    setCarrierStep(0);
    setHealthyStep(0);

    try {
      // Fire both patient scenarios in parallel; share run_id so history groups them
      const ts = new Date().toISOString();
      const runId = crypto.randomUUID();
      const extra = (useCustomMarker && customMarkerValid)
        ? { custom_marker_dna: customMarkerDna.toUpperCase() }
        : {};
      const [r1, r2] = await Promise.all([
        api.post('/search/quantum-poc/bio-local', { n_codons: nCodons, has_mutation: true,  marker_offset: markerOffset, client_timestamp: ts, run_id: runId, ...extra }),
        api.post('/search/quantum-poc/bio-local', { n_codons: nCodons, has_mutation: false, marker_offset: markerOffset, client_timestamp: ts, run_id: runId, ...extra }),
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
    <div className="max-w-4xl mx-auto">
      <div className="bg-slate-900/60 rounded-3xl p-8 border border-slate-800 shadow-2xl backdrop-blur-md">

        {/* Header */}
        <div className="flex items-center gap-4 mb-6">
          <div className="w-12 h-12 rounded-2xl bg-red-500/20 flex items-center justify-center">
            <Dna className="w-6 h-6 text-red-400" />
          </div>
          <div>
            <h2 className="text-2xl font-bold text-white">BRCA1 Constrained Quantum Search</h2>
            <p className="text-slate-400 text-sm mt-1">
              Grover's algorithm — search space restricted to real NCBI patient DNA nodes only.
            </p>
          </div>
        </div>

        <div className="space-y-6">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">

            {/* LEFT: NCBI sources + codon selector */}
            <div className="space-y-5">
              <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4 space-y-3 text-sm">
                <div className="flex items-center justify-between">
                  <div className="text-xs font-medium text-slate-400 uppercase tracking-widest">NCBI Data Sources</div>
                  <button
                    onClick={() => {
                      const next = markerOffset + 1;
                      const nextOffset = next > 5 ? -5 : next;
                      setMarkerOffset(nextOffset);
                      if (nearbyWindows.length) {
                        // fast path: cycle locally through pre-fetched windows
                        setMarkerFlash(true);
                        setTimeout(() => setMarkerFlash(false), 800);
                      } else {
                        // always works — fetch the specific offset directly from the API
                        fetchMarker(nCodons, true, nextOffset);
                      }
                    }}
                    disabled={markerLoading}
                    title="Cycle to next node window around c.5266dupC hotspot"
                    className="flex items-center gap-1 text-xs text-slate-500 hover:text-amber-400 transition-colors disabled:opacity-40"
                  >
                    <RotateCcw className={`w-3 h-3 ${markerLoading ? 'animate-spin' : ''}`} />
                    {markerLoading ? 'Fetching…' : `Refresh marker${markerOffset !== 0 ? ` (${markerOffset > 0 ? '+' : ''}${markerOffset})` : ''}`}
                  </button>
                </div>
                <div className="flex items-start gap-2">
                  <span className="text-blue-400 shrink-0">🧬</span>
                  <div>
                    <span className="text-slate-300 font-medium">Patient: </span>
                    <span className="font-mono text-blue-300">NM_007294.4</span>
                    <span className="text-slate-500 ml-1 text-xs">(BRCA1 mRNA from NCBI)</span>
                  </div>
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
                      ) : (
                        <>
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

              {/* Codon selector */}
              <div>
                  <label className="block text-sm font-medium text-slate-300 mb-2">
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
                  <div className="mt-2 text-xs text-slate-500 bg-slate-950/40 border border-slate-800 rounded-lg px-3 py-2">
                    {nCodons} codon{nCodons > 1 ? 's' : ''} = {nCodons * 3} nt/node ·{' '}
                    {nQubits} qubits · up to {maxStates.toLocaleString()} states ·{' '}
                    <span className="text-amber-400">constrained to actual DNA nodes</span>
                  </div>
                </div>
            </div>

            {/* RIGHT: Custom marker + run */}
            <div className="space-y-5">

              {/* Custom marker input */}
              <div>
                <label className="flex items-center gap-2 text-sm font-medium text-slate-300 mb-3 cursor-pointer select-none">
                  <input
                    type="checkbox"
                    checked={useCustomMarker}
                    onChange={(e) => { setUseCustomMarker(e.target.checked); if (!e.target.checked) setCustomMarkerDna(''); }}
                    className="w-4 h-4 accent-amber-500"
                  />
                  Custom marker (AGCT)
                  <span className="text-slate-500 font-normal text-xs">— optional, overrides auto-computed</span>
                </label>
                {useCustomMarker && (
                  <div className="space-y-2">
                    <input
                      type="text"
                      value={customMarkerDna}
                      onChange={(e) => setCustomMarkerDna(e.target.value.toUpperCase().replace(/[^AGCT]/g, ''))}
                      maxLength={customMarkerLen}
                      placeholder={`${customMarkerLen} nucleotides (A/G/C/T)`}
                      spellCheck={false}
                      className="w-full bg-slate-950/60 border border-slate-700 rounded-xl px-4 py-2.5 font-mono text-sm text-emerald-300 uppercase outline-none focus:border-amber-500/60 transition-colors"
                    />
                    <div className={`text-xs flex items-center gap-2 ${
                      customMarkerDna.length === 0 ? 'text-slate-500'
                        : customMarkerValid ? 'text-emerald-400' : 'text-red-400'
                    }`}>
                      <span>{customMarkerDna.length}/{customMarkerLen} nt</span>
                      {customMarkerDna.length > 0 && (
                        customMarkerValid
                          ? <span>✓ valid — will use as Grover search target</span>
                          : <span>need exactly {customMarkerLen} characters</span>
                      )}
                    </div>
                  </div>
                )}
              </div>

              {/* Run button */}
              <button onClick={handleRun}
                disabled={isRunning || (useCustomMarker && customMarkerDna.length > 0 && !customMarkerValid)}
                className="w-full bg-gradient-to-r from-red-700 to-rose-600 hover:from-red-600 hover:to-rose-500 text-white py-4 text-base rounded-xl font-bold transition-all disabled:opacity-50 shadow-lg shadow-red-900/30 flex items-center justify-center gap-3">
                {isRunning
                  ? <><span className="animate-spin inline-block">⚙</span> Running search…</>
                  : <><FlaskConical className="w-5 h-5" /> Compare Both Patients
                      <span className="text-sm font-normal opacity-70">({nQubits} qubits each)</span>
                    </>
                }
              </button>
            </div>
          </div>

          {/* ── Comparison Results ─────────────────────────────────── */}
          <AnimatePresence>
            {comparisonReady && (
              <motion.div initial={{ opacity: 0, y: 16 }} animate={{ opacity: 1, y: 0 }}
                className="space-y-6 pt-6 border-t border-slate-800">

                {/* Shared marker info */}
                <div className="bg-red-950/20 border border-red-800/30 rounded-xl p-4 flex flex-wrap gap-4 items-center text-sm">
                  <div className="text-xs text-red-400 uppercase tracking-wider shrink-0">Disease Marker</div>
                  <div className="font-mono text-red-300 tracking-widest">{carrierResult.marker_dna}</div>
                  <div className="font-mono text-slate-500 text-xs">{carrierResult.marker_bits}</div>
                  <div className="text-slate-500 text-xs ml-auto">{carrierResult.marker_variant} · {carrierResult.marker_region}</div>
                </div>

                {/* Specificity warning — marker not unique enough */}
                {carrierResult.marker_in_reference && (
                  <div className="bg-amber-950/40 border border-amber-500/30 rounded-2xl p-4 text-sm space-y-2">
                    <div className="font-semibold text-amber-400 flex items-center gap-2">
                      ⚠ Marker Not Specific Enough — {nCodons}-codon ({nCodons * 3}-nt) marker appears in both sequences
                    </div>
                    <p className="text-amber-200/70 text-xs leading-relaxed">
                      With only {nCodons * 3} nucleotides, there are only 4<sup>{nCodons * 3}</sup> = {Math.pow(4, nCodons * 3).toLocaleString()} possible sequences.
                      The BRCA1 reference (~7 kb) has enough nodes that this specific {nCodons * 3}-nt pattern coincidentally appears in the healthy sequence too — so Grover detects it in both patients.
                      This is a real challenge in molecular diagnostics: <strong className="text-amber-300">markers must be long enough to be unique</strong>.
                    </p>
                    <p className="text-amber-300/80 text-xs">
                      {useCustomMarker
                        ? '→ Try a different custom marker sequence, or increase the codon count.'
                        : <><strong>→ Increase to 2 codons (12 qubits)</strong> or <strong>3 codons (18 qubits)</strong> for a marker specific enough to discriminate.</>
                      }
                    </p>
                  </div>
                )}

                {/* 100% confidence explanation */}
                {(carrierResult.confidence >= 0.995 || healthyResult.confidence >= 0.995) && (
                  <div className="bg-slate-800/40 border border-slate-700/40 rounded-xl px-4 py-3 text-xs text-slate-400 flex items-start gap-2">
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

                {/* ── Per-patient cards with tables ──────────────────────────────── */}
                <div className="text-xs text-slate-500 uppercase tracking-widest">
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
                    {/* Stats */}
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
                    {/* DNA node table */}
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
                    {/* Stats */}
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
                    {/* DNA node table */}
                    {healthyResult.nodes_preview && (
                      <NodeTable nodes={healthyResult.nodes_preview} targetBits={healthyResult.marker_bits} />
                    )}
                  </div>
                </div>

                {/* Step navigator for the selected patient */}
                {selectedResult.step_circuits && selectedResult.step_circuits.length > 0 && (
                  <GroverStepNavigator
                    key={activePatient}
                    stepCircuits={selectedResult.step_circuits}
                    activeStep={selectedStep}
                    onStepChange={setSelectedStep}
                    targetBits={selectedResult.marker_bits}
                  />
                )}

              </motion.div>
            )}
          </AnimatePresence>

        </div>
      </div>
    </div>
  );
}
