import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Save, Cpu, Dna, Search, FlaskConical, RotateCcw } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import toast from 'react-hot-toast';
import api from '../lib/api';
import GroverStepNavigator from './GroverStepNavigator';
import AmplitudeChart from './AmplitudeChart';

// ── Encoding helpers ──────────────────────────────────────────────────────────
const DEC = { '00': 'A', '01': 'C', 10: 'G', 11: 'T' };
const bitsToDna = (bits) => {
  let s = '';
  for (let i = 0; i < bits.length; i += 2) s += DEC[bits.substring(i, i + 2)] ?? '?';
  return s;
};
const countsToProbabilities = (counts) => {
  if (!counts) return null;
  const total = Object.values(counts).reduce((a, b) => a + b, 0);
  if (!total) return null;
  return Object.fromEntries(Object.entries(counts).map(([k, v]) => [k, v / total]));
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
  const [backendType, setBackendType]   = useState('simulator');
  const [ibmApiKey, setIbmApiKey]       = useState('');
  const [ibmCrn, setIbmCrn]             = useState('');
  const [credsSaved, setCredsSaved]     = useState(false);
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

  // Reset offset and fetch when codon count changes
  useEffect(() => {
    setMarkerOffset(0);
    setNearbyWindows([]);
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

  // IBM state (unchanged)
  const [ibmResult, setIbmResult]       = useState(null);
  const [activeStep, setActiveStep]     = useState(0);
  const [ibmJobId, setIbmJobId]         = useState('');
  const [ibmJobStatus, setIbmJobStatus] = useState('');
  const [ibmBackend, setIbmBackend]     = useState('');
  const ibmPollRef                      = useRef(null);
  const ibmPollCountRef                 = useRef(0);
  const IBM_MAX_POLLS                   = 50; // 50 × 8 s ≈ 7 min

  const stopPoll = useCallback(() => {
    if (ibmPollRef.current) { clearInterval(ibmPollRef.current); ibmPollRef.current = null; }
    ibmPollCountRef.current = 0;
  }, []);
  useEffect(() => () => stopPoll(), [stopPoll]);

  const nQubits   = nCodons * 6;
  const maxStates = Math.pow(2, nQubits);

  useEffect(() => {
    api.get('/credentials').then((r) => {
      if (r.data.ibm_api_key) { setIbmApiKey(r.data.ibm_api_key); setIbmCrn(r.data.ibm_crn); setCredsSaved(true); }
    }).catch(() => {});
  }, []);

  const handleSaveCredentials = async () => {
    try {
      await api.post('/credentials', { api_key: ibmApiKey, crn: ibmCrn });
      setCredsSaved(true);
      toast.success('IBM credentials saved.');
    } catch { toast.error('Failed to save. Are you logged in?'); }
  };

    const pollIbmStatus = useCallback(async (jobId, markerBits) => {
    ibmPollCountRef.current += 1;
    if (ibmPollCountRef.current > IBM_MAX_POLLS) {
      stopPoll();
      setIbmJobStatus('TIMEOUT');
      toast.error('IBM QPU job timed out after 7 minutes. Check IBM Cloud console for job status.');
      return;
    }
    try {
      const res = await api.post('/search/quantum-poc/ibm-status', { job_id: jobId });
      // Guard: component may have unmounted while request was in-flight
      if (!ibmPollRef.current && res.data.status !== 'DONE') return;
      setIbmJobStatus(res.data.status);
      if (res.data.status === 'DONE') {
        stopPoll();
        const total = Object.values(res.data.counts ?? {}).reduce((a, b) => a + b, 0);
        setIbmResult((prev) => ({
          ...prev,
          measured_state:   res.data.measured_state,
          counts:           res.data.counts,
          detection_result: res.data.measured_state === (markerBits ?? prev?.marker_bits) ? 'FOUND' : 'NOT_FOUND',
          confidence:       total > 0 ? (res.data.counts[markerBits ?? prev?.marker_bits] ?? 0) / total : 0,
        }));
        toast.success('IBM QPU job complete!');
      } else if (res.data.status === 'ERROR' || res.data.status === 'CANCELLED') {
        stopPoll();
        toast.error(`IBM job ${res.data.status.toLowerCase()}.`);
      }
    } catch { toast.error('Error checking IBM job status.'); }
  }, [stopPoll]);

  const handleRun = async () => {
    if (backendType === 'ibm_cloud' && !credsSaved) {
      toast.error('Please save your IBM credentials first.'); return;
    }
    setIsRunning(true);
    setCarrierResult(null);
    setHealthyResult(null);
    setIbmResult(null);
    setActiveStep(0);
    setCarrierStep(0);
    setHealthyStep(0);
    setIbmJobId('');
    setIbmJobStatus('');
    setIbmBackend('');
    stopPoll();

    try {
      if (backendType === 'simulator') {
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
      } else {
        const res = await api.post('/search/quantum-poc/bio-ibm-submit', { n_codons: nCodons, client_timestamp: new Date().toISOString() });
        setIbmJobId(res.data.job_id);
        setIbmJobStatus(res.data.status);
        setIbmBackend(res.data.backend ?? '');
        setIbmResult({ ...res.data, measured_state: null, counts: null, detection_result: null, confidence: 0 });
        ibmPollRef.current = setInterval(() => pollIbmStatus(res.data.job_id, res.data.marker_bits), 8000);
      }
    } catch (err) {
      stopPoll();
      toast.error(err.response?.data?.detail ?? 'Quantum execution failed.');
    }
    setIsRunning(false);
  };

  const ibmDone = ibmJobStatus === 'DONE' || ibmJobStatus === 'ERROR' || ibmJobStatus === 'CANCELLED' || ibmJobStatus === 'TIMEOUT';
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
                {backendType === 'simulator' && (
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
                )}
              </div>

              {/* Codon selector — hidden for IBM (locked to 1) */}
              {backendType === 'simulator' && (
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
              )}
            </div>

            {/* RIGHT: Backend + credentials */}
            <div className="space-y-5">
              <div>
                <label className="block text-sm font-medium text-slate-300 mb-2">Execution Backend</label>
                <select value={backendType} onChange={(e) => setBackendType(e.target.value)}
                  className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm outline-none">
                  <option value="simulator">Local Qiskit Simulator (Aer)</option>
                  <option value="ibm_cloud">IBM Cloud QPU (real hardware)</option>
                </select>
                {backendType === 'ibm_cloud' && (
                  <p className="text-xs text-blue-400/80 mt-1.5">
                    IBM free tier: fixed at 6 qubits (1 codon) for QPU compatibility.
                  </p>
                )}
              </div>

              <AnimatePresence>
                {backendType === 'ibm_cloud' && (
                  <motion.div initial={{ opacity: 0, height: 0 }} animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    className="bg-blue-900/10 border border-blue-500/20 p-4 rounded-2xl space-y-3">
                    <div className="text-xs text-blue-300/70">
                      {credsSaved ? '✅ IBM credentials loaded. API key stays server-side.'
                        : '⚠ No IBM credentials saved — enter below.'}
                    </div>
                    <div>
                      <label className="block text-xs text-blue-200/60 mb-1">IAM API Key</label>
                      <input type="password" value={ibmApiKey} onChange={(e) => setIbmApiKey(e.target.value)}
                        className="w-full bg-slate-900 border border-blue-500/30 rounded-lg px-3 py-2 text-sm outline-none" />
                    </div>
                    <div>
                      <label className="block text-xs text-blue-200/60 mb-1">Instance CRN</label>
                      <input type="text" value={ibmCrn} onChange={(e) => setIbmCrn(e.target.value)}
                        className="w-full bg-slate-900 border border-blue-500/30 rounded-lg px-3 py-2 text-sm outline-none" />
                    </div>
                    <button onClick={handleSaveCredentials} disabled={!ibmApiKey || !ibmCrn}
                      className="w-full flex items-center justify-center gap-2 bg-blue-700 hover:bg-blue-600 disabled:opacity-50 text-white text-sm py-2 rounded-lg transition-colors">
                      <Save className="w-4 h-4" /> Save Credentials
                    </button>
                  </motion.div>
                )}
              </AnimatePresence>

              {/* Run button */}
              <button onClick={handleRun}
                disabled={isRunning || (backendType === 'ibm_cloud' && !credsSaved)}
                className="w-full bg-gradient-to-r from-red-700 to-rose-600 hover:from-red-600 hover:to-rose-500 text-white py-4 text-base rounded-xl font-bold transition-all disabled:opacity-50 shadow-lg shadow-red-900/30 flex items-center justify-center gap-3">
                {isRunning
                  ? <><span className="animate-spin inline-block">⚙</span> Running search…</>
                  : backendType === 'simulator'
                    ? <><FlaskConical className="w-5 h-5" /> Compare Both Patients
                        <span className="text-sm font-normal opacity-70">({nQubits} qubits each)</span>
                      </>
                    : <><Search className="w-5 h-5" /> Run BRCA1 Detection
                        <span className="text-sm font-normal opacity-70">(6 qubits)</span>
                      </>
                }
              </button>
            </div>
          </div>

          {/* IBM job status */}
          <AnimatePresence>
            {ibmJobId && (
              <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }}
                className="bg-slate-950/80 border border-blue-900/50 p-4 rounded-2xl">
                <div className="flex items-start justify-between gap-4 flex-wrap">
                  <div className="space-y-1.5">
                    {ibmBackend && (
                      <div className="flex items-center gap-2 text-sm">
                        <Cpu className="w-4 h-4 text-blue-400" />
                        <span className="text-slate-400">Backend:</span>
                        <span className="font-mono text-blue-300">{ibmBackend}</span>
                      </div>
                    )}
                    <div className="text-sm text-slate-400">
                      Job: <span className="font-mono text-slate-300 text-xs break-all ml-1">{ibmJobId}</span>
                    </div>
                    <div className="text-sm flex items-center gap-2">
                      <span className="text-slate-400">Status:</span>
                      <span className={`font-semibold ${
                        ibmJobStatus === 'DONE' ? 'text-emerald-400'
                          : ibmJobStatus === 'ERROR' || ibmJobStatus === 'CANCELLED' ? 'text-red-400'
                          : 'text-amber-400'}`}>{ibmJobStatus}</span>
                      {!ibmDone && (
                        <span className="flex gap-1">
                          {[0,1,2].map((i) => (
                            <motion.span key={i} className="w-1.5 h-1.5 rounded-full bg-amber-400 inline-block"
                              animate={{ opacity: [0.3,1,0.3] }}
                              transition={{ duration: 1.2, repeat: Infinity, delay: i * 0.4 }} />
                          ))}
                        </span>
                      )}
                    </div>
                  </div>
                  {!ibmDone && (
                    <button onClick={() => pollIbmStatus(ibmJobId, ibmResult?.marker_bits)}
                      className="bg-slate-800 hover:bg-slate-700 px-4 py-2 rounded-xl text-sm border border-slate-700 shrink-0">
                      Refresh
                    </button>
                  )}
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* ── SIMULATOR: Comparison Results ─────────────────────────────────── */}
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
                      The BRCA1 reference (~7 kb) has enough nodes that this specific pattern coincidentally appears in the healthy sequence too — so Grover detects it in both patients.
                      This is a real challenge in molecular diagnostics: <strong className="text-amber-300">markers must be long enough to be unique</strong>.
                    </p>
                    <p className="text-amber-300/80 text-xs">
                      → Increase to <strong>2 codons (12 qubits)</strong> or <strong>3 codons (18 qubits)</strong> for a marker specific enough to discriminate.
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

          {/* ── IBM: Single Result ────────────────────────────────────────────── */}
          <AnimatePresence>
            {ibmResult && (
              <motion.div initial={{ opacity: 0, y: 16 }} animate={{ opacity: 1, y: 0 }}
                className="space-y-6 pt-6 border-t border-slate-800">

                {/* Disclaimer: IBM path uses standard (unconstrained) Grover */}
                <div className="bg-blue-950/30 border border-blue-700/30 rounded-xl px-4 py-2.5 text-xs text-blue-300/80 flex items-start gap-2">
                  <span className="text-blue-400 mt-0.5 shrink-0">ℹ</span>
                  <span>
                    <strong className="text-blue-200">Standard Grover on hardware.</strong>{' '}
                    Due to QPU gate-depth limits, the IBM path uses H⊗ⁿ initialisation — uniform over all
                    2<sup>{ibmResult.n_qubits}</sup> states rather than the DNA-constrained |ψ₀⟩ used by the
                    local simulator. The oracle still targets the correct BRCA1 marker.
                  </span>
                </div>

                <div className="grid grid-cols-1 sm:grid-cols-3 gap-4 text-sm">
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4">
                    <div className="text-xs text-slate-500 uppercase tracking-wider mb-1">Patient</div>
                    <div className="font-mono text-blue-300 text-xs">{ibmResult.patient_accession}</div>
                    <div className="text-slate-500 text-xs mt-0.5">{ibmResult.total_nodes} nodes · {ibmResult.n_qubits} qubits</div>
                  </div>
                  <div className="bg-red-950/30 border border-red-800/40 rounded-xl p-4">
                    <div className="text-xs text-red-400 uppercase tracking-wider mb-1">Disease Marker</div>
                    <div className="font-mono text-red-300 tracking-widest text-xs">{ibmResult.marker_dna}</div>
                    <div className="font-mono text-slate-400 text-xs">{ibmResult.marker_bits}</div>
                  </div>
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4">
                    <div className="text-xs text-slate-500 uppercase tracking-wider mb-1">Search Space</div>
                    <div className="font-mono text-amber-400 font-bold">{ibmResult.n_unique}</div>
                    <div className="text-slate-500 text-xs">of {ibmResult.n_unconstrained} states · {ibmResult.iterations} iterations</div>
                  </div>
                </div>

                {ibmResult.counts && (
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-5">
                    <div className="text-xs font-medium text-slate-400 uppercase tracking-widest mb-3">
                      QPU Measurement Distribution
                    </div>
                    <AmplitudeChart
                      probabilities={countsToProbabilities(ibmResult.counts)}
                      targetBits={ibmResult.marker_bits}
                      stepIndex={3}
                    />
                  </div>
                )}

                {ibmResult.circuit_diagram && (
                  <div className="border border-slate-800 rounded-xl overflow-hidden">
                    <div className="px-4 py-2.5 bg-slate-900/80 border-b border-slate-800 text-xs font-medium text-slate-400 uppercase tracking-widest">
                      Transpiled Circuit
                    </div>
                    <div className="p-5 overflow-x-auto">
                      <pre className="text-[10px] text-amber-500/70 leading-tight">{ibmResult.circuit_diagram}</pre>
                    </div>
                  </div>
                )}

                {ibmJobId && !ibmDone && !ibmResult.counts && (
                  <motion.div animate={{ opacity: [0.4,1,0.4] }} transition={{ duration: 1.5, repeat: Infinity }}
                    className="text-center py-4 text-slate-500 text-sm">
                    ⏳ Waiting for IBM QPU — auto-checking every 8s (max ~7 min)
                  </motion.div>
                )}
                {ibmJobStatus === 'TIMEOUT' && (
                  <div className="text-center py-4 text-amber-500/80 text-sm">
                    ⏱ Job timed out. Check the <strong>IBM Cloud console</strong> for final status using job ID: <span className="font-mono text-xs">{ibmJobId}</span>
                  </div>
                )}

                {ibmResult.detection_result && ibmResult.measured_state && (
                  <motion.div
                    initial={{ scale: 0.85, opacity: 0 }}
                    animate={{ scale: 1, opacity: 1 }}
                    className={`rounded-2xl border p-6 text-center ${
                      ibmResult.detection_result === 'FOUND'
                        ? 'bg-red-900/30 border-red-500/50'
                        : 'bg-emerald-900/20 border-emerald-500/40'
                    }`}>
                    <div className={`text-3xl font-bold mb-1 ${ibmResult.detection_result === 'FOUND' ? 'text-red-400' : 'text-emerald-400'}`}>
                      {ibmResult.detection_result === 'FOUND' ? '🧬 BRCA1 DETECTED' : '✅ NOT FOUND'}
                    </div>
                    <div className="mt-2 text-xs text-slate-400">
                      Confidence: <span className="font-mono text-white">{(ibmResult.confidence * 100).toFixed(1)}%</span>
                    </div>
                  </motion.div>
                )}

              </motion.div>
            )}
          </AnimatePresence>

        </div>
      </div>
    </div>
  );
}
