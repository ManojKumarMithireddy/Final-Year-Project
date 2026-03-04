import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Save, Cpu, Dna, Search } from 'lucide-react';
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

// ── Detection badge ───────────────────────────────────────────────────────────
function DetectionBadge({ result, confidence }) {
  const found = result === 'FOUND';
  return (
    <motion.div
      initial={{ scale: 0.85, opacity: 0 }}
      animate={{ scale: 1, opacity: 1 }}
      className={`rounded-2xl border p-6 text-center ${
        found ? 'bg-red-900/30 border-red-500/50' : 'bg-emerald-900/20 border-emerald-500/40'
      }`}>
      <div className={`text-3xl font-bold mb-1 ${found ? 'text-red-400' : 'text-emerald-400'}`}>
        {found ? '🧬 BRCA1 DETECTED' : '✅ NOT FOUND'}
      </div>
      <div className={`text-sm ${found ? 'text-red-300' : 'text-emerald-300'}`}>
        {found ? 'Disease marker c.5266dupC present in patient DNA' : 'No matching disease marker found'}
      </div>
      <div className="mt-2 text-xs text-slate-400">
        Confidence: <span className="font-mono text-white">{(confidence * 100).toFixed(1)}%</span>
        <span className="ml-2 text-slate-600">({Math.round(confidence * 1024)}/1024 shots matched target)</span>
      </div>
    </motion.div>
  );
}

// ── Compact node table ────────────────────────────────────────────────────────
function NodeTable({ nodes, targetBits }) {
  return (
    <div className="border border-slate-800 rounded-xl overflow-hidden">
      <div className="px-4 py-2.5 bg-slate-900/80 border-b border-slate-800 flex items-center justify-between">
        <span className="text-xs font-medium text-slate-400 uppercase tracking-widest">Patient DNA Node Table</span>
        <span className="text-xs text-slate-500">{nodes.length} nodes shown</span>
      </div>
      <div className="max-h-52 overflow-y-auto">
        <table className="w-full text-left text-xs">
          <thead className="sticky top-0 bg-slate-900 border-b border-slate-800">
            <tr>
              <th className="px-4 py-2 text-slate-500 font-medium">#</th>
              <th className="px-4 py-2 text-slate-500 font-medium">DNA</th>
              <th className="px-4 py-2 text-slate-500 font-medium">Bits</th>
              <th className="px-4 py-2 text-slate-500 font-medium">Target</th>
            </tr>
          </thead>
          <tbody>
            {nodes.map((nd) => (
              <tr key={nd.position}
                className={`border-b border-slate-800/30 hover:bg-slate-800/20 ${nd.is_target ? 'bg-amber-900/25' : ''}`}>
                <td className="px-4 py-1.5 font-mono text-slate-500">{nd.position}</td>
                <td className={`px-4 py-1.5 font-mono tracking-widest ${nd.is_target ? 'text-amber-300 font-bold' : 'text-emerald-400'}`}>
                  {nd.dna}
                </td>
                <td className={`px-4 py-1.5 font-mono ${nd.is_target ? 'text-amber-400 font-bold' : 'text-slate-400'}`}>
                  {nd.bits}
                </td>
                <td className="px-4 py-1.5">{nd.is_target && <span className="text-amber-400">🎯</span>}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

// ── Main component ────────────────────────────────────────────────────────────
export default function GroverPOC() {
  const [nCodons, setNCodons]           = useState(1);
  const [backendType, setBackendType]   = useState('simulator');
  const [ibmApiKey, setIbmApiKey]       = useState('');
  const [ibmCrn, setIbmCrn]             = useState('');
  const [credsSaved, setCredsSaved]     = useState(false);
  const [isRunning, setIsRunning]       = useState(false);
  const [result, setResult]             = useState(null);
  const [activeStep, setActiveStep]     = useState(0);
  const [ibmJobId, setIbmJobId]         = useState('');
  const [ibmJobStatus, setIbmJobStatus] = useState('');
  const [ibmBackend, setIbmBackend]     = useState('');
  const ibmPollRef                      = useRef(null);

  const stopPoll = useCallback(() => {
    if (ibmPollRef.current) { clearInterval(ibmPollRef.current); ibmPollRef.current = null; }
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
    try {
      const res = await api.post('/search/quantum-poc/ibm-status', { job_id: jobId });
      setIbmJobStatus(res.data.status);
      if (res.data.status === 'DONE') {
        stopPoll();
        const total = Object.values(res.data.counts ?? {}).reduce((a, b) => a + b, 0);
        setResult((prev) => ({
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
    setResult(null);
    setActiveStep(0);
    setIbmJobId('');
    setIbmJobStatus('');
    setIbmBackend('');
    stopPoll();

    try {
      if (backendType === 'simulator') {
        const res = await api.post('/search/quantum-poc/bio-local', { n_codons: nCodons });
        setResult(res.data);
      } else {
        // IBM path always uses n_codons=1 (6 qubits — free tier)
        const res = await api.post('/search/quantum-poc/bio-ibm-submit', { n_codons: 1 });
        setIbmJobId(res.data.job_id);
        setIbmJobStatus(res.data.status);
        setIbmBackend(res.data.backend ?? '');
        setResult({ ...res.data, measured_state: null, counts: null, detection_result: null, confidence: 0 });
        ibmPollRef.current = setInterval(() => pollIbmStatus(res.data.job_id, res.data.marker_bits), 8000);
      }
    } catch (err) {
      toast.error(err.response?.data?.detail ?? 'Quantum execution failed.');
    }
    setIsRunning(false);
  };

  const ibmDone = ibmJobStatus === 'DONE' || ibmJobStatus === 'ERROR' || ibmJobStatus === 'CANCELLED';

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

        {/* Algorithm note */}
        <div className="bg-amber-950/30 border border-amber-500/20 rounded-2xl p-4 mb-6 text-sm text-amber-200/70 leading-relaxed">
          <span className="font-semibold text-amber-400">Why constrained? </span>
          Standard Grover uses H⊗ⁿ — uniform over <em>all</em> 2ⁿ bit combinations, most of which
          are biologically impossible. This PoC replaces H⊗ⁿ with a state preparation
          <span className="font-mono bg-amber-900/30 px-1 mx-1 rounded">DNA|ψ₀⟩</span>
          that loads only the real DNA nodes from the patient table. The diffusion operator
          <span className="font-mono bg-amber-900/30 px-1 mx-1 rounded">2|ψ₀⟩⟨ψ₀|−I</span>
          keeps all amplitude flow within that valid subset.
        </div>

        <div className="space-y-6">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">

            {/* LEFT: NCBI sources + codon selector */}
            <div className="space-y-5">
              <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4 space-y-3 text-sm">
                <div className="text-xs font-medium text-slate-400 uppercase tracking-widest">NCBI Data Sources</div>
                <div className="flex items-start gap-2">
                  <span className="text-blue-400 shrink-0">🧬</span>
                  <div>
                    <span className="text-slate-300 font-medium">Patient: </span>
                    <span className="font-mono text-blue-300">NM_007294.4</span>
                    <span className="text-slate-500 ml-1 text-xs">(BRCA1 mRNA + c.5266dupC simulated carrier)</span>
                  </div>
                </div>
                <div className="flex items-start gap-2">
                  <span className="text-red-400 shrink-0">🎯</span>
                  <div>
                    <span className="text-slate-300 font-medium">Marker: </span>
                    <span className="font-mono text-red-300">c.5266dupC</span>
                    <span className="text-slate-500 ml-1 text-xs">pathogenic — exon 20</span>
                  </div>
                </div>
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
                  : <><Search className="w-5 h-5" /> Run BRCA1 Detection
                    <span className="text-sm font-normal opacity-70">
                      ({backendType === 'ibm_cloud' ? 6 : nQubits} qubits)
                    </span>
                  </>}
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
                    <button onClick={() => pollIbmStatus(ibmJobId, result?.marker_bits)}
                      className="bg-slate-800 hover:bg-slate-700 px-4 py-2 rounded-xl text-sm border border-slate-700 shrink-0">
                      Refresh
                    </button>
                  )}
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Results */}
          <AnimatePresence>
            {result && (
              <motion.div initial={{ opacity: 0, y: 16 }} animate={{ opacity: 1, y: 0 }}
                className="space-y-6 pt-6 border-t border-slate-800">

                {/* Data provenance row */}
                <div className="grid grid-cols-1 sm:grid-cols-3 gap-4 text-sm">
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4">
                    <div className="text-xs text-slate-500 uppercase tracking-wider mb-1">Patient</div>
                    <div className="font-mono text-blue-300 text-xs">{result.patient_accession}</div>
                    <div className="text-slate-500 text-xs mt-0.5">{result.total_nodes} nodes · {result.n_qubits} qubits</div>
                  </div>
                  <div className="bg-red-950/30 border border-red-800/40 rounded-xl p-4">
                    <div className="text-xs text-red-400 uppercase tracking-wider mb-1">Disease Marker</div>
                    <div className="font-mono text-red-300 tracking-widest text-xs">{result.marker_dna}</div>
                    <div className="font-mono text-slate-400 text-xs">{result.marker_bits}</div>
                  </div>
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4">
                    <div className="text-xs text-slate-500 uppercase tracking-wider mb-1">Search Space</div>
                    <div className="font-mono text-amber-400 font-bold">{result.n_unique}</div>
                    <div className="text-slate-500 text-xs">of {result.n_unconstrained} states · {result.iterations} iterations</div>
                  </div>
                </div>

                {/* Node table */}
                {result.nodes_preview && (
                  <NodeTable nodes={result.nodes_preview} targetBits={result.marker_bits} />
                )}

                {/* Step navigator — simulator path */}
                {result.step_circuits && result.step_circuits.length > 0 && (
                  <GroverStepNavigator
                    stepCircuits={result.step_circuits}
                    activeStep={activeStep}
                    onStepChange={setActiveStep}
                    targetBits={result.marker_bits}
                  />
                )}

                {/* IBM: counts histogram when done */}
                {!result.step_circuits && result.counts && (
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-5">
                    <div className="text-xs font-medium text-slate-400 uppercase tracking-widest mb-3">
                      QPU Measurement Distribution
                      <span className="ml-2 text-slate-600 font-normal normal-case">
                        {Object.values(result.counts).reduce((a, b) => a + b, 0).toLocaleString()} shots
                      </span>
                    </div>
                    <AmplitudeChart
                      probabilities={countsToProbabilities(result.counts)}
                      targetBits={result.marker_bits}
                      stepIndex={3}
                    />
                  </div>
                )}

                {/* IBM circuit diagram */}
                {!result.step_circuits && result.circuit_diagram && (
                  <div className="border border-slate-800 rounded-xl overflow-hidden">
                    <div className="px-4 py-2.5 bg-slate-900/80 border-b border-slate-800 text-xs font-medium text-slate-400 uppercase tracking-widest">
                      Transpiled Circuit
                    </div>
                    <div className="p-5 overflow-x-auto">
                      <pre className="text-[10px] text-amber-500/70 leading-tight">{result.circuit_diagram}</pre>
                    </div>
                  </div>
                )}

                {/* IBM waiting */}
                {ibmJobId && !ibmDone && !result.counts && (
                  <motion.div animate={{ opacity: [0.4,1,0.4] }} transition={{ duration: 1.5, repeat: Infinity }}
                    className="text-center py-4 text-slate-500 text-sm">
                    ⏳ Waiting for IBM QPU — auto-checking every 8s
                  </motion.div>
                )}

                {/* Detection badge */}
                {result.detection_result && result.measured_state && (
                  <DetectionBadge result={result.detection_result} confidence={result.confidence} />
                )}

              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </div>
    </div>
  );
}
