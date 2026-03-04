import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Activity, Save, Cpu, Dna, FlaskConical, Search, ChevronDown, ChevronUp } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import toast from 'react-hot-toast';
import api from '../lib/api';
import AmplitudeChart from './AmplitudeChart';

// ── Encoding helpers (2-bit: A=00 C=01 G=10 T=11) ─────────────────────────
const ENC = { A: '00', C: '01', G: '10', T: '11' };
const DEC = { '00': 'A', '01': 'C', 10: 'G', 11: 'T' };
const dnaToBits = (dna) => dna.split('').map((c) => ENC[c] ?? '00').join('');
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

// ── Detection badge ─────────────────────────────────────────────────────────
function DetectionBadge({ result, confidence }) {
  const found = result === 'FOUND';
  return (
    <motion.div
      initial={{ scale: 0.8, opacity: 0 }}
      animate={{ scale: 1, opacity: 1 }}
      className={`flex flex-col items-center justify-center rounded-2xl border p-6 text-center ${
        found
          ? 'bg-red-900/30 border-red-500/50'
          : 'bg-emerald-900/20 border-emerald-500/40'
      }`}>
      <div className={`text-4xl mb-2 ${found ? 'text-red-400' : 'text-emerald-400'}`}>
        {found ? '🧬 DETECTED' : '✅ NOT FOUND'}
      </div>
      <div className={`text-sm font-semibold ${found ? 'text-red-300' : 'text-emerald-300'}`}>
        {found ? 'Disease marker present in search space' : 'No matching disease marker found'}
      </div>
      <div className="mt-2 text-xs text-slate-400">
        Confidence: <span className="font-mono text-white">{(confidence * 100).toFixed(1)}%</span>
      </div>
    </motion.div>
  );
}

// ── Search space comparison bar ─────────────────────────────────────────────
function SearchSpaceBar({ nUnique, nUnconstrained }) {
  const pct = nUnconstrained > 0 ? (nUnique / nUnconstrained) * 100 : 0;
  return (
    <div className="space-y-2">
      <div className="flex justify-between text-xs text-slate-400">
        <span>Constrained search space</span>
        <span className="font-mono text-amber-400">{nUnique} / {nUnconstrained} states</span>
      </div>
      <div className="h-3 bg-slate-800 rounded-full overflow-hidden">
        <motion.div
          initial={{ width: 0 }}
          animate={{ width: `${Math.max(pct, 1)}%` }}
          transition={{ duration: 0.8, ease: 'easeOut' }}
          className="h-full bg-gradient-to-r from-amber-500 to-orange-400 rounded-full"
        />
      </div>
      <div className="text-xs text-slate-500">
        {pct.toFixed(2)}% of the unconstrained 2<sup>{Math.log2(nUnconstrained) | 0}</sup> Hilbert space searched
      </div>
    </div>
  );
}

// ── Node table (scrollable) ─────────────────────────────────────────────────
function NodeTable({ nodes, targetBits, title }) {
  const [expanded, setExpanded] = useState(false);
  const shown = expanded ? nodes : nodes.slice(0, 8);
  return (
    <div className="border border-slate-800 rounded-xl overflow-hidden bg-slate-950/50">
      <div
        className="flex items-center justify-between px-4 py-2.5 bg-slate-900/80 border-b border-slate-800 cursor-pointer select-none"
        onClick={() => setExpanded((p) => !p)}>
        <span className="text-xs font-medium text-slate-400 uppercase tracking-widest">{title}</span>
        <span className="flex items-center gap-1.5 text-xs text-slate-500">
          {nodes.length} nodes {expanded ? <ChevronUp className="w-3.5 h-3.5" /> : <ChevronDown className="w-3.5 h-3.5" />}
        </span>
      </div>
      <div className="max-h-56 overflow-y-auto">
        <table className="w-full text-left text-xs">
          <thead className="sticky top-0 bg-slate-900 border-b border-slate-800">
            <tr>
              <th className="px-4 py-2 text-slate-500 font-medium">#</th>
              <th className="px-4 py-2 text-slate-500 font-medium">DNA</th>
              <th className="px-4 py-2 text-slate-500 font-medium">Bits</th>
              <th className="px-4 py-2 text-slate-500 font-medium">Match</th>
            </tr>
          </thead>
          <tbody>
            {shown.map((nd) => (
              <tr
                key={nd.position}
                className={`border-b border-slate-800/30 hover:bg-slate-800/30 transition-colors ${
                  nd.is_target ? 'bg-amber-900/25' : ''
                }`}>
                <td className="px-4 py-2 font-mono text-slate-500">{nd.position}</td>
                <td className={`px-4 py-2 font-mono tracking-widest ${nd.is_target ? 'text-amber-300 font-bold' : 'text-emerald-400'}`}>
                  {nd.dna}
                </td>
                <td className={`px-4 py-2 font-mono ${nd.is_target ? 'text-amber-400 font-bold' : 'text-slate-400'}`}>
                  {nd.bits}
                </td>
                <td className="px-4 py-2">
                  {nd.is_target && <span className="text-amber-400 text-base">🎯</span>}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {nodes.length > 8 && (
        <button
          onClick={() => setExpanded((p) => !p)}
          className="w-full py-2 text-xs text-slate-500 hover:text-slate-300 transition-colors bg-slate-900/50">
          {expanded ? 'Show less' : `Show all ${nodes.length} nodes`}
        </button>
      )}
    </div>
  );
}

// ── Main component ──────────────────────────────────────────────────────────
export default function GroverPOC() {
  // Controls
  const [nCodons, setNCodons]         = useState(1);
  const [backendType, setBackendType] = useState('simulator');

  // IBM credentials
  const [ibmApiKey, setIbmApiKey]         = useState('');
  const [ibmCrn, setIbmCrn]               = useState('');
  const [credsSaved, setCredsSaved]       = useState(false);

  // Run state
  const [isRunning, setIsRunning]         = useState(false);
  const [result, setResult]               = useState(null);

  // IBM job polling
  const [ibmJobId, setIbmJobId]           = useState('');
  const [ibmJobStatus, setIbmJobStatus]   = useState('');
  const [ibmBackend, setIbmBackend]       = useState('');
  const ibmPollRef                        = useRef(null);

  const stopPoll = useCallback(() => {
    if (ibmPollRef.current) { clearInterval(ibmPollRef.current); ibmPollRef.current = null; }
  }, []);
  useEffect(() => () => stopPoll(), [stopPoll]);

  // Qubit count derived from codons
  const nQubits         = nCodons * 6;
  const maxStates       = Math.pow(2, nQubits);
  const CODON_OPTIONS   = [1, 2, 3];

  // Load saved IBM credentials on mount
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
    } catch { toast.error('Failed to save credentials. Are you logged in?'); }
  };

  // IBM job status polling callback
  const pollIbmStatus = useCallback(async (jobId) => {
    try {
      const res = await api.post('/search/quantum-poc/ibm-status', { job_id: jobId });
      setIbmJobStatus(res.data.status);
      if (res.data.status === 'DONE') {
        stopPoll();
        setResult((prev) => ({
          ...prev,
          measured_state: res.data.measured_state,
          counts: res.data.counts,
          detection_result: res.data.measured_state === prev?.marker_bits ? 'FOUND' : 'NOT_FOUND',
          confidence: res.data.counts
            ? (res.data.counts[prev?.marker_bits] ?? 0) / Object.values(res.data.counts).reduce((a, b) => a + b, 0)
            : 0,
          execution_time_ms: res.data.execution_time_ms ?? 0,
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
      toast.error('Please save your IBM Cloud credentials first.'); return;
    }
    if (backendType === 'ibm_cloud' && nCodons > 2) {
      toast.error('IBM QPU supports up to 2 codons (12 qubits) in this PoC.'); return;
    }
    setIsRunning(true);
    setResult(null);
    setIbmJobId('');
    setIbmJobStatus('');
    setIbmBackend('');
    stopPoll();

    try {
      if (backendType === 'simulator') {
        const res = await api.post('/search/quantum-poc/bio-local', { n_codons: nCodons });
        setResult(res.data);
      } else {
        const res = await api.post('/search/quantum-poc/bio-ibm-submit', { n_codons: nCodons });
        setIbmJobId(res.data.job_id);
        setIbmJobStatus(res.data.status);
        setIbmBackend(res.data.backend ?? '');
        setResult({
          ...res.data,
          measured_state: null,
          counts: null,
          detection_result: null,
          confidence: 0,
        });
        ibmPollRef.current = setInterval(() => pollIbmStatus(res.data.job_id), 8000);
      }
    } catch (err) {
      toast.error(err.response?.data?.detail ?? 'Quantum execution failed.');
    }
    setIsRunning(false);
  };

  const probabilities = result?.counts
    ? countsToProbabilities(result.counts)
    : result?.init_probs
    ? countsToProbabilities(result.init_probs)
    : null;

  return (
    <div className="max-w-4xl mx-auto">
      <div className="bg-slate-900/60 rounded-3xl p-8 border border-slate-800 shadow-2xl backdrop-blur-md">

        {/* ── Header ── */}
        <div className="flex items-center gap-4 mb-6">
          <div className="w-12 h-12 rounded-2xl bg-red-500/20 flex items-center justify-center">
            <Dna className="w-6 h-6 text-red-400" />
          </div>
          <div>
            <h2 className="text-2xl font-bold text-white">BRCA1 Constrained Quantum Search PoC</h2>
            <p className="text-slate-400 text-sm mt-1">
              Grover's algorithm with constrained state preparation — search space restricted to real NCBI patient DNA nodes.
            </p>
          </div>
        </div>

        {/* ── Algorithm callout ── */}
        <div className="bg-amber-950/30 border border-amber-500/20 rounded-2xl p-4 mb-6 text-sm leading-relaxed space-y-1.5">
          <div className="font-semibold text-amber-400 flex items-center gap-2">
            <FlaskConical className="w-4 h-4" /> Why constrained state preparation?
          </div>
          <p className="text-amber-200/70">
            Standard Grover uses H⊗ⁿ to create a <em>uniform superposition over all 2ⁿ bit combinations</em>.
            For genomic data, most of those states correspond to no real DNA sequence — wasting amplitude on
            biologically impossible nodes.
          </p>
          <p className="text-amber-200/70">
            This PoC replaces H⊗ⁿ with <span className="font-mono bg-amber-900/30 px-1 rounded">StatePrep(|ψ₀⟩)</span> that
            loads only the DNA nodes present in the patient data table. The constrained diffusion operator
            is then <span className="font-mono bg-amber-900/30 px-1 rounded">2|ψ₀⟩⟨ψ₀| − I</span>, keeping
            all amplitude flow within the valid search space.
          </p>
        </div>

        <div className="space-y-8">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-8">

            {/* ── LEFT: NCBI data info ── */}
            <div className="space-y-5">
              <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-5 space-y-3">
                <div className="text-xs font-medium text-slate-400 uppercase tracking-widest mb-1">NCBI Data Sources</div>
                <div className="space-y-2 text-sm">
                  <div className="flex items-start gap-2">
                    <span className="text-blue-400 shrink-0 mt-0.5">🧬</span>
                    <div>
                      <span className="text-slate-300 font-medium">Patient sequence: </span>
                      <span className="font-mono text-blue-300">NM_007294.4</span>
                      <span className="text-slate-500 ml-1">(BRCA1 mRNA reference)</span>
                    </div>
                  </div>
                  <div className="flex items-start gap-2">
                    <span className="text-red-400 shrink-0 mt-0.5">🎯</span>
                    <div>
                      <span className="text-slate-300 font-medium">Disease marker: </span>
                      <span className="font-mono text-red-300">c.5266dupC</span>
                      <span className="text-slate-500 ml-1">pathogenic hotspot, exon 20</span>
                    </div>
                  </div>
                </div>
              </div>

              {/* Codon selector */}
              <div>
                <label className="block text-sm font-medium text-slate-300 mb-2">
                  Codons per node
                  <span className="ml-2 text-slate-500 font-normal text-xs">
                    (controls qubits + data volume)
                  </span>
                </label>
                <div className="grid grid-cols-3 gap-2">
                  {CODON_OPTIONS.map((c) => (
                    <button
                      key={c}
                      onClick={() => setNCodons(c)}
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
                  {nCodons} codon{nCodons > 1 ? 's' : ''} = {nCodons * 3} nt/node
                  · {nQubits} qubits · up to {maxStates.toLocaleString()} states · {' '}
                  <span className="text-amber-400">constrained to actual DNA nodes only</span>
                </div>
              </div>

              {/* Backend selector */}
              <div className="pt-4 border-t border-slate-800">
                <label className="block text-sm font-medium text-slate-300 mb-2">Execution Backend</label>
                <select
                  value={backendType}
                  onChange={(e) => setBackendType(e.target.value)}
                  className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm outline-none">
                  <option value="simulator">Local Qiskit Simulator (Aer)</option>
                  <option value="ibm_cloud">IBM Cloud QPU (real hardware)</option>
                </select>
                {backendType === 'ibm_cloud' && nCodons > 2 && (
                  <p className="text-xs text-red-400/80 mt-1.5">
                    ⚠ IBM QPU path is limited to 2 codons (12 qubits). Switch to 1 or 2 codons.
                  </p>
                )}
              </div>

              {/* IBM credentials */}
              <AnimatePresence>
                {backendType === 'ibm_cloud' && (
                  <motion.div
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    className="bg-blue-900/10 border border-blue-500/20 p-4 rounded-2xl space-y-3">
                    <div className="text-xs text-blue-300/70">
                      {credsSaved
                        ? '✅ IBM credentials loaded. API key never leaves the server.'
                        : '⚠ No IBM credentials saved. Enter and save below.'}
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
            </div>

            {/* ── RIGHT: How it works summary ── */}
            <div className="bg-slate-950/40 border border-slate-800 rounded-2xl p-5 text-sm space-y-4 h-fit">
              <div className="text-xs font-medium text-slate-400 uppercase tracking-widest">Pipeline</div>
              {[
                { n: 1, color: 'text-blue-400', bg: 'bg-blue-500/10', label: 'Fetch patient DNA', desc: 'NM_007294.4 from NCBI; split into ' + nCodons * 3 + '-nt nodes' },
                { n: 2, color: 'text-purple-400', bg: 'bg-purple-500/10', label: 'Fetch disease marker', desc: 'c.5266dupC hotspot from NCBI → target bitstring' },
                { n: 3, color: 'text-amber-400', bg: 'bg-amber-500/10', label: 'State preparation', desc: '|ψ₀⟩ = uniform over valid DNA nodes only (not all 2ⁿ)' },
                { n: 4, color: 'text-orange-400', bg: 'bg-orange-500/10', label: 'Constrained Grover', desc: 'Oracle + 2|ψ₀⟩⟨ψ₀|−I diffusion × k iterations' },
                { n: 5, color: 'text-emerald-400', bg: 'bg-emerald-500/10', label: 'Measure & report', desc: 'FOUND if top state = marker bits; confidence from shot distribution' },
              ].map((s) => (
                <div key={s.n} className="flex items-start gap-3">
                  <div className={`w-6 h-6 rounded ${s.bg} flex items-center justify-center text-xs font-bold ${s.color} shrink-0 mt-0.5`}>{s.n}</div>
                  <div>
                    <div className={`font-medium ${s.color}`}>{s.label}</div>
                    <div className="text-slate-500 text-xs mt-0.5">{s.desc}</div>
                  </div>
                </div>
              ))}
            </div>
          </div>

          {/* ── Run button ── */}
          <button
            onClick={handleRun}
            disabled={isRunning || (backendType === 'ibm_cloud' && (!credsSaved || nCodons > 2))}
            className="w-full bg-gradient-to-r from-red-700 to-rose-600 hover:from-red-600 hover:to-rose-500 text-white py-4 text-lg rounded-xl font-bold transition-all disabled:opacity-50 shadow-lg shadow-red-900/30 flex items-center justify-center gap-3">
            {isRunning
              ? <><span className="animate-spin">⚙</span> Running constrained Grover search...</>
              : <><Search className="w-5 h-5" /> Run BRCA1 Detection ({nQubits} qubits)</>}
          </button>

          {/* ── IBM job status ── */}
          <AnimatePresence>
            {ibmJobId && (
              <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }}
                className="bg-slate-950/80 border border-blue-900/50 p-5 rounded-2xl">
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
                      Job ID: <span className="font-mono text-slate-300 text-xs break-all ml-1">{ibmJobId}</span>
                    </div>
                    <div className="text-sm flex items-center gap-2">
                      <span className="text-slate-400">Status:</span>
                      <span className={`font-semibold ${
                        ibmJobStatus === 'DONE' ? 'text-emerald-400'
                          : ibmJobStatus === 'ERROR' || ibmJobStatus === 'CANCELLED' ? 'text-red-400'
                          : 'text-amber-400'
                      }`}>{ibmJobStatus}</span>
                      {ibmJobStatus !== 'DONE' && ibmJobStatus !== 'ERROR' && ibmJobStatus !== 'CANCELLED' && (
                        <span className="flex gap-1">
                          {[0, 1, 2].map((i) => (
                            <motion.span key={i} className="w-1.5 h-1.5 rounded-full bg-amber-400 inline-block"
                              animate={{ opacity: [0.3, 1, 0.3] }}
                              transition={{ duration: 1.2, repeat: Infinity, delay: i * 0.4 }} />
                          ))}
                        </span>
                      )}
                    </div>
                  </div>
                  {ibmJobStatus !== 'DONE' && ibmJobStatus !== 'ERROR' && ibmJobStatus !== 'CANCELLED' && (
                    <button onClick={() => pollIbmStatus(ibmJobId)}
                      className="bg-slate-800 hover:bg-slate-700 px-4 py-2 rounded-xl text-sm border border-slate-700 shrink-0">
                      Refresh now
                    </button>
                  )}
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* ── Results ── */}
          <AnimatePresence>
            {result && (
              <motion.div initial={{ opacity: 0, y: 16 }} animate={{ opacity: 1, y: 0 }}
                className="space-y-6 pt-6 border-t border-slate-800">

                {/* Data provenance */}
                <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-sm">
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-4 space-y-1.5">
                    <div className="text-xs text-slate-500 uppercase tracking-wider font-medium">Patient Data</div>
                    <div className="font-mono text-blue-300">{result.patient_accession}</div>
                    <div className="text-slate-400 text-xs">{result.total_nodes} nodes loaded · {result.n_qubits} qubits</div>
                  </div>
                  <div className="bg-red-950/30 border border-red-800/40 rounded-xl p-4 space-y-1.5">
                    <div className="text-xs text-red-400 uppercase tracking-wider font-medium">Disease Marker (Target)</div>
                    <div className="font-mono text-red-300 tracking-widest">{result.marker_dna}</div>
                    <div className="text-slate-400 text-xs font-mono">{result.marker_bits}</div>
                    <div className="text-slate-500 text-xs">{result.marker_variant} · {result.marker_gene}</div>
                  </div>
                </div>

                {/* Search space comparison */}
                <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-5">
                  <div className="text-xs font-medium text-slate-400 uppercase tracking-widest mb-3">Search Space</div>
                  <SearchSpaceBar nUnique={result.n_unique} nUnconstrained={result.n_unconstrained} />
                  <div className="grid grid-cols-3 gap-3 mt-4 text-center text-sm">
                    <div>
                      <div className="text-slate-500 text-xs uppercase tracking-wider mb-1">Constrained</div>
                      <div className="font-mono text-amber-400 text-lg font-bold">{result.n_unique}</div>
                      <div className="text-slate-600 text-xs">DNA nodes</div>
                    </div>
                    <div className="border-x border-slate-800">
                      <div className="text-slate-500 text-xs uppercase tracking-wider mb-1">Iterations k</div>
                      <div className="font-mono text-white text-lg font-bold">{result.iterations}</div>
                      <div className="text-slate-600 text-xs">⌊π/4·√N⌋</div>
                    </div>
                    <div>
                      <div className="text-slate-500 text-xs uppercase tracking-wider mb-1">Time</div>
                      <div className="font-mono text-white text-lg font-bold">
                        {result.execution_time_ms > 0 ? `${result.execution_time_ms.toFixed(0)}ms` : '—'}
                      </div>
                      <div className="text-slate-600 text-xs">simulation</div>
                    </div>
                  </div>
                </div>

                {/* Patient node table */}
                {result.nodes_preview && (
                  <NodeTable nodes={result.nodes_preview} targetBits={result.marker_bits} title="Patient DNA Node Table" />
                )}

                {/* Probability distribution */}
                {probabilities && (
                  <div className="bg-slate-950/60 border border-slate-800 rounded-xl p-5">
                    <div className="text-xs font-medium text-slate-400 uppercase tracking-widest mb-1">
                      {result.counts ? 'Measurement Distribution (post-Grover)' : 'Initial State Distribution (DNA|ψ₀⟩)'}
                    </div>
                    <div className="text-xs text-slate-600 mb-4">
                      {result.counts
                        ? 'Amplitude-amplified counts after constrained Grover — target state should dominate.'
                        : 'Uniform distribution over valid DNA nodes only (not all 2ⁿ states).'}
                    </div>
                    <AmplitudeChart
                      probabilities={probabilities}
                      targetBits={result.marker_bits}
                      stepIndex={result.counts ? 3 : 0}
                    />
                  </div>
                )}

                {/* Detection result badge — only after measurement */}
                {result.detection_result && result.measured_state && (
                  <DetectionBadge result={result.detection_result} confidence={result.confidence} />
                )}

                {/* IBM waiting state */}
                {!result.detection_result && ibmJobId && (
                  <div className="text-center py-6 text-slate-500 text-sm">
                    <motion.div animate={{ opacity: [0.4, 1, 0.4] }} transition={{ duration: 1.5, repeat: Infinity }}>
                      ⏳ Waiting for IBM QPU job to complete…
                    </motion.div>
                    <div className="text-xs mt-1 text-slate-600">Auto-refreshing every 8s</div>
                  </div>
                )}

                {/* Raw counts table (IBM) */}
                {result.counts && (
                  <div className="border border-slate-800 rounded-xl overflow-hidden">
                    <div className="px-4 py-2.5 bg-slate-900/80 border-b border-slate-800 text-xs font-medium text-slate-400 uppercase tracking-widest">
                      Raw Measurement Counts
                    </div>
                    <div className="max-h-40 overflow-y-auto">
                      <table className="w-full text-xs text-left">
                        <thead className="sticky top-0 bg-slate-900 border-b border-slate-800">
                          <tr>
                            <th className="px-4 py-2 text-slate-500">State</th>
                            <th className="px-4 py-2 text-slate-500">Shots</th>
                            <th className="px-4 py-2 text-slate-500">Probability</th>
                            <th className="px-4 py-2 text-slate-500">DNA</th>
                          </tr>
                        </thead>
                        <tbody>
                          {Object.entries(result.counts)
                            .sort(([, a], [, b]) => b - a)
                            .map(([state, shots]) => {
                              const total = Object.values(result.counts).reduce((a, b) => a + b, 0);
                              const isTop = state === result.measured_state;
                              const isMark = state === result.marker_bits;
                              return (
                                <tr key={state} className={`border-b border-slate-800/30 ${isTop ? 'bg-amber-900/20' : ''}`}>
                                  <td className={`px-4 py-1.5 font-mono tracking-widest ${isTop ? 'text-amber-300 font-bold' : 'text-slate-400'}`}>
                                    |{state}⟩ {isMark && '🎯'}
                                  </td>
                                  <td className="px-4 py-1.5 font-mono text-slate-400">{shots.toLocaleString()}</td>
                                  <td className="px-4 py-1.5 font-mono text-slate-400">{(shots / total * 100).toFixed(1)}%</td>
                                  <td className="px-4 py-1.5 font-mono text-emerald-400">{bitsToDna(state)}</td>
                                </tr>
                              );
                            })}
                        </tbody>
                      </table>
                    </div>
                  </div>
                )}

                {/* Circuit diagram */}
                {result.circuit_diagram && (
                  <div className="border border-slate-800 rounded-xl overflow-hidden">
                    <div className="px-4 py-2.5 bg-slate-900/80 border-b border-slate-800 text-xs font-medium text-slate-400 uppercase tracking-widest">
                      Transpiled Circuit (StatePrep / Oracle / Constrained Diffusion)
                    </div>
                    <div className="p-5 overflow-x-auto">
                      <pre className="text-[10px] text-amber-500/70 leading-tight">{result.circuit_diagram}</pre>
                    </div>
                  </div>
                )}

              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </div>
    </div>
  );
}
