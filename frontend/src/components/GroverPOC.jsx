import React, { useState, useEffect } from 'react';
import { Activity, Save } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import toast from 'react-hot-toast';
import api from '../lib/api';
import AmplitudeChart from './AmplitudeChart';

// ── Crypto helpers ────────────────────────────────────────────────────────────

const dnaToBits = (dna) => {
  const map = { A: '00', C: '01', G: '10', T: '11' };
  return dna.split('').map(b => map[b] || '00').join('');
};
const bitsToDna = (bits) => {
  const map = { '00': 'A', '01': 'C', '10': 'G', '11': 'T' };
  let dna = '';
  for (let i = 0; i < bits.length; i += 2) dna += map[bits.substring(i, i + 2)] || '?';
  return dna;
};
const generateRandomDNA = (len) => {
  const bases = ['A', 'C', 'G', 'T'];
  return Array.from({ length: len }, () => bases[Math.floor(Math.random() * 4)]).join('');
};
const generateOTPKey = (len) => Array.from({ length: len }, () => Math.round(Math.random())).join('');
const xorBits = (a, b) => a.split('').map((bit, i) => bit === b[i] ? '0' : '1').join('');

export default function GroverPOC() {
  const [pocDataset, setPocDataset] = useState(generateRandomDNA(32));
  const [numQubits, setNumQubits] = useState(4);
  const [toyBitstring, setToyBitstring] = useState('0000');

  // Credentials state
  const [ibmApiKey, setIbmApiKey] = useState('');
  const [ibmCrn, setIbmCrn] = useState('');
  const [backendType, setBackendType] = useState('simulator');
  const [saveCredentials, setSaveCredentials] = useState(false);
  const [credsSaved, setCredsSaved] = useState(false);

  // Noise model
  const [noiseLevel, setNoiseLevel] = useState(0.0);

  // Simulation state
  const [toyResult, setToyResult] = useState(null);
  const [isToyRunning, setIsToyRunning] = useState(false);
  const [pocEncryptionState, setPocEncryptionState] = useState(null);
  const [ibmJobId, setIbmJobId] = useState('');
  const [ibmJobStatus, setIbmJobStatus] = useState('');
  const [activeStep, setActiveStep] = useState(0);   // index into step_circuits

  // Derived: chunks for display
  const pocBits = dnaToBits(pocDataset);
  const pocChunks = [];
  for (let i = 0; i < pocBits.length; i += numQubits) {
    if (i + numQubits <= pocBits.length)
      pocChunks.push({ index: pocChunks.length, bits: pocBits.substring(i, i + numQubits) });
  }

  // Auto-fetch saved IBM credentials on mount
  useEffect(() => {
    const fetchCreds = async () => {
      try {
        const res = await api.get('/credentials');
        if (res.data.ibm_api_key) {
          setIbmApiKey(res.data.ibm_api_key);
          setIbmCrn(res.data.ibm_crn);
          setCredsSaved(true);
        }
      } catch (err) { /* not logged in, skip */ }
    };
    fetchCreds();
  }, []);

  const handleSaveCredentials = async () => {
    try {
      await api.post('/credentials', { api_key: ibmApiKey, crn: ibmCrn });
      setCredsSaved(true);
      toast.success('IBM credentials saved successfully.');
    } catch (err) {
      toast.error('Failed to save credentials. Are you logged in?');
    }
  };

  const handleRunToy = async () => {
    if (toyBitstring.length !== numQubits || !/^[01]+$/.test(toyBitstring)) {
      toast.error(`Target bitstring must be exactly ${numQubits} bits long and contain only 0s and 1s.`);
      return;
    }
    if (backendType === 'ibm_cloud' && !credsSaved) {
      toast.error('Please save your IBM Cloud credentials before running on the QPU.');
      return;
    }

    setIsToyRunning(true);
    setToyResult(null);
    setIbmJobStatus('');
    setIbmJobId('');

    // Client-side OTP encryption (PIR demo)
    const otpKey = generateOTPKey(toyBitstring.length);
    const encryptedPayload = xorBits(toyBitstring, otpKey);
    setPocEncryptionState({ bits: toyBitstring, key: otpKey, encrypted: encryptedPayload });

    try {
      if (backendType === 'simulator') {
        const res = await api.post('/search/quantum-simulation-poc',
          { target_bits: encryptedPayload, noise_level: noiseLevel },
        );
        setToyResult(res.data);
        setActiveStep(0);
      } else {
        // IBM: credentials fetched server-side — only target_bits sent
        const res = await api.post('/search/quantum-poc/ibm-submit',
          { target_bits: encryptedPayload },
        );
        setIbmJobId(res.data.job_id);
        setIbmJobStatus(res.data.status);
        setToyResult({
          circuit_diagram: res.data.circuit_diagram,
          iterations: res.data.iterations,
          execution_time_ms: res.data.execution_time_ms,
          measured_state: 'PENDING...'
        });
      }
    } catch (err) {
      const msg = err.response?.data?.detail || 'Error running quantum POC circuit.';
      toast.error(msg);
    }
    setIsToyRunning(false);
  };

  const handleRefreshIbmJob = async () => {
    if (!ibmJobId) return;
    try {
      const res = await api.post('/search/quantum-poc/ibm-status',
        { job_id: ibmJobId },
      );
      setIbmJobStatus(res.data.status);
      if (res.data.status === 'DONE') {
        setToyResult(prev => ({
          ...prev,
          measured_state: res.data.measured_state,
          execution_time_ms: res.data.execution_time_ms
        }));
      }
    } catch (err) {
      toast.error('Error checking IBM job status.');
    }
  };

  return (
    <div className="max-w-4xl mx-auto">
      <div className="bg-slate-900/60 rounded-3xl p-8 border border-slate-800 shadow-2xl backdrop-blur-md">
        {/* Header */}
        <div className="flex items-center gap-4 mb-6">
          <div className="w-12 h-12 rounded-2xl bg-amber-500/20 flex items-center justify-center">
            <Activity className="w-6 h-6 text-amber-400" />
          </div>
          <div>
            <h2 className="text-2xl font-bold text-white">Quantum Sim Proof of Concept</h2>
            <p className="text-slate-400 text-sm mt-1">Blind Grover circuit on local simulator or IBM Cloud QPU.</p>
          </div>
        </div>

        {/* PIR Privacy Callout */}
        <div className="bg-blue-950/40 border border-blue-500/20 rounded-2xl p-4 mb-6 text-sm text-blue-200/80 leading-relaxed">
          <span className="font-semibold text-blue-400">📋 PIR Demo Note: </span>
          Your query is XOR-encrypted with a random OTP key <em>before</em> leaving your browser.
          The server runs Grover's on the encrypted bits — it cannot determine your original query.
          You decrypt the result locally. <span className="text-amber-400/80">⚠ Caveat:</span> A
          cryptographically rigorous PIR requires an oblivious oracle; this is a pedagogical approximation.
        </div>

        <div className="space-y-8">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
            {/* LEFT: Dataset */}
            <div>
              <div className="flex justify-between items-end mb-2">
                <label className="block text-sm font-medium text-slate-300">Random DNA Dataset (Mock)</label>
                <button
                  onClick={() => setPocDataset(generateRandomDNA(32))}
                  className="text-xs bg-slate-800 hover:bg-slate-700 px-3 py-1.5 rounded-lg text-slate-300 transition-colors"
                >
                  Regenerate
                </button>
              </div>
              <div className="bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm text-slate-300 font-mono tracking-widest break-all">
                {pocDataset}
              </div>

              <div className="mt-6">
                <label className="block text-sm font-medium text-slate-300 mb-2">
                  Dataset as {numQubits}-qubit Chunks
                  <span className="text-slate-500 font-normal ml-2">(searching 2<sup>{numQubits}</sup>={2 ** numQubits} states)</span>
                </label>
                <div className="max-h-56 overflow-y-auto border border-slate-800 rounded-xl bg-slate-950/50">
                  <table className="w-full text-left text-sm text-slate-300">
                    <thead className="text-xs text-slate-400 bg-slate-900 sticky top-0 border-b border-slate-800">
                      <tr>
                        <th className="px-5 py-3 font-medium">Chunk</th>
                        <th className="px-5 py-3 font-medium">Bits</th>
                        <th className="px-5 py-3 font-medium">DNA</th>
                      </tr>
                    </thead>
                    <tbody>
                      {pocChunks.map((chunk) => {
                        const isMatch = chunk.bits === toyBitstring;
                        return (
                          <tr key={chunk.index} className={`border-b border-slate-800/30 hover:bg-slate-800/50 transition-colors ${isMatch ? 'bg-amber-900/20' : ''}`}>
                            <td className="px-5 py-2.5 font-mono text-slate-500">{chunk.index}</td>
                            <td className={`px-5 py-2.5 font-mono tracking-widest ${isMatch ? 'text-amber-400 font-bold' : ''}`}>{chunk.bits}</td>
                            <td className="px-5 py-2.5 font-mono text-emerald-400">{bitsToDna(chunk.bits)}</td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
              </div>
            </div>

            {/* RIGHT: Controls */}
            <div className="space-y-5">
              {/* Qubit count */}
              <div>
                <label className="block text-sm font-medium text-slate-300 mb-2">Number of Qubits</label>
                <select
                  value={numQubits}
                  onChange={(e) => { const v = parseInt(e.target.value); setNumQubits(v); setToyBitstring('0'.repeat(v)); }}
                  className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm focus:border-amber-500/50 outline-none transition-all"
                >
                  {[2, 4, 6].map(n => (
                    <option key={n} value={n}>{n} Qubits (searches {2 ** n} states)</option>
                  ))}
                </select>
              </div>

              {/* Target state */}
              <div>
                <label className="block text-sm font-medium text-slate-300 mb-2">Target State (Query)</label>
                <input
                  type="text"
                  value={toyBitstring}
                  onChange={(e) => { if (/^[01]*$/.test(e.target.value)) setToyBitstring(e.target.value); }}
                  maxLength={numQubits}
                  className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm focus:border-amber-500/50 outline-none font-mono transition-all"
                  placeholder={`${numQubits}-bit string`}
                />
                <div className="mt-2 text-xs text-slate-400 flex items-center justify-between">
                  <span>Mapped DNA: <span className="font-mono text-emerald-400 bg-emerald-900/20 px-1 py-0.5 rounded">{bitsToDna(toyBitstring)}</span></span>
                  <span className="opacity-60">(A:00 C:01 G:10 T:11)</span>
                </div>
              </div>

              {/* Backend selector */}
              <div className="pt-4 border-t border-slate-800">
                <label className="block text-sm font-medium text-slate-300 mb-2">Execution Backend</label>
                <select
                  value={backendType}
                  onChange={(e) => setBackendType(e.target.value)}
                  className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm outline-none transition-all"
                >
                  <option value="simulator">Local Qiskit Simulator (Aer)</option>
                  <option value="ibm_cloud">IBM Cloud QPU (real hardware)</option>
                </select>
              </div>

              {/* Noise level — simulator only */}
              <AnimatePresence>
                {backendType === 'simulator' && (
                  <motion.div
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    className="bg-slate-950/40 border border-slate-800 p-4 rounded-2xl"
                  >
                    <label className="block text-sm font-medium text-slate-300 mb-3">
                      Depolarizing Noise Level
                      <span className="ml-2 text-amber-400 font-mono">{(noiseLevel * 100).toFixed(1)}%</span>
                      <span className="ml-2 text-xs text-slate-500">(0% = ideal, 5% = noisy NISQ)</span>
                    </label>
                    <input
                      type="range"
                      min="0" max="0.05" step="0.001"
                      value={noiseLevel}
                      onChange={(e) => setNoiseLevel(parseFloat(e.target.value))}
                      className="w-full accent-amber-500"
                    />
                    {noiseLevel > 0 && (
                      <p className="text-xs text-amber-500/70 mt-2">
                        ⚠ At {(noiseLevel * 100).toFixed(1)}% error rate, the target state may no longer be the most probable measurement outcome — simulating real NISQ hardware degradation.
                      </p>
                    )}
                  </motion.div>
                )}
              </AnimatePresence>

              {/* IBM credentials panel */}
              <AnimatePresence>
                {backendType === 'ibm_cloud' && (
                  <motion.div
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    exit={{ opacity: 0, height: 0 }}
                    className="bg-blue-900/10 border border-blue-500/20 p-5 rounded-2xl"
                  >
                    <div className="text-xs text-blue-300/70 mb-3 font-medium">
                      {credsSaved
                        ? '✅ IBM credentials loaded from your account. The server uses them directly — your API key is never sent in circuit requests.'
                        : '⚠ No IBM credentials saved. Enter and save them below before running.'}
                    </div>
                    <div className="space-y-3">
                      <div>
                        <label className="block text-xs font-medium text-blue-200/60 mb-1">IAM API Key</label>
                        <input type="password" value={ibmApiKey} onChange={(e) => setIbmApiKey(e.target.value)}
                          className="w-full bg-slate-900 border border-blue-500/30 rounded-lg px-3 py-2 text-sm focus:border-blue-400 outline-none" />
                      </div>
                      <div>
                        <label className="block text-xs font-medium text-blue-200/60 mb-1">Instance CRN</label>
                        <input type="text" value={ibmCrn} onChange={(e) => setIbmCrn(e.target.value)}
                          className="w-full bg-slate-900 border border-blue-500/30 rounded-lg px-3 py-2 text-sm focus:border-blue-400 outline-none" />
                      </div>
                      <button
                        onClick={handleSaveCredentials}
                        disabled={!ibmApiKey || !ibmCrn}
                        className="w-full flex items-center justify-center gap-2 bg-blue-700 hover:bg-blue-600 disabled:opacity-50 text-white text-sm py-2 rounded-lg transition-colors"
                      >
                        <Save className="w-4 h-4" /> Save Credentials to Account
                      </button>
                    </div>
                  </motion.div>
                )}
              </AnimatePresence>

              <button
                onClick={handleRunToy}
                disabled={isToyRunning || (backendType === 'ibm_cloud' && !credsSaved)}
                className="w-full bg-gradient-to-r from-amber-600 to-orange-500 hover:from-amber-500 hover:to-orange-400 text-white py-3.5 text-lg rounded-xl font-bold transition-all disabled:opacity-50 mt-2 shadow-lg shadow-amber-500/20"
              >
                {isToyRunning ? 'Executing...' : 'Run Simulation'}
              </button>
            </div>
          </div>

          {/* IBM Job Status */}
          <AnimatePresence>
            {backendType === 'ibm_cloud' && ibmJobId && (
              <motion.div
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                className="bg-slate-950/80 border border-slate-800 p-5 rounded-2xl flex items-center justify-between"
              >
                <div>
                  <div className="text-sm text-slate-400">Job ID: <span className="font-mono text-slate-300 ml-2">{ibmJobId}</span></div>
                  <div className="text-sm mt-1">Status: <span className={`font-semibold ml-2 ${ibmJobStatus === 'DONE' ? 'text-emerald-400' : 'text-amber-400 animate-pulse'}`}>{ibmJobStatus}</span></div>
                </div>
                {ibmJobStatus !== 'DONE' && (
                  <button onClick={handleRefreshIbmJob}
                    className="bg-slate-800 hover:bg-slate-700 px-4 py-2 rounded-xl text-sm transition-colors border border-slate-700">
                    Refresh
                  </button>
                )}
              </motion.div>
            )}
          </AnimatePresence>

          {/* PIR Pipeline Visualization */}
          <AnimatePresence>
            {pocEncryptionState && (
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: 'auto' }}
                className="space-y-6 pt-6 border-t border-slate-800"
              >
                {/* Step 1: Client Encrypt */}
                <div className="bg-slate-900 border border-slate-800 p-6 rounded-2xl shadow-inner">
                  <h3 className="text-sm font-bold text-blue-400 mb-5 flex items-center gap-2">
                    <div className="w-6 h-6 rounded bg-blue-500/20 flex items-center justify-center text-xs">1</div>
                    Client-Side Encryption (OTP-XOR)
                  </h3>
                  <div className="grid grid-cols-1 sm:grid-cols-3 gap-5 text-center font-mono text-sm">
                    <div className="bg-slate-950 p-4 rounded-xl border border-slate-800 shadow-sm">
                      <div className="text-slate-500 text-xs mb-2">Query Bits</div>
                      <div className="text-white text-lg tracking-widest">{pocEncryptionState.bits}</div>
                    </div>
                    <div className="bg-slate-950 p-4 rounded-xl border border-slate-800 text-amber-500 flex flex-col justify-center shadow-sm relative">
                      <div className="absolute top-1/2 -left-6 transform -translate-y-1/2 text-slate-500 text-xl font-sans">⨁</div>
                      <div className="text-amber-600/70 text-xs mb-2 font-sans font-medium uppercase tracking-wider">OTP Key</div>
                      <div className="text-lg tracking-widest">{pocEncryptionState.key}</div>
                      <div className="absolute top-1/2 -right-6 transform -translate-y-1/2 text-slate-500 text-xl font-sans">=</div>
                    </div>
                    <div className="bg-blue-900/20 p-4 rounded-xl border border-blue-500/30 text-blue-400 flex flex-col justify-center shadow-sm">
                      <div className="text-blue-500/70 text-xs mb-2 font-sans font-medium uppercase tracking-wider">Encrypted → Server</div>
                      <div className="text-lg tracking-widest font-bold">{pocEncryptionState.encrypted}</div>
                    </div>
                  </div>
                </div>

                {/* Step 2: Server Execution */}
                {toyResult && (
                  <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }}
                    className="border-l-2 border-dashed border-slate-700 ml-8 pl-8 relative">
                    <div className="absolute -left-[5px] top-1/2 transform -translate-y-1/2 w-2 h-8 bg-slate-700 rounded-full" />
                    <div className="bg-slate-900 border border-slate-800 p-6 rounded-2xl shadow-inner">
                      <h3 className="text-sm font-bold text-amber-400 mb-5 flex items-center gap-2">
                        <div className="w-6 h-6 rounded bg-amber-500/20 flex items-center justify-center text-xs">2</div>
                        Server Execution (Blind — sees encrypted bits only)
                        {toyResult.noise_level > 0 && (
                          <span className="ml-auto text-xs text-amber-500/60 font-normal">
                            noise={(toyResult.noise_level * 100).toFixed(1)}%
                          </span>
                        )}
                      </h3>
                      <div className="grid grid-cols-3 gap-4 text-center text-sm mb-6 bg-slate-950 p-5 rounded-xl border border-slate-800">
                        <div>
                          <span className="block text-slate-500 text-xs uppercase tracking-wider mb-1">State Config</span>
                          <span className="text-white font-mono tracking-widest text-lg">{toyResult.measured_state}</span>
                        </div>
                        <div className="border-x border-slate-800">
                          <span className="block text-slate-500 text-xs uppercase tracking-wider mb-1">Iterations</span>
                          <span className="text-white font-mono text-lg">{toyResult.iterations}</span>
                        </div>
                        <div>
                          <span className="block text-slate-500 text-xs uppercase tracking-wider mb-1">Time (ms)</span>
                          <span className="text-white font-mono text-lg">{toyResult.execution_time_ms.toFixed(1)}</span>
                        </div>
                      </div>
                      {/* ── Step-by-step circuit navigator ── */}
                      {toyResult.step_circuits && toyResult.step_circuits.length > 0 ? (
                        <div className="border border-slate-800 rounded-xl overflow-hidden bg-slate-950">
                          {/* Step tab bar */}
                          <div className="flex border-b border-slate-800 bg-slate-900/60 overflow-x-auto">
                            {toyResult.step_circuits.map((step, idx) => {
                              const colors = [
                                'text-blue-400 border-blue-500 bg-blue-500/10',
                                'text-purple-400 border-purple-500 bg-purple-500/10',
                                'text-amber-400 border-amber-500 bg-amber-500/10',
                                'text-emerald-400 border-emerald-500 bg-emerald-500/10',
                              ];
                              const inactiveColors = [
                                'text-slate-400 hover:text-blue-400',
                                'text-slate-400 hover:text-purple-400',
                                'text-slate-400 hover:text-amber-400',
                                'text-slate-400 hover:text-emerald-400',
                              ];
                              const isActive = activeStep === idx;
                              return (
                                <button
                                  key={idx}
                                  onClick={() => setActiveStep(idx)}
                                  className={`flex-1 min-w-[120px] px-4 py-3 text-xs font-semibold border-b-2 transition-all flex flex-col items-center gap-1 ${isActive
                                    ? colors[idx]
                                    : `border-transparent ${inactiveColors[idx]}`
                                    }`}
                                >
                                  <span className="opacity-60 text-[10px] uppercase tracking-wider">Step {step.step}</span>
                                  <span>{step.label}</span>
                                </button>
                              );
                            })}
                          </div>

                          {/* Active step content */}
                          <AnimatePresence mode="wait">
                            {toyResult.step_circuits.map((step, idx) =>
                              activeStep === idx ? (
                                <motion.div
                                  key={idx}
                                  initial={{ opacity: 0, y: 6 }}
                                  animate={{ opacity: 1, y: 0 }}
                                  exit={{ opacity: 0, y: -6 }}
                                  transition={{ duration: 0.18 }}
                                >
                                  {/* Description banner */}
                                  <div className="px-5 pt-4 pb-3 text-xs text-slate-400 leading-relaxed border-b border-slate-800/60">
                                    {step.description}
                                  </div>

                                  {/* Stack vertically: circuit top, gate cards bottom */}
                                  <div className="flex flex-col">
                                    <div className="p-5 overflow-x-auto border-b border-slate-800/60 pb-8">
                                      <div className="text-[10px] text-slate-500 uppercase tracking-widest mb-4 font-medium flex items-center justify-between">
                                        <span>Composed Circuit Diagram</span>
                                        <span className="text-amber-500/80 normal-case font-normal bg-amber-500/10 px-2 py-0.5 rounded-md border border-amber-500/20">Active step highlighted</span>
                                      </div>
                                      <div className="bg-slate-950/80 p-4 rounded-xl border border-slate-800 shadow-inner overflow-x-auto">
                                        <div className="flex flex-row items-center pointer-events-none w-max">
                                          {toyResult.step_circuits.map((s, idx) => (
                                            <pre
                                              key={`diagram-part-${idx}`}
                                              className={`text-[10px] leading-[11px] m-0 p-0 transition-all duration-300 ${idx === activeStep
                                                ? 'text-amber-400 drop-shadow-[0_0_12px_rgba(251,191,36,0.8)]'
                                                : 'text-slate-300 opacity-30'
                                                }`}
                                            >
                                              {s.diagram}
                                            </pre>
                                          ))}
                                        </div>
                                      </div>
                                    </div>

                                    {/* Gate explanation cards */}
                                    {step.gates && step.gates.length > 0 && (
                                      <div className="p-5 space-y-4 bg-slate-900/40">
                                        <div className="text-[10px] text-slate-500 uppercase tracking-widest font-medium mb-2">
                                          Gates Used in This Step
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                                          {step.gates.map((gate, gi) => {
                                            const gateColors = [
                                              { badge: 'bg-blue-500/20 text-blue-400 border-blue-500/30', card: 'border-blue-500/20 bg-blue-950/20' },
                                              { badge: 'bg-purple-500/20 text-purple-400 border-purple-500/30', card: 'border-purple-500/20 bg-purple-950/20' },
                                              { badge: 'bg-amber-500/20 text-amber-400 border-amber-500/30', card: 'border-amber-500/20 bg-amber-950/20' },
                                              { badge: 'bg-emerald-500/20 text-emerald-400 border-emerald-500/30', card: 'border-emerald-500/20 bg-emerald-950/20' },
                                            ];
                                            const c = gateColors[gi % gateColors.length];
                                            return (
                                              <div key={gi} className={`rounded-xl border p-3 ${c.card}`}>
                                                <div className="flex items-center gap-2 mb-2">
                                                  <span className={`text-xs font-mono font-bold px-2 py-0.5 rounded-md border ${c.badge}`}>
                                                    {gate.symbol}
                                                  </span>
                                                  <div className="flex-1 min-w-0">
                                                    <div className="text-xs font-semibold text-slate-200 truncate">{gate.name}</div>
                                                    <div className="text-[10px] text-slate-500">×{gate.count} application{gate.count !== 1 ? 's' : ''}</div>
                                                  </div>
                                                </div>
                                                <p className="text-[10px] text-slate-400 leading-relaxed">
                                                  {gate.explanation}
                                                </p>
                                              </div>
                                            );
                                          })}
                                        </div>
                                      </div>
                                    )}
                                  </div>

                                  {/* Dynamic amplitude chart */}
                                  {step.probabilities && (
                                    <div className="px-5 py-4 border-t border-slate-800/60">
                                      <AmplitudeChart
                                        probabilities={step.probabilities}
                                        targetBits={pocEncryptionState?.encrypted || toyBitstring}
                                        stepIndex={idx}
                                      />
                                    </div>
                                  )}
                                </motion.div>
                              ) : null
                            )}
                          </AnimatePresence>
                        </div>
                      ) : (
                        /* Fallback: full circuit (IBM path) */
                        toyResult.circuit_diagram && (
                          <div className="border border-slate-800 rounded-xl overflow-hidden bg-slate-950">
                            <div className="p-3 text-slate-500 text-xs flex items-center gap-2 border-b border-slate-800 bg-slate-900/50 uppercase tracking-widest font-medium">
                              Transpiled Circuit
                            </div>
                            <div className="p-5 overflow-x-auto">
                              <pre className="text-[10px] text-amber-500/70 leading-tight">
                                {toyResult.circuit_diagram}
                              </pre>
                            </div>
                          </div>
                        )
                      )}
                    </div>
                  </motion.div>
                )}

                {/* Step 3: Client Decrypt */}
                {toyResult && (
                  <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }}
                    className="bg-slate-900 border border-slate-800 p-6 rounded-2xl shadow-inner">
                    <h3 className="text-sm font-bold text-emerald-400 mb-5 flex items-center gap-2">
                      <div className="w-6 h-6 rounded bg-emerald-500/20 flex items-center justify-center text-xs">3</div>
                      Client-Side Decryption (XOR with OTP key)
                    </h3>
                    <div className="grid grid-cols-1 sm:grid-cols-3 gap-5 text-center font-mono text-sm">
                      <div className="bg-slate-950 p-4 rounded-xl border border-slate-800 shadow-sm">
                        <div className="text-slate-500 text-xs mb-2">Returned State</div>
                        <div className="text-white text-lg tracking-widest">{toyResult.measured_state}</div>
                      </div>
                      <div className="bg-slate-950 p-4 rounded-xl border border-slate-800 text-amber-500 flex flex-col justify-center shadow-sm relative">
                        <div className="absolute top-1/2 -left-6 transform -translate-y-1/2 text-slate-500 text-xl font-sans">⨁</div>
                        <div className="text-amber-600/70 text-xs mb-2 font-sans font-medium uppercase tracking-wider">OTP Key</div>
                        <div className="text-lg tracking-widest">{pocEncryptionState.key}</div>
                        <div className="absolute top-1/2 -right-6 transform -translate-y-1/2 text-slate-500 text-xl font-sans">=</div>
                      </div>
                      <div className="bg-emerald-900/20 p-4 rounded-xl border border-emerald-500/30 text-emerald-400 flex flex-col justify-center shadow-sm relative overflow-hidden">
                        <div className="text-emerald-500/70 text-xs mb-2 font-sans font-medium uppercase tracking-wider">Decrypted</div>
                        <div className="text-lg tracking-widest font-bold">
                          {xorBits(toyResult.measured_state, pocEncryptionState.key)}
                        </div>
                        {xorBits(toyResult.measured_state, pocEncryptionState.key) === pocEncryptionState.bits && (
                          <div className="text-xs text-emerald-500/80 mt-1">✅ Matches original query</div>
                        )}
                        {noiseLevel > 0 && xorBits(toyResult.measured_state, pocEncryptionState.key) !== pocEncryptionState.bits && (
                          <div className="text-xs text-amber-500/80 mt-1">⚠ Noise degraded result</div>
                        )}
                      </div>
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
