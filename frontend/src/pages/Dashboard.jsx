import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { Search, Server, Activity, Dna, Lock, Save } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

const API_BASE = 'http://localhost:8000/api';

// --- Cryptography Helpers ---
const dnaToBits = (dna) => {
  const map = { A: '00', C: '01', G: '10', T: '11' };
  return dna.split('').map(b => map[b] || '00').join('');
};

const bitsToDna = (bits) => {
  const map = { '00': 'A', '01': 'C', '10': 'G', '11': 'T' };
  let dna = '';
  for (let i = 0; i < bits.length; i += 2) {
    const pair = bits.substring(i, i + 2);
    dna += map[pair] || '?';
  }
  return dna;
};

const generateRandomDNA = (len) => {
  const bases = ['A', 'C', 'G', 'T'];
  return Array.from({length: len}, () => bases[Math.floor(Math.random() * bases.length)]).join('');
};

const generateOTPKey = (len) => {
  return Array.from({ length: len }, () => Math.round(Math.random())).join('');
};

const xorBits = (a, b) => {
  return a.split('').map((bit, i) => bit === b[i] ? '0' : '1').join('');
};

export default function Dashboard() {
  const [activeTab, setActiveTab] = useState('search');
  const getAuthHeaders = () => ({ headers: { Authorization: `Bearer ${localStorage.getItem('token')}` } });

  // Tab 1: Hybrid Search State
  const [accession, setAccession] = useState('NM_007294');
  const [targetKmer, setTargetKmer] = useState('ATGCT');
  const [datasetData, setDatasetData] = useState(null);
  
  const [isFetching, setIsFetching] = useState(false);
  const [isSearching, setIsSearching] = useState(false);
  
  const [classicalResult, setClassicalResult] = useState(null);
  const [quantumMetrics, setQuantumMetrics] = useState(null);
  
  // Tab 2: Quantum POC State
  const [pocDataset, setPocDataset] = useState(generateRandomDNA(32));
  const [numQubits, setNumQubits] = useState(4);
  const [toyBitstring, setToyBitstring] = useState('0000');
  
  const pocBits = dnaToBits(pocDataset);
  const pocChunks = [];
  for (let i = 0; i < pocBits.length; i += numQubits) {
    if (i + numQubits <= pocBits.length) {
      pocChunks.push({
        index: pocChunks.length,
        bits: pocBits.substring(i, i + numQubits)
      });
    }
  }
  
  const [toyResult, setToyResult] = useState(null);
  const [isToyRunning, setIsToyRunning] = useState(false);
  const [pocEncryptionState, setPocEncryptionState] = useState(null);

  // IBM Cloud State
  const [backendType, setBackendType] = useState('simulator');
  const [ibmApiKey, setIbmApiKey] = useState('');
  const [ibmCrn, setIbmCrn] = useState('');
  const [saveCredentials, setSaveCredentials] = useState(false);
  const [ibmJobStatus, setIbmJobStatus] = useState('');
  const [ibmJobId, setIbmJobId] = useState('');

  // Auto-fetch credentials
  useEffect(() => {
    const fetchCreds = async () => {
      try {
        const res = await axios.get(`${API_BASE}/credentials`, getAuthHeaders());
        if (res.data.ibm_api_key) {
          setIbmApiKey(res.data.ibm_api_key);
          setIbmCrn(res.data.ibm_crn);
        }
      } catch (err) {
        console.error("Could not load saved credentials");
      }
    };
    fetchCreds();
  }, []);

  const handleFetchDataset = async () => {
    setIsFetching(true);
    setDatasetData(null);
    setClassicalResult(null);
    setQuantumMetrics(null);
    try {
      const res = await axios.get(`${API_BASE}/sequence/${accession}`, getAuthHeaders());
      setDatasetData(res.data);
    } catch (err) {
      alert("Error fetching standard dataset. Check ID or network connection.");
    }
    setIsFetching(false);
  };

  const handleSearch = async () => {
    if (!datasetData) return alert("Fetch a dataset first!");
    
    setIsSearching(true);
    
    try {
      const cRes = await axios.post(`${API_BASE}/search/classical`, {
        dataset: datasetData.sequence,
        target: targetKmer
      }, getAuthHeaders());
      setClassicalResult(cRes.data);
      
      const qRes = await axios.post(`${API_BASE}/search/quantum-simulation`, {
         n_windows: cRes.data.n_windows
      }, getAuthHeaders());
      setQuantumMetrics(qRes.data);
      
    } catch (err) {
      alert("Error running search.");
    }
    setIsSearching(false);
  };

  const handleRunToy = async () => {
    if (toyBitstring.length !== numQubits || !/^[01]+$/.test(toyBitstring)) {
      alert(`Target bitstring must be exactly ${numQubits} bits long and contain only 0s and 1s.`);
      return;
    }

    // Save credentials if checked
    if (backendType === 'ibm_cloud' && saveCredentials) {
      try {
        await axios.post(`${API_BASE}/credentials`, {
          api_key: ibmApiKey,
          crn: ibmCrn
        }, getAuthHeaders());
      } catch (err) {
        console.error("Failed to save credentials", err);
      }
    }

    setIsToyRunning(true);
    setToyResult(null);
    setIbmJobStatus('');
    setIbmJobId('');
    
    const otpKey = generateOTPKey(toyBitstring.length);
    const encryptedPayload = xorBits(toyBitstring, otpKey);
    setPocEncryptionState({ bits: toyBitstring, key: otpKey, encrypted: encryptedPayload });

    try {
      if (backendType === 'simulator') {
        const res = await axios.post(`${API_BASE}/search/quantum-simulation-poc`, {
          target_bits: encryptedPayload
        }, getAuthHeaders());
        setToyResult(res.data);
      } else {
        const res = await axios.post(`${API_BASE}/search/quantum-poc/ibm-submit`, {
          target_bits: encryptedPayload,
          api_key: ibmApiKey,
          crn: ibmCrn
        }, getAuthHeaders());
        setIbmJobId(res.data.job_id);
        setIbmJobStatus(res.data.status);
        setToyResult({
          circuit_diagram: res.data.circuit_diagram,
          iterations: res.data.iterations,
          execution_time_ms: res.data.execution_time_ms,
          measured_state: "PENDING..."
        });
      }
    } catch (err) {
      alert("Error running quantum POC circuit.");
    }
    setIsToyRunning(false);
  };

  const handleRefreshIbmJob = async () => {
    if (!ibmJobId) return;
    try {
      const res = await axios.post(`${API_BASE}/search/quantum-poc/ibm-status`, {
        job_id: ibmJobId,
        api_key: ibmApiKey,
        crn: ibmCrn
      }, getAuthHeaders());
      setIbmJobStatus(res.data.status);
      if (res.data.status === 'DONE') {
        setToyResult(prev => ({
          ...prev,
          measured_state: res.data.measured_state,
          execution_time_ms: res.data.execution_time_ms
        }));
      }
    } catch (err) {
      alert("Error checking IBM job status.");
    }
  };

  const generateChartData = () => {
    if (!quantumMetrics || !classicalResult) return [];
    const data = [];
    const maxN = classicalResult.n_windows;
    const step = Math.max(1, Math.floor(maxN / 10));
    for (let i = step; i <= maxN; i += step) {
      data.push({ windows: i, classical: i, quantum: Math.floor((Math.PI / 4) * Math.sqrt(i)) });
    }
    if (data.length === 0 || data[data.length-1].windows !== maxN) {
        data.push({ windows: maxN, classical: maxN, quantum: Math.floor((Math.PI / 4) * Math.sqrt(maxN)) });
    }
    return data;
  };

  return (
    <div className="max-w-7xl mx-auto p-6 md:p-8">
      {/* Header Tabs with Framer Motion Animation */}
      <div className="flex justify-center mb-8">
        <div className="flex gap-2 bg-slate-900/50 p-1.5 rounded-2xl border border-slate-800 backdrop-blur-md relative">
          {['search', 'poc'].map((tabId) => (
            <button
              key={tabId}
              onClick={() => setActiveTab(tabId)}
              className={`relative px-6 py-2.5 rounded-xl font-medium text-sm transition-colors z-10 ${activeTab === tabId ? 'text-white' : 'text-slate-400 hover:text-white'}`}
            >
              {activeTab === tabId && (
                <motion.div
                  layoutId="activeTabBadge"
                  className="absolute inset-0 bg-blue-600 rounded-xl"
                  initial={false}
                  transition={{ type: "spring", stiffness: 400, damping: 30 }}
                />
              )}
              <span className="relative z-10">{tabId === 'search' ? 'Large-Scale Hybrid Search' : 'Quantum Sim Proof of Concept'}</span>
            </button>
          ))}
        </div>
      </div>

      <AnimatePresence mode="wait">
        <motion.div
          key={activeTab}
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          transition={{ duration: 0.2 }}
        >
          {/* --- TAB 1: LARGE SCALE HYBRID SEARCH --- */}
          {activeTab === 'search' && (
            <div className="grid grid-cols-1 xl:grid-cols-12 gap-8">
              
              {/* LEFT COLUMN: Data Control */}
              <div className="xl:col-span-4 flex flex-col gap-6">
                <div className="bg-slate-900/60 rounded-3xl p-6 border border-slate-800 shadow-2xl backdrop-blur-md">
                  <h2 className="text-xl font-semibold mb-6 flex items-center gap-3">
                    <div className="w-10 h-10 rounded-xl bg-blue-500/20 flex items-center justify-center">
                      <Server className="w-5 h-5 text-blue-400" />
                    </div>
                    Data Configuration
                  </h2>
                  
                  <div className="space-y-5">
                    <div>
                      <label className="block text-sm font-medium text-slate-400 mb-1.5">Standard Dataset ID</label>
                      <div className="flex gap-2">
                        <input 
                          type="text" 
                          value={accession}
                          onChange={(e) => setAccession(e.target.value)}
                          className="w-full bg-slate-950/50 border border-slate-700/50 rounded-xl px-4 py-2.5 text-sm focus:border-blue-500/50 focus:ring-2 focus:ring-blue-500/20 outline-none transition-all"
                          placeholder="e.g. NM_007294"
                        />
                        <button 
                          onClick={handleFetchDataset}
                          disabled={isFetching}
                          className="bg-blue-600 hover:bg-blue-500 px-5 py-2.5 rounded-xl font-medium transition-colors disabled:opacity-50"
                        >
                          {isFetching ? 'Fetching...' : 'Fetch'}
                        </button>
                      </div>
                    </div>

                    <AnimatePresence>
                      {datasetData && (
                        <motion.div 
                          initial={{ opacity: 0, height: 0 }} 
                          animate={{ opacity: 1, height: 'auto' }} 
                          exit={{ opacity: 0, height: 0 }}
                          className="bg-emerald-900/10 p-4 rounded-xl border border-emerald-500/20"
                        >
                          <div className="text-sm text-emerald-400 mb-1.5 font-medium flex items-center gap-2">
                            <span className="w-2 h-2 rounded-full bg-emerald-400 animate-pulse"></span> Data Loaded
                          </div>
                          <div className="text-xs text-slate-400">Total Elements: <span className="text-white font-mono">{datasetData.length.toLocaleString()}</span> units</div>
                          <div className="text-xs text-slate-500 mt-2 truncate font-mono bg-black/20 p-2 rounded">
                            {datasetData.sequence.substring(0, 40)}...
                          </div>
                        </motion.div>
                      )}
                    </AnimatePresence>

                    <div>
                      <label className="block text-sm font-medium text-slate-400 mb-1.5">Target Sequence Pattern</label>
                      <input 
                        type="text" 
                        value={targetKmer}
                        onChange={(e) => setTargetKmer(e.target.value.toUpperCase())}
                        className="w-full bg-slate-950/50 border border-slate-700/50 rounded-xl px-4 py-2.5 text-sm focus:border-purple-500/50 focus:ring-2 focus:ring-purple-500/20 outline-none font-mono transition-all"
                        placeholder="e.g. ATGCT"
                      />
                    </div>
                  </div>
                  
                  <button 
                    onClick={handleSearch}
                    disabled={isSearching || !datasetData}
                    className="w-full mt-6 bg-gradient-to-r from-purple-600 to-blue-600 hover:from-purple-500 hover:to-blue-500 py-3.5 rounded-xl font-bold flex items-center justify-center gap-2 transition-all disabled:opacity-50 disabled:grayscale relative overflow-hidden group"
                  >
                    <div className="absolute inset-0 bg-white/20 translate-y-full group-hover:translate-y-0 transition-transform"></div>
                    <Search className="w-5 h-5 relative z-10" />
                    <span className="relative z-10">{isSearching ? 'Computing...' : 'Run Hybrid Search'}</span>
                  </button>
                </div>
              </div>

              {/* RIGHT COLUMN: Results & Dashboards */}
              <div className="xl:col-span-8 flex flex-col gap-6">
                
                <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                  {/* Classical Reality Panel */}
                  <div className="bg-slate-900/60 rounded-3xl p-6 border border-slate-800 shadow-2xl relative overflow-hidden group">
                    <div className="absolute top-0 right-0 w-32 h-32 bg-blue-500/10 rounded-full blur-[40px] group-hover:bg-blue-500/20 transition-colors" />
                    <h3 className="text-lg font-semibold text-white mb-1">Classical Reality</h3>
                    <p className="text-xs text-slate-400 mb-6">Actual Python sliding window execution</p>
                    
                    {classicalResult ? (
                      <div className="space-y-4 relative z-10">
                        <div>
                          <div className="text-sm text-slate-400 mb-1">Search Space (Windows)</div>
                          <div className="text-3xl font-mono text-white tracking-tight">{classicalResult.n_windows.toLocaleString()}</div>
                        </div>
                        <div>
                          <div className="text-sm text-slate-400 mb-1">Required Queries O(N)</div>
                          <div className="text-2xl font-mono text-blue-400 tracking-tight">{classicalResult.classical_queries.toLocaleString()}</div>
                        </div>
                        <div className="grid grid-cols-2 gap-4 mt-4">
                           <div className="bg-slate-950/50 border border-slate-800/50 p-4 rounded-2xl">
                              <div className="text-xs text-slate-500 mb-1">Matches Found</div>
                              <div className="text-xl font-semibold text-emerald-400">{classicalResult.matches.length}</div>
                           </div>
                           <div className="bg-slate-950/50 border border-slate-800/50 p-4 rounded-2xl">
                              <div className="text-xs text-slate-500 mb-1">Exec Time</div>
                              <div className="text-xl font-semibold text-slate-300">{classicalResult.execution_time_ms.toFixed(1)} <span className="text-sm text-slate-500">ms</span></div>
                           </div>
                        </div>
                      </div>
                    ) : (
                      <div className="h-48 flex items-center justify-center text-slate-600/70 text-sm">
                        Run a search to view classical metrics
                      </div>
                    )}
                  </div>

                  {/* Quantum Theory Panel */}
                  <div className="bg-slate-900/60 rounded-3xl p-6 border border-slate-800 shadow-2xl relative overflow-hidden group">
                    <div className="absolute -bottom-10 -right-10 w-48 h-48 bg-purple-500/10 rounded-full blur-[50px] group-hover:bg-purple-500/20 transition-colors" />
                    <Lock className="absolute -bottom-4 -right-4 w-32 h-32 text-slate-800/30 rotate-12" />
                    <h3 className="text-lg font-semibold text-white mb-1">Quantum Theory</h3>
                    <p className="text-xs text-slate-400 mb-6">Theoretical Grover's requirements</p>

                    {quantumMetrics ? (
                      <div className="space-y-5 relative z-10">
                        <div>
                          <div className="text-sm text-slate-400 mb-1">Required Oracle Calls O(√N)</div>
                          <div className="text-4xl font-mono text-purple-400 tracking-tight">{quantumMetrics.quantum_queries.toLocaleString()}</div>
                        </div>
                        <div className="bg-gradient-to-br from-purple-900/40 to-blue-900/20 border border-purple-500/30 p-5 rounded-2xl backdrop-blur-sm">
                          <div className="text-sm text-purple-300 mb-2 font-medium">Scaling Advantage</div>
                          <div className="text-2xl font-bold text-white mb-2 flex items-baseline gap-2">
                            {quantumMetrics.advantage_ratio.toLocaleString()}x
                            <span className="text-sm font-normal text-purple-200/60">fewer queries</span>
                          </div>
                          <div className="text-xs text-slate-400">
                            (Classical requires {classicalResult.classical_queries.toLocaleString()} checks)
                          </div>
                        </div>
                      </div>
                    ) : (
                      <div className="h-48 flex items-center justify-center text-slate-600/70 text-sm relative z-10">
                        Run a search to view theoretical advantage
                      </div>
                    )}
                  </div>
                </div>

                {/* Chart Panel */}
                <div className="bg-slate-900/60 rounded-3xl p-6 border border-slate-800 shadow-2xl flex-1 flex flex-col min-h-[400px]">
                  <h3 className="text-lg font-semibold text-white mb-6">Complexity Divergence</h3>
                  
                  {quantumMetrics ? (
                    <motion.div initial={{opacity: 0}} animate={{opacity: 1}} className="flex-1 w-full relative">
                      <ResponsiveContainer width="100%" height="100%">
                        <LineChart data={generateChartData()} margin={{ top: 15, right: 30, left: 10, bottom: 20 }}>
                          <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" vertical={false} />
                          <XAxis 
                            dataKey="windows" 
                            stroke="#64748b" 
                            tickFormatter={(val) => `${(val/1000).toFixed(0)}k`}
                            tick={{fontSize: 12}}
                            axisLine={false}
                            tickLine={false}
                            dy={10}
                          />
                          <YAxis 
                            stroke="#64748b"
                            tickFormatter={(val) => `${(val/1000).toFixed(0)}k`}
                            tick={{fontSize: 12}}
                            axisLine={false}
                            tickLine={false}
                            dx={-10}
                          />
                          <Tooltip 
                            contentStyle={{ backgroundColor: '#0f172a', borderColor: '#334155', color: '#f8fafc', borderRadius: '12px', boxShadow: '0 10px 15px -3px rgb(0 0 0 / 0.5)' }}
                            itemStyle={{ color: '#e2e8f0' }}
                            cursor={{stroke: '#334155', strokeWidth: 1, strokeDasharray: '5 5'}}
                          />
                          <Legend verticalAlign="top" height={40} iconType="circle"/>
                          <Line type="monotone" dataKey="classical" name="Classical O(N)" stroke="#3b82f6" strokeWidth={3} dot={false} activeDot={{r: 6, fill: '#3b82f6', stroke: '#0f172a', strokeWidth: 2}} />
                          <Line type="monotone" dataKey="quantum" name="Quantum O(√N)" stroke="#a855f7" strokeWidth={3} dot={false} activeDot={{r: 6, fill: '#a855f7', stroke: '#0f172a', strokeWidth: 2}} />
                        </LineChart>
                      </ResponsiveContainer>
                    </motion.div>
                  ) : (
                    <div className="flex-1 flex items-center justify-center text-slate-600/70 text-sm">
                      Chart will appear after running a search
                    </div>
                  )}
                </div>
              </div>
            </div>
          )}

          {/* --- TAB 2: QUANTUM SIMULATION POC --- */}
          {activeTab === 'poc' && (
            <div className="max-w-4xl mx-auto">
              <div className="bg-slate-900/60 rounded-3xl p-8 border border-slate-800 shadow-2xl backdrop-blur-md">
                <div className="flex items-center gap-4 mb-6">
                  <div className="w-12 h-12 rounded-2xl bg-amber-500/20 flex items-center justify-center">
                    <Activity className="w-6 h-6 text-amber-400" />
                  </div>
                  <div>
                    <h2 className="text-2xl font-bold text-white">Quantum Sim Proof of Concept</h2>
                    <p className="text-slate-400 text-sm mt-1">Execute a blind Grover circuit on local simulator or IBM Cloud.</p>
                  </div>
                </div>
                
                <div className="space-y-8">
                  {/* Setup Section */}
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                    <div>
                      <div className="flex justify-between items-end mb-2">
                        <label className="block text-sm font-medium text-slate-300">Random DNA Dataset (Mock)</label>
                        <button onClick={() => setPocDataset(generateRandomDNA(32))} className="text-xs bg-slate-800 hover:bg-slate-700 px-3 py-1.5 rounded-lg text-slate-300 transition-colors">Regenerate</button>
                      </div>
                      <div className="bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm text-slate-300 font-mono tracking-widest break-all">
                        {pocDataset}
                      </div>

                      <div className="mt-6">
                        <label className="block text-sm font-medium text-slate-300 mb-2">Dataset Broken into {numQubits}-Qubit Chunks</label>
                        <div className="max-h-56 overflow-y-auto border border-slate-800 rounded-xl bg-slate-950/50 custom-scrollbar">
                          <table className="w-full text-left text-sm text-slate-300">
                             <thead className="text-xs text-slate-400 bg-slate-900 sticky top-0 border-b border-slate-800 z-10">
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
                                   )
                                })}
                             </tbody>
                          </table>
                        </div>
                      </div>
                    </div>

                    <div className="space-y-6">
                      <div>
                        <label className="block text-sm font-medium text-slate-300 mb-2">Number of Qubits</label>
                        <select 
                          value={numQubits}
                          onChange={(e) => {
                            const val = parseInt(e.target.value);
                            setNumQubits(val);
                            setToyBitstring('0'.repeat(val));
                          }}
                          className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm focus:border-amber-500/50 focus:ring-2 focus:ring-amber-500/20 outline-none transition-all"
                        >
                          {[2, 3, 4, 5, 6].map(n => (
                            <option key={n} value={n}>{n} Qubits</option>
                          ))}
                        </select>
                      </div>

                      <div>
                        <label className="block text-sm font-medium text-slate-300 mb-2">Target State (Query)</label>
                        <input 
                          type="text"
                          value={toyBitstring}
                          onChange={(e) => {
                            const val = e.target.value;
                            if (/^[01]*$/.test(val)) setToyBitstring(val);
                          }}
                          maxLength={numQubits}
                          className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm focus:border-amber-500/50 focus:ring-2 focus:ring-amber-500/20 outline-none font-mono transition-all"
                          placeholder={`${numQubits}-bit string`}
                        />
                        <div className="mt-2 text-xs text-slate-400 flex items-center justify-between">
                           <span>Mapped DNA: <span className="font-mono text-emerald-400 bg-emerald-900/20 px-1 py-0.5 rounded">{bitsToDna(toyBitstring)}</span></span>
                           <span className="opacity-60">(A:00 C:01 G:10 T:11)</span>
                        </div>
                      </div>

                      <div className="pt-6 border-t border-slate-800">
                        <label className="block text-sm font-medium text-slate-300 mb-2">Execution Backend</label>
                        <select 
                          value={backendType}
                          onChange={(e) => setBackendType(e.target.value)}
                          className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-sm focus:border-blue-500/50 outline-none transition-all"
                        >
                          <option value="simulator">Local Qiskit Simulator</option>
                          <option value="ibm_cloud">IBM Cloud QPU</option>
                        </select>
                      </div>

                      {/* IBM Credentials Form */}
                      <AnimatePresence>
                        {backendType === 'ibm_cloud' && (
                          <motion.div 
                            initial={{ opacity: 0, height: 0 }}
                            animate={{ opacity: 1, height: 'auto' }}
                            exit={{ opacity: 0, height: 0 }}
                            className="bg-blue-900/10 border border-blue-500/20 p-5 rounded-2xl relative overflow-hidden"
                          >
                            <div className="space-y-4 relative z-10">
                              <div>
                                <label className="block text-xs font-medium text-blue-200/60 mb-1">IAM API Key</label>
                                <input 
                                  type="password" 
                                  value={ibmApiKey}
                                  onChange={(e) => setIbmApiKey(e.target.value)}
                                  className="w-full bg-slate-900 border border-blue-500/30 rounded-lg px-3 py-2 text-sm focus:border-blue-400 outline-none"
                                />
                              </div>
                              <div>
                                <label className="block text-xs font-medium text-blue-200/60 mb-1">Instance CRN</label>
                                <input 
                                  type="text" 
                                  value={ibmCrn}
                                  onChange={(e) => setIbmCrn(e.target.value)}
                                  className="w-full bg-slate-900 border border-blue-500/30 rounded-lg px-3 py-2 text-sm focus:border-blue-400 outline-none"
                                />
                              </div>
                              <label className="flex items-center gap-2 cursor-pointer group mt-2">
                                <input 
                                  type="checkbox" 
                                  checked={saveCredentials}
                                  onChange={(e) => setSaveCredentials(e.target.checked)}
                                  className="w-4 h-4 rounded border-slate-600 text-blue-500 focus:ring-blue-500 bg-slate-800"
                                />
                                <span className="text-xs text-blue-200 group-hover:text-blue-100 transition-colors flex items-center gap-1">
                                  <Save className="w-3 h-3" /> Save to database on submit
                                </span>
                              </label>
                            </div>
                          </motion.div>
                        )}
                      </AnimatePresence>

                      <button 
                        onClick={handleRunToy}
                        disabled={isToyRunning || (backendType === 'ibm_cloud' && (!ibmApiKey || !ibmCrn))}
                        className="w-full bg-gradient-to-r from-amber-600 to-orange-500 hover:from-amber-500 hover:to-orange-400 text-white py-3.5 text-lg rounded-xl font-bold transition-all disabled:opacity-50 mt-4 shadow-lg shadow-amber-500/20"
                      >
                        {isToyRunning ? 'Executing...' : `Run Simulation`}
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
                          <button 
                            onClick={handleRefreshIbmJob}
                            className="bg-slate-800 hover:bg-slate-700 px-4 py-2 rounded-xl text-sm transition-colors border border-slate-700"
                          >
                            Refresh
                          </button>
                        )}
                      </motion.div>
                    )}
                  </AnimatePresence>

                  {/* Pipeline Visualization */}
                  <AnimatePresence>
                    {pocEncryptionState && (
                      <motion.div 
                        initial={{ opacity: 0, height: 0 }}
                        animate={{ opacity: 1, height: 'auto' }}
                        className="space-y-6 pt-6 border-t border-slate-800"
                      >
                        {/* 1. Client Side */}
                        <div className="bg-slate-900 border border-slate-800 p-6 rounded-2xl relative shadow-inner">
                          <h3 className="text-sm font-bold text-blue-400 mb-5 flex items-center gap-2">
                             <div className="w-6 h-6 rounded bg-blue-500/20 flex items-center justify-center text-xs">1</div> 
                             Client-Side Encryption
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
                              <div className="text-blue-500/70 text-xs mb-2 font-sans font-medium uppercase tracking-wider">Encrypted</div>
                              <div className="text-lg tracking-widest font-bold">{pocEncryptionState.encrypted}</div>
                            </div>
                          </div>
                        </div>

                        {/* 2. Server Side */}
                        {toyResult && (
                          <motion.div initial={{opacity:0}} animate={{opacity:1}} className="border-l-2 border-dashed border-slate-700 ml-8 pl-8 relative">
                            <div className="absolute -left-[5px] top-1/2 transform -translate-y-1/2 w-2 h-8 bg-slate-700 rounded-full"></div>
                            
                            <div className="bg-slate-900 border border-slate-800 p-6 rounded-2xl shadow-inner">
                              <h3 className="text-sm font-bold text-amber-400 mb-5 flex items-center gap-2">
                                <div className="w-6 h-6 rounded bg-amber-500/20 flex items-center justify-center text-xs">2</div> 
                                Server Execution (Blind)
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
                              
                              {toyResult.circuit_diagram && (
                                <div className="border border-slate-800 rounded-xl overflow-hidden bg-slate-950">
                                  <div className="p-3 text-slate-500 text-xs flex items-center gap-2 border-b border-slate-800 bg-slate-900/50 uppercase tracking-widest font-medium">
                                    Transpiled Circuit
                                  </div>
                                  <div className="p-5 overflow-x-auto custom-scrollbar">
                                    <pre className="text-[10px] text-amber-500/70 leading-tight">
{toyResult.circuit_diagram}
                                    </pre>
                                  </div>
                                </div>
                              )}
                            </div>
                          </motion.div>
                        )}

                        {/* 3. Decryption */}
                        {toyResult && (
                          <motion.div initial={{opacity:0}} animate={{opacity:1}} className="bg-slate-900 border border-slate-800 p-6 rounded-2xl relative shadow-inner">
                            <h3 className="text-sm font-bold text-emerald-400 mb-5 flex items-center gap-2">
                              <div className="w-6 h-6 rounded bg-emerald-500/20 flex items-center justify-center text-xs">3</div> 
                              Client-Side Decryption
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
                                <div className="text-lg tracking-widest font-bold">{xorBits(toyResult.measured_state, pocEncryptionState.key)}</div>
                                
                                {xorBits(toyResult.measured_state, pocEncryptionState.key) === pocEncryptionState.bits && (
                                  <div className="absolute top-2 right-2 flex space-x-1">
                                    <span className="w-1.5 h-1.5 rounded-full bg-emerald-500 animate-ping"></span>
                                    <span className="w-1.5 h-1.5 rounded-full bg-emerald-500"></span>
                                  </div>
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
          )}
        </motion.div>
      </AnimatePresence>
    </div>
  );
}
