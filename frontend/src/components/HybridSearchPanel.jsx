import React, { useState } from 'react';
import axios from 'axios';
import { Search, Server } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import ComplexityChart from './ComplexityChart';

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000/api';
const getAuthHeaders = () => ({ headers: { Authorization: `Bearer ${localStorage.getItem('token')}` } });

export default function HybridSearchPanel() {
  const [accession, setAccession] = useState('NM_007294');
  const [targetKmer, setTargetKmer] = useState('ATGCT');
  const [datasetData, setDatasetData] = useState(null);
  const [isFetching, setIsFetching] = useState(false);
  const [isSearching, setIsSearching] = useState(false);
  const [classicalResult, setClassicalResult] = useState(null);
  const [quantumMetrics, setQuantumMetrics] = useState(null);

  const handleFetchDataset = async () => {
    setIsFetching(true);
    setDatasetData(null);
    setClassicalResult(null);
    setQuantumMetrics(null);
    try {
      const res = await axios.get(`${API_BASE}/sequence/${accession}`, getAuthHeaders());
      setDatasetData(res.data);
    } catch (err) {
      const msg = err.response?.data?.detail || 'Error fetching dataset. Check ID or network connection.';
      alert(msg);
    }
    setIsFetching(false);
  };

  const handleSearch = async () => {
    if (!datasetData) return alert('Fetch a dataset first!');
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
      alert('Error running search.');
    }
    setIsSearching(false);
  };

  return (
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
              <label className="block text-sm font-medium text-slate-400 mb-1.5">NCBI Accession ID</label>
              <div className="flex gap-2">
                <input
                  type="text"
                  value={accession}
                  onChange={(e) => setAccession(e.target.value.trim())}
                  className="w-full bg-slate-950/50 border border-slate-700/50 rounded-xl px-4 py-2.5 text-sm focus:border-blue-500/50 focus:ring-2 focus:ring-blue-500/20 outline-none transition-all font-mono"
                  placeholder="e.g. NM_007294"
                />
                <button
                  onClick={handleFetchDataset}
                  disabled={isFetching}
                  className="bg-blue-600 hover:bg-blue-500 px-5 py-2.5 rounded-xl font-medium transition-colors disabled:opacity-50 whitespace-nowrap"
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
                    <span className="w-2 h-2 rounded-full bg-emerald-400 animate-pulse" /> Data Loaded
                  </div>
                  <div className="text-xs text-slate-400">
                    Total Bases: <span className="text-white font-mono">{datasetData.length.toLocaleString()}</span> bp
                  </div>
                  <div className="text-xs text-slate-500 mt-2 truncate font-mono bg-black/20 p-2 rounded">
                    {datasetData.sequence.substring(0, 40)}...
                  </div>
                </motion.div>
              )}
            </AnimatePresence>

            <div>
              <label className="block text-sm font-medium text-slate-400 mb-1.5">Target k-mer Pattern</label>
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
            <div className="absolute inset-0 bg-white/20 translate-y-full group-hover:translate-y-0 transition-transform" />
            <Search className="w-5 h-5 relative z-10" />
            <span className="relative z-10">{isSearching ? 'Computing...' : 'Run Hybrid Search'}</span>
          </button>
        </div>
      </div>

      {/* RIGHT COLUMN: Results */}
      <div className="xl:col-span-8 flex flex-col gap-6">
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {/* Classical Result */}
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
                    <div className="text-xl font-semibold text-slate-300">
                      {classicalResult.execution_time_ms.toFixed(1)} <span className="text-sm text-slate-500">ms</span>
                    </div>
                  </div>
                </div>
              </div>
            ) : (
              <div className="h-48 flex items-center justify-center text-slate-600/70 text-sm">
                Run a search to view classical metrics
              </div>
            )}
          </div>

          {/* Quantum Theory */}
          <div className="bg-slate-900/60 rounded-3xl p-6 border border-slate-800 shadow-2xl relative overflow-hidden group">
            <div className="absolute -bottom-10 -right-10 w-48 h-48 bg-purple-500/10 rounded-full blur-[50px] group-hover:bg-purple-500/20 transition-colors" />
            <h3 className="text-lg font-semibold text-white mb-1">Quantum Theory</h3>
            <p className="text-xs text-slate-400 mb-1">Theoretical Grover's oracle calls O(√N)</p>
            <p className="text-xs text-amber-500/70 mb-6 italic">
              ⚠ Analytically derived — not measured from a circuit of size N
            </p>
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
                    <span className="text-sm font-normal text-purple-200/60">fewer oracle calls</span>
                  </div>
                  <div className="text-xs text-slate-400">
                    Classical requires {classicalResult.classical_queries.toLocaleString()} comparisons
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

        {/* Chart */}
        <div className="bg-slate-900/60 rounded-3xl p-6 border border-slate-800 shadow-2xl flex-1 flex flex-col min-h-[400px]">
          <div className="mb-6">
            <h3 className="text-lg font-semibold text-white">Complexity Divergence</h3>
            <p className="text-xs text-slate-500 mt-1">
              Classical O(N) vs. Theoretical Grover O(√N) — values are analytically computed, not experimentally measured
            </p>
          </div>
          <ComplexityChart classicalResult={classicalResult} quantumMetrics={quantumMetrics} />
        </div>
      </div>
    </div>
  );
}
