import React, { useState, useEffect } from 'react';
import {
  History as HistoryIcon,
  Activity,
  Server,
  Clock,
  ChevronDown,
  ChevronRight,
  Cpu,
} from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import api from '../lib/api';
import GroverStepNavigator from '../components/GroverStepNavigator';

export default function History() {
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(true);
  const [expandedId, setExpandedId] = useState(null);
  const [circuitCache, setCircuitCache] = useState({}); // keyed by target_bits
  const [loadingCircuit, setLoadingCircuit] = useState(null);
  const [circuitError, setCircuitError] = useState({}); // keyed by target_bits
  const [activeStep, setActiveStep] = useState({}); // { [itemId]: stepIndex }

  useEffect(() => {
    const fetchHistory = async () => {
      try {
        const res = await api.get('/history');
        setHistory(res.data);
      } catch (err) {
        console.error('Failed to fetch history');
      } finally {
        setLoading(false);
      }
    };
    fetchHistory();
  }, []);

  const formatDate = (dateString) => {
    const d = new Date(dateString);
    return new Intl.DateTimeFormat('en-US', {
      month: 'short',
      day: 'numeric',
      hour: '2-digit',
      minute: '2-digit',
    }).format(d);
  };

  const handleToggle = async (item) => {
    const id = item._id || item.timestamp;
    if (expandedId === id) {
      setExpandedId(null);
      return;
    }
    setExpandedId(id);

    const bits = item.target_bits;
    if (!bits || circuitCache[bits]) return; // already cached or no bits

    setLoadingCircuit(id);
    setCircuitError((prev) => ({ ...prev, [bits]: false }));
    try {
      const res = await api.get(`/search/circuit-info?target_bits=${bits}`);
      setCircuitCache((prev) => ({ ...prev, [bits]: res.data }));
    } catch (err) {
      console.error('Failed to fetch circuit info', err);
      setCircuitError((prev) => ({ ...prev, [bits]: true }));
    } finally {
      setLoadingCircuit(null);
    }
  };

  const getStep = (id) => activeStep[id] ?? 0;
  const setStep = (id, idx) =>
    setActiveStep((prev) => ({ ...prev, [id]: idx }));

  return (
    <div className="max-w-5xl mx-auto p-6 md:p-8">
      <div className="mb-8 border-b border-slate-800 pb-6 flex items-center gap-4">
        <div className="w-12 h-12 rounded-2xl bg-purple-500/20 flex items-center justify-center">
          <HistoryIcon className="w-6 h-6 text-purple-400" />
        </div>
        <div>
          <h1 className="text-3xl font-bold text-white mb-1">
            Simulation History
          </h1>
          <p className="text-slate-400">
            Click any entry to visualize its quantum circuit and amplitude
            chart.
          </p>
        </div>
      </div>

      <div className="space-y-3">
        <AnimatePresence>
          {loading ? (
            <motion.div
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              exit={{ opacity: 0 }}
              className="flex justify-center p-12">
              <div className="w-8 h-8 border-4 border-purple-500/30 border-t-purple-500 rounded-full animate-spin" />
            </motion.div>
          ) : history.length === 0 ? (
            <motion.div
              initial={{ opacity: 0, scale: 0.95 }}
              animate={{ opacity: 1, scale: 1 }}
              className="bg-slate-900/50 border border-slate-800 p-12 rounded-3xl text-center">
              <HistoryIcon className="w-12 h-12 text-slate-700 mx-auto mb-4" />
              <h3 className="text-lg font-medium text-slate-300 mb-2">
                No History Found
              </h3>
              <p className="text-slate-500">
                Run a quantum simulation from the Dashboard to see it here.
              </p>
            </motion.div>
          ) : (
            history.map((item, i) => {
              const id = item._id || item.timestamp;
              const isExpanded = expandedId === id;
              const bits = item.target_bits;
              const circuit = bits ? circuitCache[bits] : null;
              const stepIdx = getStep(id);

              return (
                <motion.div
                  key={id}
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: i * 0.04 }}
                  className="bg-slate-900/60 border border-slate-800 rounded-2xl overflow-hidden">
                  {/* ── Summary Row (clickable) ── */}
                  <button
                    onClick={() => handleToggle(item)}
                    className="w-full text-left p-5 hover:bg-slate-800/60 transition-colors group">
                    <div className="flex flex-col md:flex-row md:items-center justify-between gap-4">
                      <div className="flex items-center gap-4">
                        <div
                          className={`w-10 h-10 rounded-xl flex items-center justify-center shrink-0 ${
                            item.type === 'quantum_ibm_submit'
                              ? 'bg-blue-500/20 text-blue-400'
                              : 'bg-amber-500/20 text-amber-400'
                          }`}>
                          {item.type === 'quantum_ibm_submit' ? (
                            <Server className="w-5 h-5" />
                          ) : (
                            <Activity className="w-5 h-5" />
                          )}
                        </div>
                        <div>
                          <h4 className="text-white font-medium flex items-center gap-2">
                            {item.type === 'quantum_ibm_submit'
                              ? 'IBM Cloud Execution'
                              : 'Local Qiskit Simulator'}
                            {item.status && (
                              <span
                                className={`text-[10px] px-2 py-0.5 rounded-full uppercase tracking-wider font-bold ${
                                  item.status === 'DONE'
                                    ? 'bg-emerald-500/20 text-emerald-400'
                                    : 'bg-amber-500/20 text-amber-400'
                                }`}>
                                {item.status}
                              </span>
                            )}
                          </h4>
                          <p className="text-xs text-slate-500 flex items-center gap-1 mt-1">
                            <Clock className="w-3 h-3" />{' '}
                            {formatDate(item.timestamp)}
                          </p>
                        </div>
                      </div>

                      <div className="flex items-center gap-5 md:ml-auto">
                        <div className="text-right">
                          <div className="text-xs text-slate-500 mb-1">
                            Target Bits
                          </div>
                          <div className="font-mono text-white bg-slate-950 px-2 py-1 rounded inline-block tracking-widest text-sm">
                            {item.target_bits || 'N/A'}
                          </div>
                        </div>

                        {item.measured_state && (
                          <div className="text-right">
                            <div className="text-xs text-slate-500 mb-1">
                              Result State
                            </div>
                            <div className="font-mono text-emerald-400 bg-emerald-900/10 border border-emerald-500/20 px-2 py-1 rounded inline-block tracking-widest text-sm font-bold">
                              {item.measured_state}
                            </div>
                          </div>
                        )}

                        {item.job_id && (
                          <div className="text-right hidden sm:block">
                            <div className="text-xs text-slate-500 mb-1">
                              Job ID
                            </div>
                            <div
                              className="font-mono text-slate-400 text-xs w-24 truncate"
                              title={item.job_id}>
                              {item.job_id}
                            </div>
                          </div>
                        )}

                        <div
                          className={`ml-2 transition-transform duration-200 ${isExpanded ? 'rotate-90' : ''} text-slate-500 group-hover:text-slate-300`}>
                          <ChevronRight className="w-4 h-4" />
                        </div>
                      </div>
                    </div>
                  </button>

                  {/* ── Expandable Circuit Panel ── */}
                  <AnimatePresence initial={false}>
                    {isExpanded && (
                      <motion.div
                        key="panel"
                        initial={{ opacity: 0, height: 0 }}
                        animate={{ opacity: 1, height: 'auto' }}
                        exit={{ opacity: 0, height: 0 }}
                        transition={{ duration: 0.2 }}
                        className="overflow-hidden">
                        <div className="border-t border-slate-800 p-5 space-y-5 bg-slate-950/40">
                          {/* Loading state */}
                          {loadingCircuit === id && (
                            <div className="flex items-center gap-3 justify-center py-8 text-slate-400 text-sm">
                              <div className="w-5 h-5 border-2 border-amber-500/30 border-t-amber-500 rounded-full animate-spin" />
                              Building circuit visualization…
                            </div>
                          )}

                          {/* No target bits */}
                          {!bits && !loadingCircuit && (
                            <p className="text-slate-500 text-sm text-center py-6">
                              No target bits recorded for this entry — circuit
                              cannot be reconstructed.
                            </p>
                          )}

                          {/* Circuit fetch error */}
                          {circuitError[bits] &&
                            !loadingCircuit &&
                            !circuit && (
                              <p className="text-red-400 text-sm text-center py-6">
                                Failed to load circuit visualization. Please try
                                again.
                              </p>
                            )}

                          {/* Circuit data available */}
                          {circuit && (
                            <>
                              {/* Stats row */}
                              <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 text-center text-sm">
                                {[
                                  {
                                    label: 'Qubits',
                                    value: circuit.qubits,
                                    color: 'text-blue-400',
                                  },
                                  {
                                    label: 'Iterations',
                                    value: circuit.iterations,
                                    color: 'text-purple-400',
                                  },
                                  {
                                    label: 'Noise',
                                    value:
                                      item.noise_level != null
                                        ? `${(item.noise_level * 100).toFixed(1)}%`
                                        : '—',
                                    color: 'text-amber-400',
                                  },
                                  {
                                    label: 'Time (ms)',
                                    value:
                                      item.execution_time_ms != null
                                        ? item.execution_time_ms.toFixed(1)
                                        : '—',
                                    color: 'text-emerald-400',
                                  },
                                ].map((s) => (
                                  <div
                                    key={s.label}
                                    className="bg-slate-900 border border-slate-800 rounded-xl p-3">
                                    <div className="text-slate-500 text-xs uppercase tracking-wider mb-1">
                                      {s.label}
                                    </div>
                                    <div
                                      className={`font-mono font-bold ${s.color}`}>
                                      {s.value}
                                    </div>
                                  </div>
                                ))}
                              </div>

                              {/* Step navigator — shared with Quantum Sim POC */}
                              <GroverStepNavigator
                                stepCircuits={circuit.step_circuits}
                                activeStep={stepIdx}
                                onStepChange={(idx) => setStep(id, idx)}
                                targetBits={bits}
                              />

                              {/* Full circuit diagram (collapsible) */}
                              <details className="group">
                                <summary className="cursor-pointer flex items-center gap-2 text-xs text-slate-500 hover:text-slate-300 transition-colors select-none list-none">
                                  <Cpu className="w-3.5 h-3.5" />
                                  Full circuit diagram
                                  <ChevronDown className="w-3.5 h-3.5 transition-transform group-open:rotate-180" />
                                </summary>
                                <div className="mt-3 overflow-x-auto border border-slate-800 rounded-xl bg-slate-950 p-4">
                                  <pre className="text-[10px] text-emerald-500/60 leading-tight whitespace-pre">
                                    {circuit.circuit_diagram}
                                  </pre>
                                </div>
                              </details>
                            </>
                          )}
                        </div>
                      </motion.div>
                    )}
                  </AnimatePresence>
                </motion.div>
              );
            })
          )}
        </AnimatePresence>
      </div>
    </div>
  );
}
