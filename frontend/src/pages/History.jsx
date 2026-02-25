import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { History as HistoryIcon, Activity, Server, Clock, Search, ChevronRight } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

const API_BASE = 'http://localhost:8000/api';

export default function History() {
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const fetchHistory = async () => {
      try {
        const token = localStorage.getItem('token');
        const res = await axios.get(`${API_BASE}/history`, {
          headers: { Authorization: `Bearer ${token}` }
        });
        setHistory(res.data);
      } catch (err) {
        console.error("Failed to fetch history");
      } finally {
        setLoading(false);
      }
    };
    fetchHistory();
  }, []);

  const formatDate = (dateString) => {
    const d = new Date(dateString);
    return new Intl.DateTimeFormat('en-US', {
      month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit'
    }).format(d);
  };

  return (
    <div className="max-w-5xl mx-auto p-6 md:p-8">
      <div className="mb-8 border-b border-slate-800 pb-6 flex items-center gap-4">
        <div className="w-12 h-12 rounded-2xl bg-purple-500/20 flex items-center justify-center">
          <HistoryIcon className="w-6 h-6 text-purple-400" />
        </div>
        <div>
          <h1 className="text-3xl font-bold text-white mb-1">Simulation History</h1>
          <p className="text-slate-400">View past quantum proofs of concept and classical searches.</p>
        </div>
      </div>

      <div className="space-y-4">
        <AnimatePresence>
          {loading ? (
            <motion.div 
              initial={{opacity: 0}} animate={{opacity: 1}} exit={{opacity: 0}}
              className="flex justify-center p-12"
            >
              <div className="w-8 h-8 border-4 border-purple-500/30 border-t-purple-500 rounded-full animate-spin"></div>
            </motion.div>
          ) : history.length === 0 ? (
            <motion.div 
              initial={{opacity: 0, scale: 0.95}} animate={{opacity: 1, scale: 1}}
              className="bg-slate-900/50 border border-slate-800 p-12 rounded-3xl text-center"
            >
              <HistoryIcon className="w-12 h-12 text-slate-700 mx-auto mb-4" />
              <h3 className="text-lg font-medium text-slate-300 mb-2">No History Found</h3>
              <p className="text-slate-500">Run a quantum simulation from the Dashboard to see it here.</p>
            </motion.div>
          ) : (
            history.map((item, i) => (
              <motion.div
                key={item._id || i}
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: i * 0.05 }}
                className="bg-slate-900/60 border border-slate-800 rounded-2xl p-5 hover:bg-slate-800/80 transition-colors group cursor-default"
              >
                <div className="flex flex-col md:flex-row md:items-center justify-between gap-4">
                  <div className="flex items-center gap-4">
                    <div className={`w-10 h-10 rounded-xl flex items-center justify-center shrink-0 ${
                      item.type === 'quantum_ibm_submit' ? 'bg-blue-500/20 text-blue-400' : 'bg-amber-500/20 text-amber-400'
                    }`}>
                      {item.type === 'quantum_ibm_submit' ? <Server className="w-5 h-5" /> : <Activity className="w-5 h-5" />}
                    </div>
                    <div>
                      <h4 className="text-white font-medium flex items-center gap-2">
                        {item.type === 'quantum_ibm_submit' ? 'IBM Cloud Execution' : 'Local Qiskit Simulator'}
                        {item.status && (
                          <span className={`text-[10px] px-2 py-0.5 rounded-full uppercase tracking-wider font-bold ${
                            item.status === 'DONE' ? 'bg-emerald-500/20 text-emerald-400' : 'bg-amber-500/20 text-amber-400'
                          }`}>
                            {item.status}
                          </span>
                        )}
                      </h4>
                      <p className="text-xs text-slate-500 flex items-center gap-1 mt-1">
                        <Clock className="w-3 h-3" /> {formatDate(item.timestamp)}
                      </p>
                    </div>
                  </div>

                  <div className="flex items-center gap-6 md:ml-auto">
                    <div className="text-right">
                      <div className="text-xs text-slate-500 mb-1">Target Bits</div>
                      <div className="font-mono text-white bg-slate-950 px-2 py-1 rounded inline-block tracking-widest text-sm">
                        {item.target_bits || 'N/A'}
                      </div>
                    </div>
                    
                    {item.measured_state && (
                      <div className="text-right">
                        <div className="text-xs text-slate-500 mb-1">Result State</div>
                        <div className="font-mono text-emerald-400 bg-emerald-900/10 border border-emerald-500/20 px-2 py-1 rounded inline-block tracking-widest text-sm font-bold">
                          {item.measured_state}
                        </div>
                      </div>
                    )}
                    
                    {item.job_id && (
                      <div className="text-right hidden sm:block">
                        <div className="text-xs text-slate-500 mb-1">Job ID</div>
                        <div className="font-mono text-slate-400 text-xs w-24 truncate" title={item.job_id}>
                          {item.job_id}
                        </div>
                      </div>
                    )}
                  </div>
                </div>
              </motion.div>
            ))
          )}
        </AnimatePresence>
      </div>
    </div>
  );
}
