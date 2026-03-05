import React, { useState, useEffect, useRef, useCallback, useMemo } from 'react';
import {
  History as HistoryIcon,
  Activity,
  Clock,
  ChevronDown,
  ChevronRight,
  Cpu,
  Dna,
  Loader2,
} from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';
import api from '../lib/api';
import GroverStepNavigator from '../components/GroverStepNavigator';

// ─── Helpers ────────────────────────────────────────────────────────────────

/**
 * Groups raw history items so that bio_grover_local pairs sharing the same
 * run_id are merged into a single bio_pair display group.
 * Legacy bio items (no run_id) and all other types pass through as singles.
 */
function groupHistory(rawItems) {
  const result  = [];
  const pending = new Map(); // run_id → { idx, data }

  for (const item of rawItems) {
    if (item.type === 'bio_grover_local' && item.run_id) {
      if (pending.has(item.run_id)) {
        const { idx, data } = pending.get(item.run_id);
        if (item.has_mutation) data.carrier = item;
        else data.healthy = item;
        // Prefer whichever item has target_bits (both should match, but guard for safety)
        if (!data.target_bits && item.target_bits) data.target_bits = item.target_bits;
        if (!data.marker_variant && item.marker_variant) data.marker_variant = item.marker_variant;
        if (!data.marker_dna && item.marker_dna) data.marker_dna = item.marker_dna;
        result[idx] = data;
        pending.delete(item.run_id);
      } else {
        const data = {
          type:           'bio_pair',
          id:             item.run_id,
          run_id:         item.run_id,
          timestamp:      item.timestamp,
          carrier:        item.has_mutation ? item : null,
          healthy:        item.has_mutation ? null  : item,
          n_codons:       item.n_codons,
          marker_variant: item.marker_variant,
          marker_dna:     item.marker_dna,
          target_bits:    item.target_bits,
        };
        pending.set(item.run_id, { idx: result.length, data });
        result.push(data);
      }
    } else {
      result.push({ type: 'single', id: item._id || item.timestamp, item });
    }
  }
  return result;
}

const bioCacheKey = (target_bits, has_mutation) =>
  `${target_bits}:${has_mutation}`;

// ─── Component ───────────────────────────────────────────────────────────────

export default function History() {
  const [rawItems,     setRawItems]    = useState([]);
  const [loading,      setLoading]     = useState(true);
  const [loadingMore,  setLoadingMore] = useState(false);
  const [hasMore,      setHasMore]     = useState(false);
  const [skip,         setSkip]        = useState(0);
  const [expandedId,   setExpandedId]  = useState(null);

  // displayGroups is derived from rawItems — avoids nested setState bugs
  const displayGroups = useMemo(() => groupHistory(rawItems), [rawItems]);

  // Standard Grover circuit cache: keyed by target_bits
  const [circuitCache,    setCircuitCache]    = useState({});
  // Bio circuit cache: keyed by bioCacheKey(target_bits, has_mutation)
  const [bioCircuitCache, setBioCircuitCache] = useState({});
  const [loadingCircuit,  setLoadingCircuit]  = useState(null);
  const [circuitError,    setCircuitError]    = useState({});

  // Step index per card id
  const [activeStep,    setActiveStep]    = useState({});
  // Active patient per bio_pair group id: 'carrier' | 'healthy'
  const [activePatient, setActivePatient] = useState({});

  const sentinelRef = useRef(null);

  // ── Data fetching ──────────────────────────────────────────────────────────

  const fetchPage = useCallback(async (skipVal) => {
    try {
      const res = await api.get(`/history?skip=${skipVal}&limit=20`);
      const data = res.data;
      // Guard: backend may still be old (returns flat array) during rolling deploy
      const items    = Array.isArray(data) ? data : (Array.isArray(data?.items) ? data.items : []);
      const has_more = Array.isArray(data) ? false : Boolean(data?.has_more);
      setRawItems((prev) => skipVal === 0 ? items : [...prev, ...items]);
      setHasMore(has_more);
      setSkip(skipVal + items.length);
    } catch (err) {
      console.error('Failed to fetch history', err);
    }
  }, []);

  useEffect(() => {
    fetchPage(0).finally(() => setLoading(false));
  }, [fetchPage]);

  // IntersectionObserver — load more when sentinel scrolls into view
  useEffect(() => {
    const el = sentinelRef.current;
    if (!el) return;
    const observer = new IntersectionObserver(
      ([entry]) => {
        if (entry.isIntersecting && hasMore && !loadingMore) {
          setLoadingMore(true);
          fetchPage(skip).finally(() => setLoadingMore(false));
        }
      },
      { rootMargin: '200px' },
    );
    observer.observe(el);
    return () => observer.disconnect();
  }, [hasMore, loadingMore, skip, fetchPage]);

  // ── Utilities ──────────────────────────────────────────────────────────────

  const formatDate = (dateString) => {
    if (!dateString) return '—';
    const hasZone = /Z$|[+-]\d{2}:?\d{2}$/.test(dateString);
    const utcStr  = hasZone ? dateString : dateString + 'Z';
    const d       = new Date(utcStr);
    if (isNaN(d.getTime())) return dateString;
    return new Intl.DateTimeFormat(undefined, {
      year: 'numeric', month: 'short', day: 'numeric',
      hour: '2-digit', minute: '2-digit', timeZoneName: 'short',
    }).format(d);
  };

  const getStep    = (id)       => activeStep[id] ?? 0;
  const setStep    = (id, idx)  => setActiveStep((p) => ({ ...p, [id]: idx }));
  const getPatient = (id)       => activePatient[id] ?? 'carrier';
  const setPatient = (id, val)  => setActivePatient((p) => ({ ...p, [id]: val }));

  // ── Expand / circuit-fetch ─────────────────────────────────────────────────

  const handleToggle = async (group) => {
    const id = group.id;
    if (expandedId === id) { setExpandedId(null); return; }
    setExpandedId(id);

    if (group.type === 'bio_pair') {
      const bits    = group.target_bits;
      const nCodons = group.n_codons ?? 1;
      if (!bits) return;
      const toFetch = [true, false].filter(
        (hm) => !bioCircuitCache[bioCacheKey(bits, hm)],
      );
      if (toFetch.length === 0) return;
      setLoadingCircuit(id);
      await Promise.all(
        toFetch.map(async (hasMutation) => {
          const ck = bioCacheKey(bits, hasMutation);
          try {
            const res = await api.get(
              `/search/quantum-poc/bio-circuit-info?target_bits=${bits}&n_codons=${nCodons}&has_mutation=${hasMutation}`,
            );
            setBioCircuitCache((prev) => ({ ...prev, [ck]: res.data.step_circuits }));
          } catch {
            setCircuitError((prev) => ({ ...prev, [ck]: true }));
          }
        }),
      );
      setLoadingCircuit(null);
      return;
    }

    // Single item
    const item   = group.item;
    const isBio  = item.type === 'bio_grover_local';
    const bits   = item.target_bits;

    if (isBio && bits) {
      const ck = bioCacheKey(bits, item.has_mutation ?? true);
      if (bioCircuitCache[ck]) return;
      setLoadingCircuit(id);
      try {
        const res = await api.get(
          `/search/quantum-poc/bio-circuit-info?target_bits=${bits}&n_codons=${item.n_codons ?? 1}&has_mutation=${item.has_mutation ?? true}`,
        );
        setBioCircuitCache((prev) => ({ ...prev, [ck]: res.data.step_circuits }));
      } catch {
        setCircuitError((prev) => ({ ...prev, [ck]: true }));
      } finally {
        setLoadingCircuit(null);
      }
      return;
    }

    if (!isBio && bits) {
      if (circuitCache[bits]) return;
      setLoadingCircuit(id);
      setCircuitError((prev) => ({ ...prev, [bits]: false }));
      try {
        const res = await api.get(`/search/circuit-info?target_bits=${bits}`);
        setCircuitCache((prev) => ({ ...prev, [bits]: res.data }));
      } catch {
        setCircuitError((prev) => ({ ...prev, [bits]: true }));
      } finally {
        setLoadingCircuit(null);
      }
    }
  };

  // ── Sub-renderers ──────────────────────────────────────────────────────────

  function StatsGrid({ stats }) {
    return (
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 text-center text-sm">
        {stats.map((s) => (
          <div key={s.label} className="bg-slate-900 border border-slate-800 rounded-xl p-3">
            <div className="text-slate-500 text-xs uppercase tracking-wider mb-1">{s.label}</div>
            <div className={`font-mono font-bold text-sm ${s.color}`}>{s.value}</div>
          </div>
        ))}
      </div>
    );
  }

  function MarkerInfo({ item }) {
    return (
      <div className="bg-slate-900/60 border border-slate-800 rounded-xl p-4 space-y-2 text-sm">
        <div className="text-xs text-slate-500 uppercase tracking-widest mb-2">BRCA1 Marker</div>
        <div className="flex justify-between">
          <span className="text-slate-400">Variant</span>
          <span className="font-mono text-amber-400 font-bold">{item?.marker_variant || '—'}</span>
        </div>
        {item?.marker_dna && (
          <div className="flex justify-between">
            <span className="text-slate-400">Marker DNA</span>
            <span className="font-mono text-emerald-400 tracking-widest">{item.marker_dna}</span>
          </div>
        )}
        {item?.target_bits && (
          <div className="flex justify-between">
            <span className="text-slate-400">Marker Bits</span>
            <span className="font-mono text-slate-300 text-xs tracking-widest">{item.target_bits}</span>
          </div>
        )}
      </div>
    );
  }

  function BioCircuitSection({ bits, nCodons, hasMutation, groupId }) {
    const ck      = bioCacheKey(bits, hasMutation);
    const circuits = bioCircuitCache[ck];
    const stepIdx  = getStep(`${groupId}:${hasMutation}`);
    const isLoading = loadingCircuit === groupId;
    const hasError  = circuitError[ck];

    if (isLoading && !circuits) {
      return (
        <div className="flex items-center gap-3 justify-center py-6 text-slate-400 text-sm">
          <Loader2 className="w-5 h-5 animate-spin text-emerald-500" />
          Building circuit visualization…
        </div>
      );
    }
    if (hasError && !circuits) {
      return (
        <p className="text-red-400 text-sm text-center py-4">
          Failed to load circuit. Please try again.
        </p>
      );
    }
    if (!circuits) return null;

    return (
      <>
        <GroverStepNavigator
          stepCircuits={circuits}
          activeStep={stepIdx}
          onStepChange={(idx) => setStep(`${groupId}:${hasMutation}`, idx)}
          targetBits={bits}
        />
        <details className="group">
          <summary className="cursor-pointer flex items-center gap-2 text-xs text-slate-500 hover:text-slate-300 transition-colors select-none list-none">
            <Cpu className="w-3.5 h-3.5" />
            Full circuit diagram
            <ChevronDown className="w-3.5 h-3.5 transition-transform group-open:rotate-180" />
          </summary>
          <div className="mt-3 overflow-x-auto border border-slate-800 rounded-xl bg-slate-950 p-4">
            <pre className="text-[10px] text-emerald-500/60 leading-tight whitespace-pre">
              {circuits[stepIdx]?.diagram || '—'}
            </pre>
          </div>
        </details>
      </>
    );
  }

  // ── Render bio_pair expanded panel ─────────────────────────────────────────

  function BioPairPanel({ group }) {
    const pat     = getPatient(group.id);
    const item    = pat === 'carrier' ? group.carrier : group.healthy;
    const nCodons = group.n_codons ?? 1;
    const bits    = group.target_bits;

    return (
      <div className="space-y-4">
        {/* Patient selector tabs */}
        <div className="grid grid-cols-2 gap-3">
          {[
            { key: 'carrier', label: 'Patient A — Carrier',       data: group.carrier, found: group.carrier?.detection_result === 'FOUND', accent: 'red' },
            { key: 'healthy', label: 'Patient B — Healthy Control', data: group.healthy, found: group.healthy?.detection_result === 'FOUND', accent: 'emerald' },
          ].map(({ key, label, data, found, accent }) => (
            <button
              key={key}
              onClick={() => setPatient(group.id, key)}
              className={`rounded-2xl border p-4 text-left transition-all ${
                pat === key ? 'ring-2 ring-offset-2 ring-offset-slate-900' : 'opacity-70 hover:opacity-100'
              } ${found
                ? `bg-red-900/20 border-red-500/40 ${pat === key ? 'ring-red-500' : ''}`
                : `bg-emerald-900/10 border-emerald-500/30 ${pat === key ? 'ring-emerald-500' : ''}`
              }`}>
              <div className="text-xs font-medium text-slate-400 mb-1 uppercase tracking-wider">{label}</div>
              {data ? (
                <>
                  <div className={`text-base font-bold ${found ? 'text-red-400' : 'text-emerald-400'}`}>
                    {found ? '🧬 DETECTED' : '✅ NOT FOUND'}
                  </div>
                  <div className="text-xs text-slate-500 mt-1">
                    Confidence: <span className="font-mono text-white">{((data.confidence ?? 0) * 100).toFixed(1)}%</span>
                  </div>
                </>
              ) : (
                <div className="text-xs text-slate-600 italic">No data</div>
              )}
            </button>
          ))}
        </div>

        {/* Stats grid for selected patient */}
        {item && (
          <StatsGrid stats={[
            { label: 'Patient',   value: item.has_mutation ? 'Carrier' : 'Healthy Control', color: item.has_mutation ? 'text-red-400' : 'text-emerald-400' },
            { label: 'Qubits',    value: item.n_qubits ?? (nCodons * 6), color: 'text-blue-400' },
            { label: 'n_codons',  value: item.n_codons ?? nCodons,         color: 'text-purple-400' },
            { label: 'Time (ms)', value: item.execution_time_ms != null ? item.execution_time_ms.toFixed(1) : '—', color: 'text-amber-400' },
          ]} />
        )}

        {/* Marker info */}
        <MarkerInfo item={group.carrier ?? group.healthy} />

        {/* Circuit visualization */}
        {bits && (
          <div>
            <div className="text-xs text-slate-500 uppercase tracking-widest mb-2 flex items-center gap-1">
              <Dna className="w-3.5 h-3.5" />
              {pat === 'carrier' ? 'Carrier' : 'Healthy'} Patient — Grover Step Circuits
            </div>
            <BioCircuitSection
              bits={bits}
              nCodons={nCodons}
              hasMutation={pat === 'carrier'}
              groupId={group.id}
            />
          </div>
        )}
      </div>
    );
  }

  // ── Render single bio item expanded panel ───────────────────────────────────

  function SingleBioPanel({ item, groupId }) {
    const bits    = item.target_bits;
    const nCodons = item.n_codons ?? 1;
    const hm      = item.has_mutation ?? true;

    return (
      <div className="space-y-4">
        {item.detection_result && (
          <div className={`rounded-2xl border p-4 flex items-center gap-4 ${
            item.detection_result === 'FOUND'
              ? 'bg-red-900/20 border-red-500/40'
              : 'bg-emerald-900/20 border-emerald-500/40'
          }`}>
            <div className={`text-3xl font-black ${item.detection_result === 'FOUND' ? 'text-red-400' : 'text-emerald-400'}`}>
              {item.detection_result === 'FOUND' ? '🧬 DETECTED' : '✅ NOT FOUND'}
            </div>
            {item.confidence != null && (
              <div className="ml-auto text-right">
                <div className="text-xs text-slate-400 mb-0.5">Confidence</div>
                <div className="font-mono font-bold text-white text-lg">
                  {(item.confidence * 100).toFixed(1)}%
                </div>
              </div>
            )}
          </div>
        )}
        <StatsGrid stats={[
          { label: 'Patient',   value: hm ? 'Carrier' : 'Healthy Control', color: hm ? 'text-red-400' : 'text-emerald-400' },
          { label: 'Qubits',    value: item.n_qubits ?? (nCodons * 6), color: 'text-blue-400' },
          { label: 'n_codons',  value: item.n_codons ?? nCodons,         color: 'text-purple-400' },
          { label: 'Time (ms)', value: item.execution_time_ms != null ? item.execution_time_ms.toFixed(1) : '—', color: 'text-amber-400' },
        ]} />
        <MarkerInfo item={item} />
        {bits && (
          <div>
            <div className="text-xs text-slate-500 uppercase tracking-widest mb-2 flex items-center gap-1">
              <Dna className="w-3.5 h-3.5" /> Grover Step Circuits
            </div>
            <BioCircuitSection bits={bits} nCodons={nCodons} hasMutation={hm} groupId={groupId} />
          </div>
        )}
      </div>
    );
  }

  // ── Main render ────────────────────────────────────────────────────────────

  return (
    <div className="max-w-5xl mx-auto p-6 md:p-8">
      <div className="mb-8 border-b border-slate-800 pb-6 flex items-center gap-4">
        <div className="w-12 h-12 rounded-2xl bg-purple-500/20 flex items-center justify-center">
          <HistoryIcon className="w-6 h-6 text-purple-400" />
        </div>
        <div>
          <h1 className="text-3xl font-bold text-white mb-1">Simulation History</h1>
          <p className="text-slate-400">
            Click any entry to visualize its quantum circuit and amplitude chart.
          </p>
        </div>
      </div>

      <div className="space-y-3">
        <AnimatePresence>
          {loading ? (
            <motion.div
              initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}
              className="flex justify-center p-12">
              <div className="w-8 h-8 border-4 border-purple-500/30 border-t-purple-500 rounded-full animate-spin" />
            </motion.div>
          ) : displayGroups.length === 0 ? (
            <motion.div
              initial={{ opacity: 0, scale: 0.95 }} animate={{ opacity: 1, scale: 1 }}
              className="bg-slate-900/50 border border-slate-800 p-12 rounded-3xl text-center">
              <HistoryIcon className="w-12 h-12 text-slate-700 mx-auto mb-4" />
              <h3 className="text-lg font-medium text-slate-300 mb-2">No History Found</h3>
              <p className="text-slate-500">Run a quantum simulation from the Dashboard to see it here.</p>
            </motion.div>
          ) : (
            displayGroups.map((group, i) => {
              const id         = group.id;
              const isExpanded = expandedId === id;

              // ── bio_pair group ───────────────────────────────────────────
              if (group.type === 'bio_pair') {
                const carrierFound  = group.carrier?.detection_result === 'FOUND';
                const healthyFound  = group.healthy?.detection_result === 'FOUND';
                return (
                  <motion.div
                    key={id}
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: Math.min(i * 0.04, 0.5) }}
                    className="bg-slate-900/60 border border-slate-800 rounded-2xl overflow-hidden">
                    <button
                      onClick={() => handleToggle(group)}
                      className="w-full text-left p-5 hover:bg-slate-800/60 transition-colors group">
                      <div className="flex flex-col md:flex-row md:items-center justify-between gap-4">
                        <div className="flex items-center gap-4">
                          <div className="w-10 h-10 rounded-xl flex items-center justify-center shrink-0 bg-emerald-500/20 text-emerald-400">
                            <Dna className="w-5 h-5" />
                          </div>
                          <div>
                            <h4 className="text-white font-medium flex items-center gap-2">
                              BRCA1 Bio Simulation
                              <span className="text-[10px] px-2 py-0.5 rounded-full uppercase tracking-wider font-bold bg-emerald-500/20 text-emerald-400">
                                2 patients
                              </span>
                            </h4>
                            <p className="text-xs text-slate-500 flex items-center gap-1 mt-1">
                              <Clock className="w-3 h-3" /> {formatDate(group.timestamp)}
                            </p>
                          </div>
                        </div>
                        <div className="flex items-center gap-4 md:ml-auto">
                          {/* Carrier badge */}
                          <div className="text-right">
                            <div className="text-xs text-slate-500 mb-1">Patient A</div>
                            <div className={`text-xs font-bold px-2 py-1 rounded-full ${
                              carrierFound ? 'bg-red-500/20 text-red-400' : 'bg-emerald-500/20 text-emerald-400'
                            }`}>
                              {group.carrier ? (carrierFound ? '🧬 DETECTED' : '✅ CLEAN') : '—'}
                            </div>
                          </div>
                          {/* Healthy badge */}
                          <div className="text-right">
                            <div className="text-xs text-slate-500 mb-1">Patient B</div>
                            <div className={`text-xs font-bold px-2 py-1 rounded-full ${
                              healthyFound ? 'bg-red-500/20 text-red-400' : 'bg-emerald-500/20 text-emerald-400'
                            }`}>
                              {group.healthy ? (healthyFound ? '🧬 DETECTED' : '✅ CLEAN') : '—'}
                            </div>
                          </div>
                          {/* Marker variant */}
                          {group.marker_variant && (
                            <div className="text-right hidden sm:block">
                              <div className="text-xs text-slate-500 mb-1">Variant</div>
                              <div className="font-mono text-amber-400 text-xs">{group.marker_variant}</div>
                            </div>
                          )}
                          <div className={`ml-2 transition-transform duration-200 ${isExpanded ? 'rotate-90' : ''} text-slate-500 group-hover:text-slate-300`}>
                            <ChevronRight className="w-4 h-4" />
                          </div>
                        </div>
                      </div>
                    </button>

                    <AnimatePresence initial={false}>
                      {isExpanded && (
                        <motion.div
                          key="panel"
                          initial={{ opacity: 0, height: 0 }}
                          animate={{ opacity: 1, height: 'auto' }}
                          exit={{ opacity: 0, height: 0 }}
                          transition={{ duration: 0.2 }}
                          className="overflow-hidden">
                          <div className="border-t border-slate-800 p-5 bg-slate-950/40">
                            <BioPairPanel group={group} />
                          </div>
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </motion.div>
                );
              }

              // ── single item ──────────────────────────────────────────────
              const item   = group.item;
              const isBio  = item.type === 'bio_grover_local';
              const bits   = item.target_bits;
              const circuit    = !isBio && bits ? circuitCache[bits] : null;
              const stepIdx    = getStep(id);

              const typeLabel = {
                quantum_local_sim:     'Local Qiskit Simulator',
                bio_grover_local:      'BRCA1 Bio Simulation',
              }[item.type] ?? item.type;

              return (
                <motion.div
                  key={id}
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: Math.min(i * 0.04, 0.5) }}
                  className="bg-slate-900/60 border border-slate-800 rounded-2xl overflow-hidden">
                  {/* Summary row */}
                  <button
                    onClick={() => handleToggle(group)}
                    className="w-full text-left p-5 hover:bg-slate-800/60 transition-colors group">
                    <div className="flex flex-col md:flex-row md:items-center justify-between gap-4">
                      <div className="flex items-center gap-4">
                        <div className={`w-10 h-10 rounded-xl flex items-center justify-center shrink-0 ${
                           isBio ? 'bg-emerald-500/20 text-emerald-400'
                                : 'bg-amber-500/20 text-amber-400'
                        }`}>
                          {isBio ? <Dna className="w-5 h-5" /> : <Activity className="w-5 h-5" />}
                        </div>
                        <div>
                          <h4 className="text-white font-medium flex items-center gap-2">
                            {typeLabel}
                            {item.status && (
                              <span className={`text-[10px] px-2 py-0.5 rounded-full uppercase tracking-wider font-bold ${
                                item.status === 'DONE' ? 'bg-emerald-500/20 text-emerald-400' : 'bg-amber-500/20 text-amber-400'
                              }`}>{item.status}</span>
                            )}
                          </h4>
                          <p className="text-xs text-slate-500 flex items-center gap-1 mt-1">
                            <Clock className="w-3 h-3" /> {formatDate(item.timestamp)}
                          </p>
                        </div>
                      </div>

                      <div className="flex items-center gap-5 md:ml-auto">
                        {isBio && item.detection_result && (
                          <div className="text-right">
                            <div className="text-xs text-slate-500 mb-1">Detection</div>
                            <div className={`text-xs font-bold px-2 py-1 rounded-full ${
                              item.detection_result === 'FOUND' ? 'bg-red-500/20 text-red-400' : 'bg-emerald-500/20 text-emerald-400'
                            }`}>
                              {item.detection_result === 'FOUND' ? '🧬 DETECTED' : '✅ NOT FOUND'}
                            </div>
                            {item.confidence != null && (
                              <div className="text-xs text-slate-500 mt-0.5 text-right">{(item.confidence * 100).toFixed(1)}% conf</div>
                            )}
                          </div>
                        )}
                        {isBio && item.marker_variant && (
                          <div className="text-right hidden sm:block">
                            <div className="text-xs text-slate-500 mb-1">Variant</div>
                            <div className="font-mono text-amber-400 text-xs">{item.marker_variant}</div>
                          </div>
                        )}
                        {!isBio && (
                          <div className="text-right">
                            <div className="text-xs text-slate-500 mb-1">Target Bits</div>
                            <div className="font-mono text-white bg-slate-950 px-2 py-1 rounded inline-block tracking-widest text-sm">
                              {item.target_bits || 'N/A'}
                            </div>
                          </div>
                        )}
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
                            <div className="font-mono text-slate-400 text-xs w-24 truncate" title={item.job_id}>{item.job_id}</div>
                          </div>
                        )}
                        <div className={`ml-2 transition-transform duration-200 ${isExpanded ? 'rotate-90' : ''} text-slate-500 group-hover:text-slate-300`}>
                          <ChevronRight className="w-4 h-4" />
                        </div>
                      </div>
                    </div>
                  </button>

                  {/* Expandable panel */}
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
                          {isBio ? (
                            <SingleBioPanel item={item} groupId={id} />
                          ) : (<>
                            {loadingCircuit === id && (
                              <div className="flex items-center gap-3 justify-center py-8 text-slate-400 text-sm">
                                <div className="w-5 h-5 border-2 border-amber-500/30 border-t-amber-500 rounded-full animate-spin" />
                                Building circuit visualization…
                              </div>
                            )}
                            {!bits && !loadingCircuit && (
                              <p className="text-slate-500 text-sm text-center py-6">
                                No target bits recorded — circuit cannot be reconstructed.
                              </p>
                            )}
                            {circuitError[bits] && !loadingCircuit && !circuit && (
                              <p className="text-red-400 text-sm text-center py-6">
                                Failed to load circuit visualization. Please try again.
                              </p>
                            )}
                            {circuit && (<>
                              <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 text-center text-sm">
                                {[
                                  { label: 'Qubits',     value: circuit.qubits,    color: 'text-blue-400' },
                                  { label: 'Iterations', value: circuit.iterations, color: 'text-purple-400' },
                                  { label: 'Noise',      value: item.noise_level != null ? `${(item.noise_level * 100).toFixed(1)}%` : '—', color: 'text-amber-400' },
                                  { label: 'Time (ms)',  value: item.execution_time_ms != null ? item.execution_time_ms.toFixed(1) : '—', color: 'text-emerald-400' },
                                ].map((s) => (
                                  <div key={s.label} className="bg-slate-900 border border-slate-800 rounded-xl p-3">
                                    <div className="text-slate-500 text-xs uppercase tracking-wider mb-1">{s.label}</div>
                                    <div className={`font-mono font-bold ${s.color}`}>{s.value}</div>
                                  </div>
                                ))}
                              </div>
                              <GroverStepNavigator
                                stepCircuits={circuit.step_circuits}
                                activeStep={stepIdx}
                                onStepChange={(idx) => setStep(id, idx)}
                                targetBits={bits}
                              />
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
                            </>)}
                          </>)}
                        </div>
                      </motion.div>
                    )}
                  </AnimatePresence>
                </motion.div>
              );
            })
          )}
        </AnimatePresence>

        {/* Infinite-scroll sentinel */}
        <div ref={sentinelRef} className="h-4" />
        {loadingMore && (
          <div className="flex justify-center py-4">
            <Loader2 className="w-6 h-6 text-slate-500 animate-spin" />
          </div>
        )}
        {!hasMore && displayGroups.length > 0 && !loading && (
          <p className="text-center text-xs text-slate-600 py-2">
            — end of history —
          </p>
        )}
      </div>
    </div>
  );
}