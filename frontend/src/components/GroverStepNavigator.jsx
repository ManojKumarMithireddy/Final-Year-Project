/**
 * GroverStepNavigator
 * ─────────────────────────────────────────────────────────────────────────────
 * Shared step-by-step Grover circuit visualiser used by both the Quantum Sim
 * POC panel and the History page.
 *
 * Props
 *   stepCircuits  – array returned by the /circuit-info or /quantum-simulation-poc
 *                   endpoint (step.label, step.description, step.diagram,
 *                   step.gates, step.probabilities)
 *   activeStep    – currently selected step index (controlled)
 *   onStepChange  – (idx: number) => void
 *   targetBits    – bitstring used as the search target (for AmplitudeChart)
 */

import React from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import AmplitudeChart from './AmplitudeChart';

const STEP_COLORS = [
  {
    active: 'text-blue-400 border-blue-500 bg-blue-500/10',
    inactive: 'text-slate-400 hover:text-blue-400',
    badge: 'bg-blue-500/20 text-blue-400 border-blue-500/30',
    card: 'border-blue-500/20 bg-blue-950/20',
  },
  {
    active: 'text-purple-400 border-purple-500 bg-purple-500/10',
    inactive: 'text-slate-400 hover:text-purple-400',
    badge: 'bg-purple-500/20 text-purple-400 border-purple-500/30',
    card: 'border-purple-500/20 bg-purple-950/20',
  },
  {
    active: 'text-amber-400 border-amber-500 bg-amber-500/10',
    inactive: 'text-slate-400 hover:text-amber-400',
    badge: 'bg-amber-500/20 text-amber-400 border-amber-500/30',
    card: 'border-amber-500/20 bg-amber-950/20',
  },
  {
    active: 'text-emerald-400 border-emerald-500 bg-emerald-500/10',
    inactive: 'text-slate-400 hover:text-emerald-400',
    badge: 'bg-emerald-500/20 text-emerald-400 border-emerald-500/30',
    card: 'border-emerald-500/20 bg-emerald-950/20',
  },
];

export default function GroverStepNavigator({
  stepCircuits,
  activeStep,
  onStepChange,
  targetBits,
}) {
  if (!stepCircuits || stepCircuits.length === 0) return null;

  const currentStep = stepCircuits[activeStep];

  return (
    <div className="border border-slate-800 rounded-xl overflow-hidden bg-slate-950">
      {/* ── Tab bar ── */}
      <div className="flex border-b border-slate-800 bg-slate-900/60 overflow-x-auto">
        {stepCircuits.map((step, idx) => {
          const c = STEP_COLORS[idx % STEP_COLORS.length];
          const isActive = activeStep === idx;
          return (
            <button
              key={idx}
              onClick={() => onStepChange(idx)}
              className={`flex-1 min-w-[120px] px-4 py-3 text-xs font-semibold border-b-2 transition-all flex flex-col items-center gap-1 ${
                isActive
                  ? `${c.active} border-current`
                  : `border-transparent ${c.inactive}`
              }`}>
              <span className="opacity-60 text-[10px] uppercase tracking-wider">
                Step {step.step}
              </span>
              <span>{step.label}</span>
            </button>
          );
        })}
      </div>

      {/* ── Active step content ── */}
      <AnimatePresence mode="wait">
        {currentStep && (
          <motion.div
            key={activeStep}
            initial={{ opacity: 0, y: 6 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -6 }}
            transition={{ duration: 0.18 }}>
            {/* Description banner */}
            <div className="px-5 pt-4 pb-3 text-xs text-slate-400 leading-relaxed border-b border-slate-800/60">
              {currentStep.description}
            </div>

            {/* Composed circuit diagram — all parts inline, active one highlighted */}
            <div className="flex flex-col">
              <div className="p-5 overflow-x-auto border-b border-slate-800/60 pb-8">
                <div className="text-[10px] text-slate-500 uppercase tracking-widest mb-4 font-medium flex items-center justify-between">
                  <span>Composed Circuit Diagram</span>
                  <span className="text-amber-500/80 normal-case font-normal bg-amber-500/10 px-2 py-0.5 rounded-md border border-amber-500/20">
                    Active step highlighted
                  </span>
                </div>
                <div className="bg-slate-950/80 p-4 rounded-xl border border-slate-800 shadow-inner overflow-x-auto">
                  <div className="flex flex-row items-start pointer-events-none w-max">
                    {stepCircuits.map((s, idx) => (
                      <pre
                        key={`diagram-part-${idx}`}
                        className={`text-[10px] leading-[11px] m-0 p-0 transition-all duration-300 ${
                          idx === activeStep
                            ? 'text-amber-400 drop-shadow-[0_0_12px_rgba(251,191,36,0.8)]'
                            : 'text-slate-300 opacity-30'
                        }`}>
                        {s.diagram}
                      </pre>
                    ))}
                  </div>
                </div>
              </div>

              {/* Gate explanation cards */}
              {currentStep.gates && currentStep.gates.length > 0 && (
                <div className="p-5 space-y-4 bg-slate-900/40">
                  <div className="text-[10px] text-slate-500 uppercase tracking-widest font-medium mb-2">
                    Gates Used in This Step
                  </div>
                  <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                    {currentStep.gates.map((gate, gi) => {
                      const gc = STEP_COLORS[gi % STEP_COLORS.length];
                      return (
                        <div
                          key={gi}
                          className={`rounded-xl border p-3 ${gc.card}`}>
                          <div className="flex items-center gap-2 mb-2">
                            <span
                              className={`text-xs font-mono font-bold px-2 py-0.5 rounded-md border ${gc.badge}`}>
                              {gate.symbol}
                            </span>
                            <div className="flex-1 min-w-0">
                              <div className="text-xs font-semibold text-slate-200 truncate">
                                {gate.name}
                              </div>
                              <div className="text-[10px] text-slate-500">
                                ×{gate.count} application
                                {gate.count !== 1 ? 's' : ''}
                              </div>
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

            {/* Amplitude chart */}
            {currentStep.probabilities && (
              <div className="px-5 py-4 border-t border-slate-800/60">
                <AmplitudeChart
                  probabilities={currentStep.probabilities}
                  targetBits={targetBits}
                  stepIndex={activeStep}
                />
              </div>
            )}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
