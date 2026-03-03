/**
 * Shows a sticky banner when the HF Spaces backend is offline / waking.
 * Disappears automatically once the backend responds (with a 3 s "online" flash).
 */
import React, { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { WifiOff, RefreshCw, CheckCircle } from 'lucide-react';

export default function BackendSleepBanner({ status, onRetry }) {
  const [showOnline, setShowOnline] = useState(false);

  // Flash "online" confirmation for 3 s after waking up
  useEffect(() => {
    if (status === 'online') {
      setShowOnline(true);
      const t = setTimeout(() => setShowOnline(false), 3000);
      return () => clearTimeout(t);
    }
  }, [status]);

  const sleeping = status === 'offline' || status === 'waking' || status === 'checking';

  return (
    <AnimatePresence>
      {sleeping && (
        <motion.div
          key="sleep-banner"
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          transition={{ duration: 0.3 }}
          className="sticky top-16 z-40 w-full"
        >
          <div className={`
            flex items-center justify-between gap-4 px-4 py-2.5 text-sm
            border-b
            ${status === 'waking'
              ? 'bg-amber-950/80 border-amber-800/60 text-amber-200 backdrop-blur-md'
              : status === 'checking'
              ? 'bg-slate-900/80 border-slate-700/60 text-slate-300 backdrop-blur-md'
              : 'bg-red-950/80 border-red-900/60 text-red-200 backdrop-blur-md'}
          `}>
            <div className="flex items-center gap-3 min-w-0">
              {status === 'waking' || status === 'checking' ? (
                <RefreshCw className="w-4 h-4 shrink-0 animate-spin" />
              ) : (
                <WifiOff className="w-4 h-4 shrink-0" />
              )}
              <span className="truncate">
                {status === 'checking' && 'Checking backend connection…'}
                {status === 'offline'  && (
                  <><strong>Backend is sleeping</strong>{' '}— Hugging Face Spaces pauses free-tier apps after inactivity.</>
                )}
                {status === 'waking' && (
                  <><strong>Waking up backend…</strong>{' '}This usually takes 30–60 s. Hang tight!</>
                )}
              </span>
            </div>

            {status === 'offline' && (
              <button
                onClick={onRetry}
                className="shrink-0 flex items-center gap-1.5 bg-red-800/60 hover:bg-red-700/60 border border-red-700/50 px-3 py-1 rounded-lg transition-colors text-xs font-medium"
              >
                <RefreshCw className="w-3.5 h-3.5" /> Wake up
              </button>
            )}

            {status === 'waking' && (
              <span className="shrink-0 flex gap-1">
                {[0, 1, 2].map((i) => (
                  <motion.span
                    key={i}
                    className="w-1.5 h-1.5 rounded-full bg-amber-400"
                    animate={{ opacity: [0.3, 1, 0.3] }}
                    transition={{ duration: 1.2, repeat: Infinity, delay: i * 0.4 }}
                  />
                ))}
              </span>
            )}
          </div>
        </motion.div>
      )}

      {showOnline && (
        <motion.div
          key="online-flash"
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          transition={{ duration: 0.3 }}
          className="sticky top-16 z-40 w-full"
        >
          <div className="flex items-center gap-3 px-4 py-2.5 text-sm bg-emerald-950/80 border-b border-emerald-800/60 text-emerald-200 backdrop-blur-md">
            <CheckCircle className="w-4 h-4 shrink-0" />
            <span><strong>Backend is online</strong> — ready to search.</span>
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
