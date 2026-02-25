import React, { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import HybridSearchPanel from '../components/HybridSearchPanel';
import GroverPOC from '../components/GroverPOC';

const TABS = [
  { id: 'search', label: 'Large-Scale Hybrid Search' },
  { id: 'poc',    label: 'Quantum Sim Proof of Concept' },
];

export default function Dashboard() {
  const [activeTab, setActiveTab] = useState('search');

  return (
    <div className="max-w-7xl mx-auto p-6 md:p-8">
      {/* Tab Bar */}
      <div className="flex justify-center mb-8">
        <div className="flex gap-2 bg-slate-900/50 p-1.5 rounded-2xl border border-slate-800 backdrop-blur-md relative">
          {TABS.map(({ id, label }) => (
            <button
              key={id}
              onClick={() => setActiveTab(id)}
              className={`relative px-6 py-2.5 rounded-xl font-medium text-sm transition-colors z-10 ${
                activeTab === id ? 'text-white' : 'text-slate-400 hover:text-white'
              }`}
            >
              {activeTab === id && (
                <motion.div
                  layoutId="activeTabBadge"
                  className="absolute inset-0 bg-blue-600 rounded-xl"
                  initial={false}
                  transition={{ type: 'spring', stiffness: 400, damping: 30 }}
                />
              )}
              <span className="relative z-10">{label}</span>
            </button>
          ))}
        </div>
      </div>

      {/* Tab Content */}
      <AnimatePresence mode="wait">
        <motion.div
          key={activeTab}
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          transition={{ duration: 0.2 }}
        >
          {activeTab === 'search' && <HybridSearchPanel />}
          {activeTab === 'poc'    && <GroverPOC />}
        </motion.div>
      </AnimatePresence>
    </div>
  );
}
