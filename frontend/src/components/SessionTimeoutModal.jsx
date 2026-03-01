import React, { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Clock, LogOut, RefreshCw } from 'lucide-react';

/**
 * SessionTimeoutModal
 *
 * Untitled-UI-inspired warning modal shown before the session expires.
 * Features a circular SVG countdown ring, glassmorphism card,
 * and extend/logout action buttons.
 *
 * Props:
 *   visible       – boolean, whether the modal is shown
 *   secondsLeft   – number, countdown seconds remaining
 *   totalWarning  – number, total warning duration in seconds (for progress calc)
 *   onExtend      – callback to reset session timer
 *   onLogout      – callback to log out immediately
 */
export default function SessionTimeoutModal({ visible, secondsLeft, totalWarning, onExtend, onLogout }) {
    // SVG circle geometry
    const RADIUS = 54;
    const CIRCUMFERENCE = 2 * Math.PI * RADIUS;
    const progress = secondsLeft / totalWarning;          // 1 → 0
    const strokeOffset = CIRCUMFERENCE * (1 - progress);  // 0 → full

    // Format mm:ss
    const mins = Math.floor(secondsLeft / 60);
    const secs = secondsLeft % 60;
    const timeStr = `${mins}:${secs.toString().padStart(2, '0')}`;

    return (
        <AnimatePresence>
            {visible && (
                <motion.div
                    key="timeout-overlay"
                    initial={{ opacity: 0 }}
                    animate={{ opacity: 1 }}
                    exit={{ opacity: 0 }}
                    transition={{ duration: 0.25 }}
                    className="fixed inset-0 z-[9999] flex items-center justify-center"
                >
                    {/* Backdrop */}
                    <div className="absolute inset-0 bg-slate-950/70 backdrop-blur-sm" />

                    {/* Modal card */}
                    <motion.div
                        initial={{ opacity: 0, scale: 0.92, y: 16 }}
                        animate={{ opacity: 1, scale: 1, y: 0 }}
                        exit={{ opacity: 0, scale: 0.92, y: 16 }}
                        transition={{ type: 'spring', damping: 24, stiffness: 300 }}
                        className="relative z-10 w-full max-w-sm mx-4"
                    >
                        <div className="bg-slate-900/80 backdrop-blur-xl border border-slate-700/50 rounded-2xl shadow-2xl shadow-black/40 p-8">
                            {/* Circular countdown */}
                            <div className="flex justify-center mb-6">
                                <div className="relative w-32 h-32">
                                    {/* Background ring */}
                                    <svg className="w-full h-full -rotate-90" viewBox="0 0 120 120">
                                        <circle
                                            cx="60" cy="60" r={RADIUS}
                                            fill="none"
                                            stroke="currentColor"
                                            strokeWidth="6"
                                            className="text-slate-800"
                                        />
                                        {/* Progress ring */}
                                        <circle
                                            cx="60" cy="60" r={RADIUS}
                                            fill="none"
                                            strokeWidth="6"
                                            strokeLinecap="round"
                                            strokeDasharray={CIRCUMFERENCE}
                                            strokeDashoffset={strokeOffset}
                                            className={`transition-[stroke-dashoffset] duration-1000 ease-linear ${secondsLeft <= 30 ? 'text-red-500' : secondsLeft <= 60 ? 'text-amber-500' : 'text-blue-500'
                                                }`}
                                            stroke="currentColor"
                                        />
                                    </svg>
                                    {/* Center text */}
                                    <div className="absolute inset-0 flex flex-col items-center justify-center">
                                        <span className={`text-2xl font-bold tabular-nums ${secondsLeft <= 30 ? 'text-red-400' : secondsLeft <= 60 ? 'text-amber-400' : 'text-white'
                                            }`}>
                                            {timeStr}
                                        </span>
                                        <span className="text-[10px] uppercase tracking-widest text-slate-500 mt-0.5">remaining</span>
                                    </div>
                                </div>
                            </div>

                            {/* Icon + heading */}
                            <div className="text-center mb-6">
                                <div className="w-10 h-10 bg-amber-500/10 border border-amber-500/20 rounded-xl flex items-center justify-center mx-auto mb-3">
                                    <Clock className="w-5 h-5 text-amber-400" />
                                </div>
                                <h2 className="text-lg font-semibold text-white mb-1.5">Session Expiring Soon</h2>
                                <p className="text-sm text-slate-400 leading-relaxed">
                                    Your session will expire due to inactivity. Would you like to continue?
                                </p>
                            </div>

                            {/* Actions */}
                            <div className="space-y-2.5">
                                <button
                                    onClick={onExtend}
                                    className="w-full bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-500 hover:to-purple-500 text-white font-medium py-2.5 rounded-xl flex items-center justify-center gap-2 transition-all text-sm shadow-lg shadow-blue-500/20 hover:shadow-blue-500/30"
                                >
                                    <RefreshCw className="w-4 h-4" /> Extend Session
                                </button>
                                <button
                                    onClick={onLogout}
                                    className="w-full border border-slate-700 hover:border-slate-600 text-slate-300 hover:text-white font-medium py-2.5 rounded-xl flex items-center justify-center gap-2 transition-all text-sm bg-slate-800/40 hover:bg-slate-800/70"
                                >
                                    <LogOut className="w-4 h-4" /> Log Out
                                </button>
                            </div>
                        </div>
                    </motion.div>
                </motion.div>
            )}
        </AnimatePresence>
    );
}
