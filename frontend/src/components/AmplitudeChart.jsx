/**
 * AmplitudeChart — Dynamic bar chart showing the probability distribution
 * of all basis states at the current Grover step.
 *
 * Visually demonstrates amplitude amplification: the target state's bar
 * grows dramatically from uniform (Step 1) → inverted phase (Step 2)
 * → amplified (Step 3) → measured (Step 4).
 */
import React, { useMemo } from 'react';
import {
    BarChart, Bar, XAxis, YAxis, CartesianGrid,
    Tooltip, ResponsiveContainer, Cell, ReferenceLine
} from 'recharts';
import { motion } from 'framer-motion';

const STEP_COLORS = ['#3b82f6', '#a855f7', '#f59e0b', '#10b981'];

export default function AmplitudeChart({ probabilities, targetBits, stepIndex = 0 }) {
    const data = useMemo(() => {
        if (!probabilities) return [];
        return Object.entries(probabilities)
            .sort(([a], [b]) => a.localeCompare(b))
            .map(([state, prob]) => ({
                state: `|${state}⟩`,
                rawState: state,
                probability: prob,
                isTarget: state === targetBits,
            }));
    }, [probabilities, targetBits]);

    if (!data.length) return null;

    const maxProb = Math.max(...data.map(d => d.probability));
    const avgProb = 1 / data.length;
    const barColor = STEP_COLORS[stepIndex % STEP_COLORS.length];

    return (
        <motion.div
            initial={{ opacity: 0, y: 8 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.3 }}
            className="w-full"
        >
            <div className="flex items-center justify-between mb-3">
                <div className="text-[10px] text-slate-500 uppercase tracking-widest font-medium">
                    Probability Distribution
                </div>
                <div className="flex items-center gap-3 text-[10px]">
                    <span className="flex items-center gap-1.5">
                        <span className="w-2.5 h-2.5 rounded-sm bg-amber-400 inline-block" />
                        <span className="text-slate-400">Target</span>
                    </span>
                    <span className="flex items-center gap-1.5">
                        <span className="w-2.5 h-2.5 rounded-sm inline-block" style={{ backgroundColor: barColor, opacity: 0.6 }} />
                        <span className="text-slate-400">Other states</span>
                    </span>
                </div>
            </div>
            <ResponsiveContainer width="100%" height={180}>
                <BarChart data={data} margin={{ top: 5, right: 10, left: -15, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" vertical={false} />
                    <XAxis
                        dataKey="state"
                        stroke="#475569"
                        tick={{ fontSize: 9, fill: '#64748b' }}
                        axisLine={false}
                        tickLine={false}
                        interval={data.length > 16 ? Math.floor(data.length / 8) : 0}
                    />
                    <YAxis
                        stroke="#475569"
                        tick={{ fontSize: 9, fill: '#64748b' }}
                        axisLine={false}
                        tickLine={false}
                        domain={[0, Math.min(1, maxProb * 1.15)]}
                        tickFormatter={(v) => `${(v * 100).toFixed(0)}%`}
                    />
                    <Tooltip
                        contentStyle={{
                            backgroundColor: '#0f172a',
                            borderColor: '#334155',
                            color: '#ffffff',
                            borderRadius: '10px',
                            fontSize: '12px',
                            boxShadow: '0 10px 15px -3px rgb(0 0 0 / 0.5)',
                        }}
                        itemStyle={{ color: '#ffffff' }}
                        labelStyle={{ color: '#ffffff' }}
                        formatter={(value, name, props) => [
                            `${(value * 100).toFixed(2)}%`,
                            props.payload.isTarget ? '🎯 Target' : 'Probability',
                        ]}
                        cursor={{ fill: 'rgba(148, 163, 184, 0.08)' }}
                    />
                    <ReferenceLine
                        y={avgProb}
                        stroke="#64748b"
                        strokeDasharray="4 4"
                        strokeWidth={1}
                        label={{
                            value: `avg ${(avgProb * 100).toFixed(1)}%`,
                            position: 'right',
                            fill: '#64748b',
                            fontSize: 9,
                        }}
                    />
                    <Bar dataKey="probability" radius={[3, 3, 0, 0]} animationDuration={600}>
                        {data.map((entry, idx) => (
                            <Cell
                                key={idx}
                                fill={entry.isTarget ? '#f59e0b' : barColor}
                                fillOpacity={entry.isTarget ? 1 : 0.5}
                                stroke={entry.isTarget ? '#fbbf24' : 'none'}
                                strokeWidth={entry.isTarget ? 1.5 : 0}
                            />
                        ))}
                    </Bar>
                </BarChart>
            </ResponsiveContainer>
        </motion.div>
    );
}
