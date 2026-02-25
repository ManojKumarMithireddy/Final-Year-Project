import React from 'react';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid,
  Tooltip, Legend, ResponsiveContainer
} from 'recharts';
import { motion } from 'framer-motion';

/**
 * Generates classical vs. quantum complexity data points for the chart.
 * NOTE: These values are analytically computed — not measured from a real quantum circuit.
 */
function generateChartData(classicalResult) {
  const maxN = classicalResult.n_windows;
  const step = Math.max(1, Math.floor(maxN / 10));
  const data = [];
  for (let i = step; i <= maxN; i += step) {
    data.push({
      windows: i,
      classical: i,
      quantum: Math.floor((Math.PI / 4) * Math.sqrt(i))
    });
  }
  if (data.length === 0 || data[data.length - 1].windows !== maxN) {
    data.push({
      windows: maxN,
      classical: maxN,
      quantum: Math.floor((Math.PI / 4) * Math.sqrt(maxN))
    });
  }
  return data;
}

export default function ComplexityChart({ classicalResult, quantumMetrics }) {
  if (!classicalResult || !quantumMetrics) {
    return (
      <div className="flex-1 flex items-center justify-center text-slate-600/70 text-sm">
        Chart will appear after running a search
      </div>
    );
  }

  const data = generateChartData(classicalResult);

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      className="flex-1 w-full relative"
    >
      <ResponsiveContainer width="100%" height="100%">
        <LineChart data={data} margin={{ top: 15, right: 30, left: 10, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" vertical={false} />
          <XAxis
            dataKey="windows"
            stroke="#64748b"
            tickFormatter={(val) => `${(val / 1000).toFixed(0)}k`}
            tick={{ fontSize: 12 }}
            axisLine={false}
            tickLine={false}
            dy={10}
          />
          <YAxis
            stroke="#64748b"
            tickFormatter={(val) => `${(val / 1000).toFixed(0)}k`}
            tick={{ fontSize: 12 }}
            axisLine={false}
            tickLine={false}
            dx={-10}
          />
          <Tooltip
            contentStyle={{
              backgroundColor: '#0f172a',
              borderColor: '#334155',
              color: '#f8fafc',
              borderRadius: '12px',
              boxShadow: '0 10px 15px -3px rgb(0 0 0 / 0.5)'
            }}
            itemStyle={{ color: '#e2e8f0' }}
            cursor={{ stroke: '#334155', strokeWidth: 1, strokeDasharray: '5 5' }}
          />
          <Legend verticalAlign="top" height={40} iconType="circle" />
          <Line
            type="monotone"
            dataKey="classical"
            name="Classical O(N)"
            stroke="#3b82f6"
            strokeWidth={3}
            dot={false}
            activeDot={{ r: 6, fill: '#3b82f6', stroke: '#0f172a', strokeWidth: 2 }}
          />
          <Line
            type="monotone"
            dataKey="quantum"
            name="Quantum O(√N) — Theoretical"
            stroke="#a855f7"
            strokeWidth={3}
            dot={false}
            activeDot={{ r: 6, fill: '#a855f7', stroke: '#0f172a', strokeWidth: 2 }}
          />
        </LineChart>
      </ResponsiveContainer>
    </motion.div>
  );
}
