import React, { useState } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import axios from 'axios';
import { Dna, Lock, CheckCircle, XCircle } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000/api';

const PASSWORD_RULES = [
  { label: 'At least 8 characters',             test: (p) => p.length >= 8 },
  { label: 'One uppercase letter',               test: (p) => /[A-Z]/.test(p) },
  { label: 'One number or special character',    test: (p) => /[0-9!@#$%^&*]/.test(p) },
];

export default function ResetPassword() {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const token = searchParams.get('token') || '';

  const [password, setPassword]       = useState('');
  const [confirm, setConfirm]         = useState('');
  const [touched, setTouched]         = useState({ password: false, confirm: false });
  const [loading, setLoading]         = useState(false);
  const [success, setSuccess]         = useState(false);
  const [serverError, setServerError] = useState('');

  const passwordValid  = PASSWORD_RULES.every((r) => r.test(password));
  const confirmValid   = confirm === password && confirm.length > 0;
  const formValid      = passwordValid && confirmValid;

  const handleSubmit = async (e) => {
    e.preventDefault();
    setTouched({ password: true, confirm: true });
    if (!formValid) return;

    setLoading(true);
    setServerError('');
    try {
      await axios.post(`${API_BASE}/auth/reset-password`, {
        token,
        new_password: password,
      });
      setSuccess(true);
    } catch (err) {
      setServerError(
        err.response?.data?.detail || 'Something went wrong. Please try again.'
      );
    } finally {
      setLoading(false);
    }
  };

  // ── Success screen ────────────────────────────────────────────────────────────
  if (success) {
    return (
      <div className="min-h-screen flex flex-col items-center justify-center p-6 bg-slate-950 relative overflow-hidden">
        <div className="absolute top-0 left-0 w-[500px] h-[500px] bg-blue-600/10 rounded-full blur-[120px] -translate-x-1/2 -translate-y-1/2 pointer-events-none" />
        <div className="absolute bottom-0 right-0 w-[500px] h-[500px] bg-purple-600/10 rounded-full blur-[120px] translate-x-1/2 translate-y-1/2 pointer-events-none" />
        <div className="relative z-10 w-full max-w-md">
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            className="bg-slate-900/60 border border-slate-800 p-8 rounded-3xl shadow-2xl backdrop-blur-xl text-center"
          >
            <div className="w-16 h-16 bg-gradient-to-br from-green-500 to-emerald-600 rounded-2xl flex items-center justify-center mx-auto mb-5 shadow-lg shadow-green-500/20">
              <CheckCircle className="w-8 h-8 text-white" />
            </div>
            <h2 className="text-2xl font-bold text-white mb-2">Password updated!</h2>
            <p className="text-slate-400 text-sm mb-6">
              Your password has been reset successfully. You can now sign in with your new password.
            </p>
            <motion.button
              whileHover={{ scale: 1.02 }}
              whileTap={{ scale: 0.98 }}
              onClick={() => navigate('/login')}
              className="w-full bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-500 hover:to-purple-500 text-white font-medium py-3 rounded-xl transition-all"
            >
              Sign In
            </motion.button>
          </motion.div>
        </div>
      </div>
    );
  }

  // ── Invalid / missing token ───────────────────────────────────────────────────
  if (!token) {
    return (
      <div className="min-h-screen flex flex-col items-center justify-center p-6 bg-slate-950">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-slate-900/60 border border-slate-800 p-8 rounded-3xl shadow-2xl backdrop-blur-xl text-center max-w-md w-full"
        >
          <div className="w-16 h-16 bg-gradient-to-br from-red-500 to-red-700 rounded-2xl flex items-center justify-center mx-auto mb-5">
            <XCircle className="w-8 h-8 text-white" />
          </div>
          <h2 className="text-xl font-bold text-white mb-2">Invalid Reset Link</h2>
          <p className="text-slate-400 text-sm mb-6">
            This link is missing a reset token. Please request a new password reset.
          </p>
          <button
            onClick={() => navigate('/login')}
            className="text-blue-400 hover:text-blue-300 text-sm transition-colors"
          >
            ← Back to Sign In
          </button>
        </motion.div>
      </div>
    );
  }

  // ── Reset form ────────────────────────────────────────────────────────────────
  return (
    <div className="min-h-screen flex flex-col items-center justify-center p-6 relative overflow-hidden bg-slate-950">
      <div className="absolute top-0 left-0 w-[500px] h-[500px] bg-blue-600/10 rounded-full blur-[120px] -translate-x-1/2 -translate-y-1/2 pointer-events-none" />
      <div className="absolute bottom-0 right-0 w-[500px] h-[500px] bg-purple-600/10 rounded-full blur-[120px] translate-x-1/2 translate-y-1/2 pointer-events-none" />

      <div className="relative z-10 w-full max-w-md">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-slate-900/60 border border-slate-800 p-8 rounded-3xl shadow-2xl backdrop-blur-xl"
        >
          {/* Logo */}
          <div className="flex justify-center mb-6">
            <div className="w-20 h-20 bg-gradient-to-br from-blue-500 to-purple-600 rounded-2xl flex items-center justify-center shadow-lg shadow-blue-500/20">
              <Dna className="w-10 h-10 text-white" />
            </div>
          </div>

          <h1 className="text-2xl font-bold text-center text-white mb-2">Set New Password</h1>
          <p className="text-slate-400 text-center text-sm mb-8">
            Choose a strong password for your BioQuantum account.
          </p>

          <form onSubmit={handleSubmit} noValidate className="space-y-5">
            {/* New password */}
            <div className="space-y-1">
              <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">
                New Password
              </label>
              <div className="relative">
                <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-600" />
                <input
                  id="reset-new-password"
                  type="password"
                  value={password}
                  onChange={(e) => { setPassword(e.target.value); setServerError(''); }}
                  onBlur={() => setTouched((t) => ({ ...t, password: true }))}
                  className={`w-full bg-slate-950/50 border rounded-xl pl-10 pr-4 py-3 text-white placeholder-slate-600 outline-none transition-colors
                    ${touched.password && !passwordValid ? 'border-red-500 focus:border-red-400' : 'border-slate-800 focus:border-purple-500'}`}
                  placeholder="Enter new password"
                  autoComplete="new-password"
                />
              </div>
              {/* Password strength checklist */}
              <ul className="mt-2 space-y-1 pl-1">
                {PASSWORD_RULES.map(({ label, test }) => (
                  <li
                    key={label}
                    className={`text-xs flex items-center gap-1.5 transition-colors ${test(password) ? 'text-green-400' : 'text-slate-500'}`}
                  >
                    <span>{test(password) ? '✓' : '○'}</span>{label}
                  </li>
                ))}
              </ul>
            </div>

            {/* Confirm password */}
            <div className="space-y-1">
              <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">
                Confirm Password
              </label>
              <div className="relative">
                <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-600" />
                <input
                  id="reset-confirm-password"
                  type="password"
                  value={confirm}
                  onChange={(e) => { setConfirm(e.target.value); setServerError(''); }}
                  onBlur={() => setTouched((t) => ({ ...t, confirm: true }))}
                  className={`w-full bg-slate-950/50 border rounded-xl pl-10 pr-4 py-3 text-white placeholder-slate-600 outline-none transition-colors
                    ${touched.confirm && !confirmValid ? 'border-red-500 focus:border-red-400' : 'border-slate-800 focus:border-purple-500'}`}
                  placeholder="Confirm new password"
                  autoComplete="new-password"
                />
              </div>
              <AnimatePresence>
                {touched.confirm && !confirmValid && (
                  <motion.p
                    initial={{ opacity: 0, y: -4 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -4 }}
                    className="text-red-400 text-xs mt-1 ml-1"
                  >
                    Passwords do not match.
                  </motion.p>
                )}
              </AnimatePresence>
            </div>

            {/* Server error */}
            {serverError && (
              <div className="text-red-400 text-sm text-center bg-red-900/20 p-3 rounded-xl border border-red-900/50">
                {serverError}
              </div>
            )}

            <motion.button
              whileHover={{ scale: 1.02 }}
              whileTap={{ scale: 0.98 }}
              disabled={loading}
              type="submit"
              className="w-full bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-500 hover:to-purple-500 text-white font-medium py-3 rounded-xl flex items-center justify-center gap-2 transition-all mt-2 disabled:opacity-50"
            >
              {loading ? (
                <div className="w-5 h-5 border-2 border-white/30 border-t-white rounded-full animate-spin" />
              ) : (
                'Reset Password'
              )}
            </motion.button>
          </form>

          <div className="text-center mt-6">
            <button
              onClick={() => navigate('/login')}
              className="text-slate-500 hover:text-slate-300 text-sm transition-colors"
            >
              ← Back to Sign In
            </button>
          </div>
        </motion.div>
      </div>
    </div>
  );
}
