import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { GoogleLogin } from '@react-oauth/google';
import axios from 'axios';
import { Dna, Lock, Mail, ArrowRight, MailCheck, RefreshCw, KeyRound } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000/api';

// ── Validation helpers ────────────────────────────────────────────────────────
const EMAIL_REGEX = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;

function validateEmail(value) {
  if (!value.trim()) return 'Email is required.';
  if (!EMAIL_REGEX.test(value)) return 'Enter a valid email address.';
  return '';
}

function validatePassword(value, isLogin) {
  if (!value) return 'Password is required.';
  if (!isLogin) {
    if (value.length < 8) return 'Password must be at least 8 characters.';
    if (!/[A-Z]/.test(value)) return 'Include at least one uppercase letter.';
    if (!/[0-9!@#$%^&*]/.test(value))
      return 'Include at least one number or special character (!@#$%^&*).';
  }
  return '';
}

function validateName(value) {
  if (value.length > 60) return 'Name must be 60 characters or fewer.';
  return '';
}

// ── Component ─────────────────────────────────────────────────────────────────
export default function Login() {
  const navigate = useNavigate();
  const [serverError, setServerError] = useState('');
  const [isLogin, setIsLogin] = useState(true);
  const [formData, setFormData] = useState({ email: '', password: '', name: '' });
  const [touched, setTouched] = useState({ email: false, password: false, name: false });
  const [loading, setLoading] = useState(false);

  // "Check your inbox" state — shown after successful registration
  const [pendingEmail, setPendingEmail] = useState('');
  const [resendLoading, setResendLoading] = useState(false);
  const [resendSent, setResendSent] = useState(false);

  // "not verified" state — shown when login returns 403
  const [notVerifiedEmail, setNotVerifiedEmail] = useState('');

  // Forgot password state
  const [forgotMode, setForgotMode]       = useState(false);
  const [forgotEmail, setForgotEmail]     = useState('');
  const [forgotSent, setForgotSent]       = useState(false);
  const [forgotLoading, setForgotLoading] = useState(false);
  const [forgotError, setForgotError]     = useState('');

  // Derived per-field errors
  const fieldErrors = {
    name: touched.name ? validateName(formData.name) : '',
    email: touched.email ? validateEmail(formData.email) : '',
    password: touched.password ? validatePassword(formData.password, isLogin) : '',
  };

  const isFormValid =
    !validateEmail(formData.email) &&
    !validatePassword(formData.password, isLogin) &&
    (!isLogin ? !validateName(formData.name) : true);

  const handleBlur = (field) => setTouched((prev) => ({ ...prev, [field]: true }));

  const handleChange = (field, value) => {
    setFormData((prev) => ({ ...prev, [field]: value }));
    setServerError('');
    setNotVerifiedEmail('');
  };

  const handleGoogleSuccess = async (credentialResponse) => {
    try {
      const res = await axios.post(`${API_BASE}/auth/google`, {
        token: credentialResponse.credential,
      });
      localStorage.setItem('token', res.data.access_token);
      localStorage.setItem('userName', res.data.user.name);
      navigate('/');
    } catch (err) {
      setServerError(
        err.response?.data?.detail
          ? `Backend Error: ${err.response.data.detail}`
          : 'Failed to authenticate with server. Is the backend running?'
      );
      console.error(err);
    }
  };

  const handleStandardAuth = async (e) => {
    e.preventDefault();
    setTouched({ email: true, password: true, name: !isLogin });
    if (!isFormValid) return;

    setServerError('');
    setNotVerifiedEmail('');
    setLoading(true);

    try {
      const endpoint = isLogin ? '/auth/login' : '/auth/register';
      const res = await axios.post(`${API_BASE}${endpoint}`, formData);

      if (!isLogin) {
        // Registration success: no JWT yet, show "check inbox" banner
        setPendingEmail(formData.email);
        return;
      }

      // Login success
      localStorage.setItem('token', res.data.access_token);
      localStorage.setItem('userName', res.data.user.name || formData.email.split('@')[0]);
      navigate('/');
    } catch (err) {
      const detail = err.response?.data?.detail;
      if (err.response?.status === 403 && detail === 'EMAIL_NOT_VERIFIED') {
        // Show the "not verified" banner with a resend option
        setNotVerifiedEmail(formData.email);
      } else {
        setServerError(detail || 'Authentication failed. Check logs.');
      }
    } finally {
      setLoading(false);
    }
  };

  const handleResend = async (emailToResend) => {
    setResendLoading(true);
    setResendSent(false);
    try {
      await axios.post(`${API_BASE}/auth/resend-verification`, { email: emailToResend });
      setResendSent(true);
    } catch (err) {
      setServerError(err.response?.data?.detail || 'Could not resend email. Try again later.');
    } finally {
      setResendLoading(false);
    }
  };

  const switchMode = () => {
    setIsLogin((v) => !v);
    setServerError('');
    setPendingEmail('');
    setNotVerifiedEmail('');
    setResendSent(false);
    setForgotMode(false);
    setForgotSent(false);
    setForgotEmail('');
    setForgotError('');
    setFormData({ email: '', password: '', name: '' });
    setTouched({ email: false, password: false, name: false });
  };

  // Inline error display
  const FieldError = ({ msg }) =>
    msg ? (
      <motion.p
        initial={{ opacity: 0, y: -4 }}
        animate={{ opacity: 1, y: 0 }}
        exit={{ opacity: 0, y: -4 }}
        className="text-red-400 text-xs mt-1 ml-1"
      >
        {msg}
      </motion.p>
    ) : null;

  // ── Forgot password handler ─────────────────────────────────────────────────
  const handleForgot = async (e) => {
    e.preventDefault();
    if (!forgotEmail.trim() || !/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(forgotEmail)) {
      setForgotError('Enter a valid email address.');
      return;
    }
    setForgotLoading(true);
    setForgotError('');
    try {
      await axios.post(`${API_BASE}/auth/forgot-password`, { email: forgotEmail });
      setForgotSent(true);
    } catch {
      setForgotError('Something went wrong. Please try again.');
    } finally {
      setForgotLoading(false);
    }
  };

  // ── "Forgot password" sent screen ──────────────────────────────────────────
  if (forgotMode && forgotSent) {
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
            <div className="w-16 h-16 bg-gradient-to-br from-blue-500 to-purple-600 rounded-2xl flex items-center justify-center mx-auto mb-5 shadow-lg shadow-blue-500/20">
              <MailCheck className="w-8 h-8 text-white" />
            </div>
            <h2 className="text-2xl font-bold text-white mb-2">Check your inbox!</h2>
            <p className="text-slate-400 text-sm mb-1">If an account exists for</p>
            <p className="text-blue-400 font-medium text-sm mb-4 break-all">{forgotEmail}</p>
            <p className="text-slate-500 text-xs mb-6">
              a password reset link has been sent. The link expires in 1 hour.
            </p>
            <button
              onClick={() => { setForgotMode(false); setForgotSent(false); setForgotEmail(''); }}
              className="text-slate-500 hover:text-slate-300 text-sm transition-colors"
            >
              ← Back to Sign In
            </button>
          </motion.div>
        </div>
      </div>
    );
  }

  // ── "Check your inbox" screen ───────────────────────────────────────────────
  if (pendingEmail) {
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
            <div className="w-16 h-16 bg-gradient-to-br from-blue-500 to-purple-600 rounded-2xl flex items-center justify-center mx-auto mb-5 shadow-lg shadow-blue-500/20">
              <MailCheck className="w-8 h-8 text-white" />
            </div>
            <h2 className="text-2xl font-bold text-white mb-2">Check your inbox!</h2>
            <p className="text-slate-400 text-sm mb-1">
              We sent a verification link to
            </p>
            <p className="text-blue-400 font-medium text-sm mb-6 break-all">{pendingEmail}</p>
            <p className="text-slate-500 text-xs mb-6">
              Click the link in the email to activate your account. The link expires in 1 hour.
            </p>

            {resendSent ? (
              <div className="text-green-400 bg-green-900/20 border border-green-900/50 rounded-xl p-3 text-sm mb-4">
                ✅ A new verification email has been sent!
              </div>
            ) : (
              <motion.button
                whileHover={{ scale: 1.02 }}
                whileTap={{ scale: 0.98 }}
                disabled={resendLoading}
                onClick={() => handleResend(pendingEmail)}
                className="w-full border border-slate-700 hover:border-purple-500/60 text-slate-300 hover:text-white font-medium py-3 rounded-xl flex items-center justify-center gap-2 transition-all mb-4 disabled:opacity-50"
              >
                {resendLoading ? (
                  <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" />
                ) : (
                  <>
                    <RefreshCw className="w-4 h-4" /> Resend verification email
                  </>
                )}
              </motion.button>
            )}

            <button onClick={switchMode} className="text-slate-500 hover:text-slate-300 text-sm transition-colors">
              ← Back to Sign In
            </button>
          </motion.div>
        </div>
      </div>
    );
  }

  // ── Main login / register form ──────────────────────────────────────────────
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

          {/* ── Forgot password form ──────────────────────────────────── */}
          {forgotMode ? (
            <>
              <h1 className="text-2xl font-bold text-center text-white mb-2">Forgot Password?</h1>
              <p className="text-slate-400 text-center text-sm mb-8">
                Enter your email and we'll send you a reset link.
              </p>
              <form onSubmit={handleForgot} noValidate className="space-y-4 mb-6">
                <div className="space-y-1">
                  <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">Email</label>
                  <div className="relative">
                    <Mail className="absolute left-3 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-600" />
                    <input
                      id="forgot-email"
                      type="email"
                      value={forgotEmail}
                      onChange={(e) => { setForgotEmail(e.target.value); setForgotError(''); }}
                      className="w-full bg-slate-950/50 border border-slate-800 focus:border-purple-500 rounded-xl pl-10 pr-4 py-3 text-white placeholder-slate-600 outline-none transition-colors"
                      placeholder="Enter your email"
                      autoComplete="email"
                      autoFocus
                    />
                  </div>
                  {forgotError && (
                    <p className="text-red-400 text-xs mt-1 ml-1">{forgotError}</p>
                  )}
                </div>

                <motion.button
                  whileHover={{ scale: 1.02 }}
                  whileTap={{ scale: 0.98 }}
                  disabled={forgotLoading}
                  type="submit"
                  className="w-full bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-500 hover:to-purple-500 text-white font-medium py-3 rounded-xl flex items-center justify-center gap-2 transition-all disabled:opacity-50"
                >
                  {forgotLoading ? (
                    <div className="w-5 h-5 border-2 border-white/30 border-t-white rounded-full animate-spin" />
                  ) : (
                    <><KeyRound className="w-4 h-4" /> Send Reset Link</>
                  )}
                </motion.button>
              </form>

              <div className="text-center">
                <button
                  onClick={() => { setForgotMode(false); setForgotError(''); }}
                  className="text-slate-400 hover:text-white text-sm transition-colors"
                >
                  ← Back to Sign In
                </button>
              </div>
            </>
          ) : (
          <>

          <h1 className="text-3xl font-bold text-center text-white mb-2">BioQuantum</h1>
          <p className="text-slate-400 text-center text-sm mb-8">
            Hybrid classical-quantum genomic search platform.
          </p>

          <form onSubmit={handleStandardAuth} noValidate className="space-y-4 mb-6">

            <AnimatePresence mode="popLayout">
              {/* Name – register only */}
              {!isLogin && (
                <motion.div
                  key="name-field"
                  initial={{ opacity: 0, height: 0 }}
                  animate={{ opacity: 1, height: 'auto' }}
                  exit={{ opacity: 0, height: 0 }}
                >
                  <div className="space-y-1">
                    <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">
                      Full Name <span className="text-slate-600 normal-case">(Optional)</span>
                    </label>
                    <input
                      type="text"
                      value={formData.name}
                      onChange={(e) => handleChange('name', e.target.value)}
                      onBlur={() => handleBlur('name')}
                      maxLength={61}
                      className={`w-full bg-slate-950/50 border rounded-xl px-4 py-3 text-white placeholder-slate-600 outline-none transition-colors
                        ${fieldErrors.name ? 'border-red-500 focus:border-red-400' : 'border-slate-800 focus:border-purple-500'}`}
                      placeholder="Enter your name"
                    />
                    <AnimatePresence><FieldError msg={fieldErrors.name} /></AnimatePresence>
                  </div>
                </motion.div>
              )}
            </AnimatePresence>

            {/* Email */}
            <div className="space-y-1">
              <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">Email</label>
              <div className="relative">
                <Mail className="absolute left-3 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-600" />
                <input
                  type="email"
                  value={formData.email}
                  onChange={(e) => handleChange('email', e.target.value)}
                  onBlur={() => handleBlur('email')}
                  className={`w-full bg-slate-950/50 border rounded-xl pl-10 pr-4 py-3 text-white placeholder-slate-600 outline-none transition-colors
                    ${fieldErrors.email ? 'border-red-500 focus:border-red-400' : 'border-slate-800 focus:border-purple-500'}`}
                  placeholder="Enter your email"
                  autoComplete="email"
                />
              </div>
              <AnimatePresence><FieldError msg={fieldErrors.email} /></AnimatePresence>
            </div>

            {/* Password */}
            <div className="space-y-1">
              <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">Password</label>
              <div className="relative">
                <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-600" />
                <input
                  type="password"
                  value={formData.password}
                  onChange={(e) => handleChange('password', e.target.value)}
                  onBlur={() => handleBlur('password')}
                  className={`w-full bg-slate-950/50 border rounded-xl pl-10 pr-4 py-3 text-white placeholder-slate-600 outline-none transition-colors
                    ${fieldErrors.password ? 'border-red-500 focus:border-red-400' : 'border-slate-800 focus:border-purple-500'}`}
                  placeholder="Enter your password"
                  autoComplete={isLogin ? 'current-password' : 'new-password'}
                />
              </div>
              <AnimatePresence><FieldError msg={fieldErrors.password} /></AnimatePresence>

              {/* Password strength hints – register only */}
              {!isLogin && (
                <ul className="mt-2 space-y-1 pl-1">
                  {[
                    { rule: formData.password.length >= 8, label: 'At least 8 characters' },
                    { rule: /[A-Z]/.test(formData.password), label: 'One uppercase letter' },
                    { rule: /[0-9!@#$%^&*]/.test(formData.password), label: 'One number or special character' },
                  ].map(({ rule, label }) => (
                    <li key={label} className={`text-xs flex items-center gap-1.5 transition-colors ${rule ? 'text-green-400' : 'text-slate-500'}`}>
                      <span>{rule ? '✓' : '○'}</span>{label}
                    </li>
                  ))}
                </ul>
              )}

              {/* Forgot password link – login mode only */}
              {isLogin && (
                <div className="text-right mt-1">
                  <button
                    type="button"
                    onClick={() => { setForgotMode(true); setForgotEmail(formData.email); }}
                    className="text-xs text-slate-500 hover:text-purple-400 transition-colors"
                  >
                    Forgot password?
                  </button>
                </div>
              )}
            </div>

            {/* Email not verified banner (login only) */}
            {notVerifiedEmail && (
              <motion.div
                initial={{ opacity: 0, y: -4 }}
                animate={{ opacity: 1, y: 0 }}
                className="bg-amber-900/20 border border-amber-800/50 rounded-xl p-4 space-y-2"
              >
                <p className="text-amber-400 text-sm font-medium">📬 Email not verified</p>
                <p className="text-amber-300/70 text-xs">
                  Please check your inbox and click the verification link before signing in.
                </p>
                {resendSent ? (
                  <p className="text-green-400 text-xs">✅ A new verification email has been sent!</p>
                ) : (
                  <button
                    type="button"
                    disabled={resendLoading}
                    onClick={() => handleResend(notVerifiedEmail)}
                    className="flex items-center gap-1 text-xs text-amber-400 hover:text-amber-300 transition-colors disabled:opacity-50"
                  >
                    <RefreshCw className="w-3 h-3" />
                    {resendLoading ? 'Sending…' : 'Resend verification email'}
                  </button>
                )}
              </motion.div>
            )}

            {/* Generic server error */}
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
                <>
                  {isLogin ? 'Sign In' : 'Create Account'} <ArrowRight className="w-4 h-4" />
                </>
              )}
            </motion.button>
          </form>

          {/* Divider */}
          <div className="flex items-center gap-4 mb-6">
            <div className="h-px bg-slate-800 flex-1" />
            <span className="text-slate-500 text-xs font-medium uppercase tracking-wider">OR</span>
            <div className="h-px bg-slate-800 flex-1" />
          </div>

          {/* Google OAuth */}
          <div className="flex justify-center mb-6">
            <GoogleLogin
              onSuccess={handleGoogleSuccess}
              onError={() => setServerError('Google Login Widget Error')}
              theme="filled_black"
              shape="pill"
              text={isLogin ? 'signin_with' : 'signup_with'}
            />
          </div>

          {/* Toggle login / register */}
          <div className="text-center">
            <button onClick={switchMode} className="text-slate-400 hover:text-white text-sm transition-colors">
              {isLogin ? "Don't have an account? Sign up" : 'Already have an account? Sign in'}
            </button>
          </div>
          </>
          )}
        </motion.div>
      </div>
    </div>
  );
}
