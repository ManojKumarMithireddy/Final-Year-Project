import React, { useEffect, useState } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import { motion } from 'framer-motion';
import { Dna, CheckCircle, XCircle, Loader } from 'lucide-react';
import axios from 'axios';

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000/api';

export default function VerifyEmail() {
  const [searchParams] = useSearchParams();
  const navigate = useNavigate();
  const [status, setStatus] = useState('loading'); // 'loading' | 'success' | 'error'
  const [message, setMessage] = useState('');
  const [resendEmail, setResendEmail] = useState('');
  const [resendSent, setResendSent] = useState(false);
  const [resendLoading, setResendLoading] = useState(false);

  useEffect(() => {
    const token = searchParams.get('token');
    if (!token) {
      setStatus('error');
      setMessage('No verification token found in URL.');
      return;
    }

    axios
      .get(`${API_BASE}/auth/verify-email`, { params: { token } })
      .then((res) => {
        localStorage.setItem('token', res.data.access_token);
        localStorage.setItem('userName', res.data.user.name || res.data.user.email);
        localStorage.setItem('loginTimestamp', Date.now().toString());
        setStatus('success');
        // Auto-redirect after 3 seconds
        setTimeout(() => navigate('/'), 3000);
      })
      .catch((err) => {
        setStatus('error');
        setMessage(
          err.response?.data?.detail || 'Verification failed. The link may have expired.'
        );
        // If server returned the email, pre-fill resend input
        if (err.response?.data?.email) {
          setResendEmail(err.response.data.email);
        }
      });
  }, []);

  const handleResend = async (e) => {
    e.preventDefault();
    if (!resendEmail) return;
    setResendLoading(true);
    try {
      await axios.post(`${API_BASE}/auth/resend-verification`, { email: resendEmail });
      setResendSent(true);
    } catch (err) {
      setMessage(err.response?.data?.detail || 'Could not resend. Please try again.');
    } finally {
      setResendLoading(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col items-center justify-center p-6 bg-slate-950 relative overflow-hidden">
      {/* Background blobs */}
      <div className="absolute top-0 left-0 w-[500px] h-[500px] bg-blue-600/10 rounded-full blur-[120px] -translate-x-1/2 -translate-y-1/2 pointer-events-none" />
      <div className="absolute bottom-0 right-0 w-[500px] h-[500px] bg-purple-600/10 rounded-full blur-[120px] translate-x-1/2 translate-y-1/2 pointer-events-none" />

      <div className="relative z-10 w-full max-w-md">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-slate-900/60 border border-slate-800 p-8 rounded-3xl shadow-2xl backdrop-blur-xl text-center"
        >
          {/* Logo */}
          <div className="flex justify-center mb-6">
            <div className="w-16 h-16 bg-gradient-to-br from-blue-500 to-purple-600 rounded-2xl flex items-center justify-center shadow-lg shadow-blue-500/20">
              <Dna className="w-8 h-8 text-white" />
            </div>
          </div>

          {/* Loading */}
          {status === 'loading' && (
            <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }}>
              <Loader className="w-10 h-10 text-blue-400 mx-auto mb-4 animate-spin" />
              <h2 className="text-xl font-bold text-white mb-2">Verifying your email…</h2>
              <p className="text-slate-400 text-sm">Please wait a moment.</p>
            </motion.div>
          )}

          {/* Success */}
          {status === 'success' && (
            <motion.div initial={{ opacity: 0, scale: 0.9 }} animate={{ opacity: 1, scale: 1 }}>
              <CheckCircle className="w-14 h-14 text-green-400 mx-auto mb-4" />
              <h2 className="text-2xl font-bold text-white mb-2">Email verified! 🎉</h2>
              <p className="text-slate-400 text-sm mb-6">
                Your account is now active. Redirecting to the dashboard…
              </p>
              <motion.button
                whileHover={{ scale: 1.02 }}
                whileTap={{ scale: 0.98 }}
                onClick={() => navigate('/')}
                className="w-full bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-500 hover:to-purple-500 text-white font-medium py-3 rounded-xl transition-all"
              >
                Go to Dashboard
              </motion.button>
            </motion.div>
          )}

          {/* Error */}
          {status === 'error' && (
            <motion.div initial={{ opacity: 0, scale: 0.9 }} animate={{ opacity: 1, scale: 1 }}>
              <XCircle className="w-14 h-14 text-red-400 mx-auto mb-4" />
              <h2 className="text-xl font-bold text-white mb-2">Verification failed</h2>
              <p className="text-red-400 text-sm mb-6 bg-red-900/20 border border-red-900/50 rounded-xl p-3">
                {message}
              </p>

              {!resendSent ? (
                <form onSubmit={handleResend} className="space-y-3">
                  <p className="text-slate-400 text-sm">
                    Need a new link? Enter your email and we'll resend it.
                  </p>
                  <input
                    type="email"
                    value={resendEmail}
                    onChange={(e) => setResendEmail(e.target.value)}
                    placeholder="your@email.com"
                    required
                    className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-white placeholder-slate-600 outline-none focus:border-purple-500 transition-colors"
                  />
                  <motion.button
                    type="submit"
                    disabled={resendLoading}
                    whileHover={{ scale: 1.02 }}
                    whileTap={{ scale: 0.98 }}
                    className="w-full bg-gradient-to-r from-blue-600 to-purple-600 text-white font-medium py-3 rounded-xl transition-all disabled:opacity-50"
                  >
                    {resendLoading ? (
                      <Loader className="w-5 h-5 animate-spin mx-auto" />
                    ) : (
                      'Resend verification email'
                    )}
                  </motion.button>
                </form>
              ) : (
                <div className="text-green-400 bg-green-900/20 border border-green-900/50 rounded-xl p-3 text-sm">
                  ✅ A new verification email has been sent. Check your inbox!
                </div>
              )}

              <button
                onClick={() => navigate('/login')}
                className="mt-4 text-slate-500 hover:text-slate-300 text-sm transition-colors"
              >
                ← Back to Login
              </button>
            </motion.div>
          )}
        </motion.div>
      </div>
    </div>
  );
}
