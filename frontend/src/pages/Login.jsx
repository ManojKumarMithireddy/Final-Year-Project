import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { GoogleLogin } from '@react-oauth/google';
import axios from 'axios';
import { Dna, Lock, Shield, Server, Mail, ArrowRight } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

const API_BASE = 'http://localhost:8000/api';

export default function Login() {
  const navigate = useNavigate();
  const [error, setError] = useState('');
  const [isLogin, setIsLogin] = useState(true);
  const [formData, setFormData] = useState({ email: '', password: '', name: '' });
  const [loading, setLoading] = useState(false);

  const handleGoogleSuccess = async (credentialResponse) => {
    try {
      const res = await axios.post(`${API_BASE}/auth/google`, {
        token: credentialResponse.credential
      });
      localStorage.setItem('token', res.data.access_token);
      localStorage.setItem('userName', res.data.user.name);
      navigate('/');
    } catch (err) {
      if (err.response && err.response.data && err.response.data.detail) {
        setError(`Backend Error: ${err.response.data.detail}`);
      } else {
        setError('Failed to authenticate with server. Is the backend running?');
      }
      console.error(err);
    }
  };

  const handleStandardAuth = async (e) => {
    e.preventDefault();
    if (!formData.email || !formData.password) {
      setError('Please fill in all required fields');
      return;
    }
    setError('');
    setLoading(true);
    
    try {
      const endpoint = isLogin ? '/auth/login' : '/auth/register';
      const res = await axios.post(`${API_BASE}${endpoint}`, formData);
      localStorage.setItem('token', res.data.access_token);
      localStorage.setItem('userName', res.data.user.name || formData.email.split('@')[0]);
      navigate('/');
    } catch (err) {
      if (err.response && err.response.data && err.response.data.detail) {
        setError(err.response.data.detail);
      } else {
        setError('Authentication failed. Check logs.');
      }
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col items-center justify-center p-6 relative overflow-hidden bg-slate-950">
      {/* Background decorations */}
      <div className="absolute top-0 left-0 w-[500px] h-[500px] bg-blue-600/10 rounded-full blur-[120px] -translate-x-1/2 -translate-y-1/2 pointer-events-none" />
      <div className="absolute bottom-0 right-0 w-[500px] h-[500px] bg-purple-600/10 rounded-full blur-[120px] translate-x-1/2 translate-y-1/2 pointer-events-none" />

      <div className="relative z-10 w-full max-w-md">
        <motion.div 
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-slate-900/60 border border-slate-800 p-8 rounded-3xl shadow-2xl backdrop-blur-xl"
        >
          <div className="flex justify-center mb-6">
            <div className="w-20 h-20 bg-gradient-to-br from-blue-500 to-purple-600 rounded-2xl flex items-center justify-center shadow-lg shadow-blue-500/20">
              <Dna className="w-10 h-10 text-white" />
            </div>
          </div>
          
          <h1 className="text-3xl font-bold text-center text-white mb-2">BioQuantum</h1>
          <p className="text-slate-400 text-center text-sm mb-8">
            Hybrid classical-quantum genomic search platform.
          </p>

          <form onSubmit={handleStandardAuth} className="space-y-4 mb-6">
            <AnimatePresence mode="popLayout">
              {!isLogin && (
                <motion.div
                  initial={{ opacity: 0, height: 0 }}
                  animate={{ opacity: 1, height: 'auto' }}
                  exit={{ opacity: 0, height: 0 }}
                  className="space-y-1"
                >
                  <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">Full Name (Optional)</label>
                  <input 
                    type="text" 
                    value={formData.name}
                    onChange={(e) => setFormData({...formData, name: e.target.value})}
                    className="w-full bg-slate-950/50 border border-slate-800 rounded-xl px-4 py-3 text-white placeholder-slate-600 outline-none focus:border-purple-500 transition-colors"
                    placeholder="Jane Doe"
                  />
                </motion.div>
              )}
            </AnimatePresence>

            <div className="space-y-1">
              <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">Email</label>
              <div className="relative">
                <Mail className="absolute left-3 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-600" />
                <input 
                  type="email" 
                  value={formData.email}
                  onChange={(e) => setFormData({...formData, email: e.target.value})}
                  className="w-full bg-slate-950/50 border border-slate-800 rounded-xl pl-10 pr-4 py-3 text-white placeholder-slate-600 outline-none focus:border-purple-500 transition-colors"
                  placeholder="user@example.com"
                />
              </div>
            </div>

            <div className="space-y-1">
              <label className="text-xs font-medium text-slate-400 uppercase tracking-wider">Password</label>
              <div className="relative">
                <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-5 h-5 text-slate-600" />
                <input 
                  type="password" 
                  value={formData.password}
                  onChange={(e) => setFormData({...formData, password: e.target.value})}
                  className="w-full bg-slate-950/50 border border-slate-800 rounded-xl pl-10 pr-4 py-3 text-white placeholder-slate-600 outline-none focus:border-purple-500 transition-colors"
                  placeholder="••••••••"
                />
              </div>
            </div>

            {error && <div className="text-red-400 text-sm text-center bg-red-900/20 p-3 rounded-xl border border-red-900/50">{error}</div>}

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

          <div className="flex items-center gap-4 mb-6">
            <div className="h-px bg-slate-800 flex-1"></div>
            <span className="text-slate-500 text-xs font-medium uppercase tracking-wider">OR</span>
            <div className="h-px bg-slate-800 flex-1"></div>
          </div>

          <div className="flex justify-center mb-6">
            <GoogleLogin
              onSuccess={handleGoogleSuccess}
              onError={() => setError('Google Login Widget Error')}
              theme="filled_black"
              shape="pill"
              text={isLogin ? 'signin_with' : 'signup_with'}
            />
          </div>

          <div className="text-center">
            <button 
              onClick={() => { setIsLogin(!isLogin); setError(''); setFormData({email:'', password:'', name:''}); }}
              className="text-slate-400 hover:text-white text-sm transition-colors"
            >
              {isLogin ? "Don't have an account? Sign up" : "Already have an account? Sign in"}
            </button>
          </div>

        </motion.div>
      </div>
    </div>
  );
}
