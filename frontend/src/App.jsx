import React, { useEffect, useCallback, useRef, useState } from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate, useNavigate, useLocation } from 'react-router-dom';
import { LogOut, LayoutDashboard, History as HistoryIcon, Dna } from 'lucide-react';
import { Toaster } from 'react-hot-toast';
import axios from 'axios';
import Login from './pages/Login';
import Dashboard from './pages/Dashboard';
import History from './pages/History';
import VerifyEmail from './pages/VerifyEmail';
import ResetPassword from './pages/ResetPassword';
import SessionTimeoutModal from './components/SessionTimeoutModal';
import ErrorBoundary from './components/ErrorBoundary';

// ── Session configuration ─────────────────────────────────────────────────────
const SESSION_TIMEOUT_MS = 30 * 60 * 1000;  // 30 minutes
const WARNING_BEFORE_MS = 2 * 60 * 1000;   // show warning 2 min before expiry
const WARNING_SECONDS = WARNING_BEFORE_MS / 1000;

function ProtectedRoute({ children }) {
  const token = localStorage.getItem('token');
  if (!token) return <Navigate to="/login" replace />;
  return children;
}

// ── Axios 401 interceptor ─────────────────────────────────────────────────────
axios.interceptors.response.use(
  (response) => response,
  (error) => {
    if (error.response?.status === 401) {
      localStorage.removeItem('token');
      localStorage.removeItem('userName');
      localStorage.removeItem('loginTimestamp');
      window.location.href = '/login';
    }
    return Promise.reject(error);
  }
);

// ── Session timeout hook ──────────────────────────────────────────────────────
function useSessionTimeout(onTimeout) {
  const [showWarning, setShowWarning] = useState(false);
  const [secondsLeft, setSecondsLeft] = useState(WARNING_SECONDS);

  const logoutTimerRef = useRef(null);
  const warningTimerRef = useRef(null);
  const countdownRef = useRef(null);

  const clearAllTimers = useCallback(() => {
    if (logoutTimerRef.current) clearTimeout(logoutTimerRef.current);
    if (warningTimerRef.current) clearTimeout(warningTimerRef.current);
    if (countdownRef.current) clearInterval(countdownRef.current);
    logoutTimerRef.current = null;
    warningTimerRef.current = null;
    countdownRef.current = null;
  }, []);

  const startTimers = useCallback(() => {
    clearAllTimers();
    if (!localStorage.getItem('token')) return;

    // Warning fires 2 min before the hard logout
    warningTimerRef.current = setTimeout(() => {
      setShowWarning(true);
      setSecondsLeft(WARNING_SECONDS);

      // Tick down every second
      countdownRef.current = setInterval(() => {
        setSecondsLeft((prev) => {
          if (prev <= 1) {
            clearInterval(countdownRef.current);
            return 0;
          }
          return prev - 1;
        });
      }, 1000);
    }, SESSION_TIMEOUT_MS - WARNING_BEFORE_MS);

    // Hard logout
    logoutTimerRef.current = setTimeout(onTimeout, SESSION_TIMEOUT_MS);
  }, [clearAllTimers, onTimeout]);

  const extendSession = useCallback(() => {
    setShowWarning(false);
    setSecondsLeft(WARNING_SECONDS);
    localStorage.setItem('loginTimestamp', Date.now().toString());
    startTimers();
  }, [startTimers]);

  useEffect(() => {
    // Check for stale login on mount
    const loginTs = parseInt(localStorage.getItem('loginTimestamp'), 10);
    if (loginTs && Date.now() - loginTs > SESSION_TIMEOUT_MS) {
      onTimeout();
      return;
    }

    const events = ['mousemove', 'keydown', 'click', 'scroll', 'touchstart'];

    const resetOnActivity = () => {
      // Only reset if the warning is NOT showing (once the warning is visible,
      // the user must explicitly click "Extend Session")
      if (!warningTimerRef.current) return; // avoid double-reset
      if (showWarning) return;
      startTimers();
    };

    events.forEach((e) => window.addEventListener(e, resetOnActivity, { passive: true }));
    startTimers();

    return () => {
      clearAllTimers();
      events.forEach((e) => window.removeEventListener(e, resetOnActivity));
    };
  }, [startTimers, clearAllTimers, onTimeout, showWarning]);

  return { showWarning, secondsLeft, extendSession };
}

function Navigation() {
  const navigate = useNavigate();
  const location = useLocation();
  const userName = localStorage.getItem('userName') || 'User';

  const handleLogout = useCallback(() => {
    localStorage.removeItem('token');
    localStorage.removeItem('userName');
    localStorage.removeItem('loginTimestamp');
    navigate('/login');
  }, [navigate]);

  const { showWarning, secondsLeft, extendSession } = useSessionTimeout(handleLogout);

  if (location.pathname === '/login' || location.pathname === '/verify-email' || location.pathname === '/reset-password') return null;

  return (
    <>
      <nav className="bg-slate-950/80 backdrop-blur-md border-b border-slate-800 sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-6 h-16 flex justify-between items-center">
          <div className="flex items-center gap-3">
            <Dna className="w-7 h-7 text-blue-500" />
            <span className="text-xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-blue-400 to-purple-500">
              BioQuantum
            </span>
          </div>
          <div className="flex items-center gap-6">
            <button
              onClick={() => navigate('/')}
              className={`flex items-center gap-2 text-sm font-medium transition-colors ${location.pathname === '/' ? 'text-blue-400' : 'text-slate-400 hover:text-white'}`}
            >
              <LayoutDashboard className="w-4 h-4" /> Dashboard
            </button>
            <button
              onClick={() => navigate('/history')}
              className={`flex items-center gap-2 text-sm font-medium transition-colors ${location.pathname === '/history' ? 'text-purple-400' : 'text-slate-400 hover:text-white'}`}
            >
              <HistoryIcon className="w-4 h-4" /> History
            </button>
            <div className="h-6 w-px bg-slate-800"></div>
            <div className="flex items-center gap-4 text-sm">
              <span className="text-slate-300 font-medium">{userName}</span>
              <button
                onClick={handleLogout}
                className="text-slate-400 hover:text-red-400 transition-colors bg-slate-900 border border-slate-800 p-2 rounded-lg"
                title="Logout"
              >
                <LogOut className="w-4 h-4" />
              </button>
            </div>
          </div>
        </div>
      </nav>

      {/* Session timeout warning modal */}
      <SessionTimeoutModal
        visible={showWarning}
        secondsLeft={secondsLeft}
        totalWarning={WARNING_SECONDS}
        onExtend={extendSession}
        onLogout={handleLogout}
      />
    </>
  );
}

function App() {
  return (
    <Router>
      <div className="min-h-screen bg-slate-950 text-slate-200 font-sans selection:bg-purple-900/50 flex flex-col">
        <Navigation />
        <main className="flex-1 w-full relative">
          <ErrorBoundary>
            <Routes>
              <Route path="/login" element={<Login />} />
              <Route path="/verify-email" element={<VerifyEmail />} />
              <Route path="/reset-password" element={<ResetPassword />} />
              <Route path="/" element={<ProtectedRoute><Dashboard /></ProtectedRoute>} />
              <Route path="/history" element={<ProtectedRoute><History /></ProtectedRoute>} />
              <Route path="*" element={<Navigate to="/" replace />} />
            </Routes>
          </ErrorBoundary>
        </main>
        <Toaster
          position="top-right"
          toastOptions={{
            className: 'bg-slate-800 text-slate-200 border border-slate-700 text-sm',
            duration: 4000,
            style: { background: '#1e293b', color: '#e2e8f0', border: '1px solid #334155' },
          }}
        />
      </div>
    </Router>
  );
}

export default App;
