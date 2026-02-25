import React from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate, useNavigate, useLocation } from 'react-router-dom';
import { LogOut, LayoutDashboard, History as HistoryIcon, Dna } from 'lucide-react';
import Login from './pages/Login';
import Dashboard from './pages/Dashboard';
import History from './pages/History';

function ProtectedRoute({ children }) {
  const token = localStorage.getItem('token');
  if (!token) return <Navigate to="/login" replace />;
  return children;
}

function Navigation() {
  const navigate = useNavigate();
  const location = useLocation();
  const userName = localStorage.getItem('userName') || 'User';

  const handleLogout = () => {
    localStorage.removeItem('token');
    localStorage.removeItem('userName');
    navigate('/login');
  };

  if (location.pathname === '/login') return null;

  return (
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
  );
}

function App() {
  return (
    <Router>
      <div className="min-h-screen bg-slate-950 text-slate-200 font-sans selection:bg-purple-900/50 flex flex-col">
        <Navigation />
        <main className="flex-1 w-full relative">
          <Routes>
            <Route path="/login" element={<Login />} />
            <Route path="/" element={<ProtectedRoute><Dashboard /></ProtectedRoute>} />
            <Route path="/history" element={<ProtectedRoute><History /></ProtectedRoute>} />
            <Route path="*" element={<Navigate to="/" replace />} />
          </Routes>
        </main>
      </div>
    </Router>
  );
}

export default App;
