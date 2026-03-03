/**
 * Polls the backend health endpoint and returns { status, retry }.
 * status: 'checking' | 'online' | 'offline' | 'waking'
 *
 * - First attempt is immediate.
 * - If it fails, flips to 'waking' and retries every POLL_MS until it succeeds.
 * - Stops polling once 'online'.
 */
import { useState, useEffect, useCallback, useRef } from 'react';
import axios from 'axios';
import { API_BASE } from './api';

const HEALTH_URL = API_BASE.replace(/\/api\/?$/, '') + '/docs';
const POLL_MS    = 8000;   // retry interval while sleeping
const TIMEOUT_MS = 6000;   // per-request timeout

export function useBackendStatus() {
  const [status, setStatus] = useState('checking');
  const intervalRef = useRef(null);

  const ping = useCallback(async () => {
    try {
      await axios.get(HEALTH_URL, { timeout: TIMEOUT_MS });
      setStatus('online');
      if (intervalRef.current) {
        clearInterval(intervalRef.current);
        intervalRef.current = null;
      }
    } catch {
      setStatus((prev) => (prev === 'checking' ? 'offline' : 'waking'));
    }
  }, []);

  const startPolling = useCallback(() => {
    if (intervalRef.current) return;
    intervalRef.current = setInterval(ping, POLL_MS);
  }, [ping]);

  // Manual retry (called from the banner button)
  const retry = useCallback(() => {
    setStatus('waking');
    ping();
    startPolling();
  }, [ping, startPolling]);

  useEffect(() => {
    ping(); // immediate first check
    return () => {
      if (intervalRef.current) clearInterval(intervalRef.current);
    };
  }, [ping]);

  // Start polling automatically once we know it's offline/waking
  useEffect(() => {
    if (status === 'offline' || status === 'waking') startPolling();
    if (status === 'online' && intervalRef.current) {
      clearInterval(intervalRef.current);
      intervalRef.current = null;
    }
  }, [status, startPolling]);

  return { status, retry };
}
