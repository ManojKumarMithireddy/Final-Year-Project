import React from 'react';
import { AlertTriangle, RefreshCw } from 'lucide-react';

/**
 * Error Boundary — catches unhandled render errors and shows a
 * friendly fallback UI instead of a white screen.
 */
class ErrorBoundary extends React.Component {
    constructor(props) {
        super(props);
        this.state = { hasError: false, error: null };
    }

    static getDerivedStateFromError(error) {
        return { hasError: true, error };
    }

    componentDidCatch(error, errorInfo) {
        console.error('ErrorBoundary caught:', error, errorInfo);
    }

    handleReset = () => {
        this.setState({ hasError: false, error: null });
    };

    render() {
        if (this.state.hasError) {
            return (
                <div className="min-h-[60vh] flex items-center justify-center p-8">
                    <div className="max-w-md text-center space-y-6">
                        <div className="w-16 h-16 rounded-2xl bg-red-500/20 flex items-center justify-center mx-auto">
                            <AlertTriangle className="w-8 h-8 text-red-400" />
                        </div>
                        <div>
                            <h2 className="text-xl font-bold text-white mb-2">Something went wrong</h2>
                            <p className="text-slate-400 text-sm leading-relaxed">
                                An unexpected error occurred. Try refreshing or navigating back.
                            </p>
                        </div>
                        {this.state.error && (
                            <pre className="text-xs text-red-400/70 bg-red-950/20 border border-red-500/20 rounded-xl p-4 text-left overflow-auto max-h-32">
                                {this.state.error.toString()}
                            </pre>
                        )}
                        <div className="flex gap-3 justify-center">
                            <button
                                onClick={this.handleReset}
                                className="flex items-center gap-2 bg-slate-800 hover:bg-slate-700 text-white px-5 py-2.5 rounded-xl text-sm font-medium transition-colors"
                            >
                                <RefreshCw className="w-4 h-4" /> Try Again
                            </button>
                            <button
                                onClick={() => window.location.reload()}
                                className="bg-blue-600 hover:bg-blue-500 text-white px-5 py-2.5 rounded-xl text-sm font-medium transition-colors"
                            >
                                Reload Page
                            </button>
                        </div>
                    </div>
                </div>
            );
        }

        return this.props.children;
    }
}

export default ErrorBoundary;
