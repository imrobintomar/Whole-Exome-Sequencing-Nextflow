'use client';

import { useEffect, useState } from 'react';
import { onAuthStateChanged, signOut } from 'firebase/auth';
import { auth } from '@/lib/firebase';
import { User } from '@/lib/api';
import { Button } from '@/components/ui/button';
import { Dna } from 'lucide-react';
import LoginForm from '@/components/LoginForm';
import RegisterForm from '@/components/RegisterForm';
import Dashboard from '@/components/Dashboard';
import LandingPage from '@/components/LandingPage';
import AboutPage from '@/components/AboutPage';
import ResearchPage from '@/components/ResearchPage';
import ContactPage from '@/components/ContactPage';

type PageType = 'home' | 'about' | 'research' | 'contact' | 'signin';

export default function Home() {
  const [user, setUser] = useState<User | null>(null);
  const [showRegister, setShowRegister] = useState(false);
  const [loading, setLoading] = useState(true);
  const [currentPage, setCurrentPage] = useState<PageType>('home');

  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (firebaseUser) => {
      if (firebaseUser) {
        setUser({
          uid: firebaseUser.uid,
          email: firebaseUser.email || '',
          displayName: firebaseUser.displayName || firebaseUser.email || '',
        });
      } else {
        setUser(null);
      }
      setLoading(false);
    });

    return () => unsubscribe();
  }, []);

  const handleLogin = () => {
    // User state will be updated by onAuthStateChanged
  };

  const handleRegister = () => {
    // User state will be updated by onAuthStateChanged
  };

  const handleLogout = async () => {
    await signOut(auth);
  };

  const renderPage = () => {
    switch (currentPage) {
      case 'about':
        return <AboutPage />;
      case 'research':
        return <ResearchPage />;
      case 'contact':
        return <ContactPage />;
      case 'signin':
        return (
          <div className="min-h-screen flex items-center justify-center bg-gradient-to-br from-blue-50 to-indigo-100">
            <div className="w-full max-w-md">
              {showRegister ? (
                <RegisterForm
                  onRegister={handleRegister}
                  onToggle={() => setShowRegister(false)}
                />
              ) : (
                <LoginForm
                  onLogin={handleLogin}
                  onToggle={() => setShowRegister(true)}
                />
              )}
            </div>
          </div>
        );
      default:
        return <LandingPage onNavigate={setCurrentPage} onSignIn={() => setCurrentPage('signin')} />;
    }
  };

  if (loading) {
    return (
      <div className="min-h-screen flex items-center justify-center">
        <div className="text-xl">Loading...</div>
      </div>
    );
  }

  // If user is logged in, show dashboard
  if (user) {
    return <Dashboard user={user} onLogout={handleLogout} />;
  }

  // Show landing page with navigation
  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-white to-slate-50 flex flex-col">
      {/* Header */}
      <header className="border-b border-slate-200 bg-white/80 backdrop-blur-sm sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-4 py-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <div className="w-10 h-10 bg-gradient-to-br from-blue-600 to-cyan-600 rounded-lg flex items-center justify-center cursor-pointer" onClick={() => setCurrentPage('home')}>
                <Dna className="h-5 w-5 text-white" />
              </div>
              <span className="text-xl font-bold text-slate-900 cursor-pointer" onClick={() => setCurrentPage('home')}>
                ATGCFlow
              </span>
            </div>

            <nav className="hidden md:flex items-center gap-8">
              <button
                onClick={() => setCurrentPage('home')}
                className={`transition ${currentPage === 'home' ? 'text-blue-600 font-semibold' : 'text-slate-600 hover:text-slate-900'}`}
              >
                Home
              </button>
              <button
                onClick={() => setCurrentPage('about')}
                className={`transition ${currentPage === 'about' ? 'text-blue-600 font-semibold' : 'text-slate-600 hover:text-slate-900'}`}
              >
                About
              </button>
              <button
                onClick={() => setCurrentPage('research')}
                className={`transition ${currentPage === 'research' ? 'text-blue-600 font-semibold' : 'text-slate-600 hover:text-slate-900'}`}
              >
                Research
              </button>
              <button
                onClick={() => setCurrentPage('contact')}
                className={`transition ${currentPage === 'contact' ? 'text-blue-600 font-semibold' : 'text-slate-600 hover:text-slate-900'}`}
              >
                Contact
              </button>
            </nav>

            <div className="flex items-center gap-4">
              <Button className="bg-blue-600 hover:bg-blue-700 text-white" onClick={() => setCurrentPage('signin')}>
                Sign In
              </Button>
            </div>
          </div>
        </div>
      </header>

      {/* Page Content */}
      <main className="flex-grow">
        {renderPage()}
      </main>

      {/* Footer */}
      <footer className="border-t border-slate-200 bg-white mt-12">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-12">
          <div className="grid md:grid-cols-4 gap-8 mb-8">
            <div>
              <div className="flex items-center gap-2 mb-4">
                <div className="w-8 h-8 bg-gradient-to-br from-blue-600 to-cyan-600 rounded-lg flex items-center justify-center">
                  <Dna className="h-4 w-4 text-white" />
                </div>
                <span className="font-bold text-slate-900">ATGCFlow</span>
              </div>
              <p className="text-sm text-slate-600">Industry-grade whole exome sequencing analysis platform</p>
            </div>

            <div>
              <h4 className="font-semibold text-slate-900 mb-4">Navigation</h4>
              <ul className="space-y-2">
                <li><button onClick={() => setCurrentPage('home')} className="text-sm text-slate-600 hover:text-slate-900">Home</button></li>
                <li><button onClick={() => setCurrentPage('about')} className="text-sm text-slate-600 hover:text-slate-900">About</button></li>
                <li><button onClick={() => setCurrentPage('research')} className="text-sm text-slate-600 hover:text-slate-900">Research</button></li>
                <li><button onClick={() => setCurrentPage('contact')} className="text-sm text-slate-600 hover:text-slate-900">Contact</button></li>
              </ul>
            </div>

            

            <div>
              <h4 className="font-semibold text-slate-900 mb-4">Legal</h4>
              <ul className="space-y-2">
                <li><a href="#" className="text-sm text-slate-600 hover:text-slate-900">Privacy Policy</a></li>
                <li><a href="#" className="text-sm text-slate-600 hover:text-slate-900">Terms of Service</a></li>
                <li><a href="#" className="text-sm text-slate-600 hover:text-slate-900">Disclaimer</a></li>
              </ul>
            </div>
          </div>

          <div className="border-t border-slate-200 pt-8">
            <div className="text-center">
              <div className="text-sm text-slate-600">
                © 2025 ATGCFlow || Research Project || Not for Clinical Use
              </div>
            </div>
          </div>
        </div>
      </footer>
    </div>
  );
}
