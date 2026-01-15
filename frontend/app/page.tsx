'use client';

import { useEffect, useState } from 'react';
import { onAuthStateChanged, signOut } from 'firebase/auth';
import { auth } from '@/lib/firebase';
import { User } from '@/lib/api';
import { Button } from '@/components/ui/button';
import { Dna, Menu, X } from 'lucide-react';
import LoginForm from '@/components/LoginForm';
import RegisterForm from '@/components/RegisterForm';
import Dashboard from '@/components/Dashboard';
import LandingPage from '@/components/LandingPage';
import AboutPage from '@/components/AboutPage';
import FeaturesPage from '@/components/FeaturesPage';
import PublicationsPage from '@/components/PublicationsPage';
import UseCasesPage from '@/components/UseCasesPage';
import PricingPage from '@/components/PricingPage';
import ResearchPage from '@/components/ResearchPage';
import ContactPage from '@/components/ContactPage';
import PrivacyPolicyPage from '@/components/PrivacyPolicyPage';
import TermsOfServicePage from '@/components/TermsOfServicePage';
import DisclaimerPage from '@/components/DisclaimerPage';
import EmailVerificationReminder from '@/components/EmailVerificationReminder';

type PageType = 'home' | 'about' | 'features' | 'usecases' | 'publications' | 'pricing' | 'research' | 'contact' | 'signin' | 'privacy' | 'terms' | 'disclaimer';

export default function Home() {
  const [user, setUser] = useState<User | null>(null);
  const [showRegister, setShowRegister] = useState(false);
  const [loading, setLoading] = useState(true);
  const [currentPage, setCurrentPage] = useState<PageType>('home');
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

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

  const handleNavigate = (page: PageType) => {
    setCurrentPage(page);
    setMobileMenuOpen(false);
    window.scrollTo({ top: 0, behavior: 'smooth' });
  };

  const renderPage = () => {
    switch (currentPage) {
      case 'about':
        return <AboutPage onNavigate={(page) => handleNavigate(page as PageType)} />;
      case 'features':
        return <FeaturesPage />;
      case 'usecases':
        return <UseCasesPage />;
      case 'publications':
        return <PublicationsPage />;
      case 'pricing':
        return <PricingPage onSignIn={() => handleNavigate('signin')} isAuthenticated={false} />;
      case 'research':
        return <ResearchPage />;
      case 'contact':
        return <ContactPage />;
      case 'privacy':
        return <PrivacyPolicyPage onBack={() => handleNavigate('home')} />;
      case 'terms':
        return <TermsOfServicePage onBack={() => handleNavigate('home')} />;
      case 'disclaimer':
        return <DisclaimerPage onBack={() => handleNavigate('home')} />;
      case 'signin':
        return (
          <div className="min-h-screen flex items-center justify-center bg-gradient-to-br from-purple-primary/10 via-slate-50 to-cyan/10">
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
        return <LandingPage onNavigate={handleNavigate} onSignIn={() => handleNavigate('signin')} />;
    }
  };

  if (loading) {
    return (
      <div className="min-h-screen flex items-center justify-center">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-purple-primary"></div>
      </div>
    );
  }

  // If user is logged in, check email verification
  if (user) {
    // Check if email is verified
    const firebaseUser = auth.currentUser;
    if (firebaseUser && !firebaseUser.emailVerified) {
      return <EmailVerificationReminder email={user.email} />;
    }
    return <Dashboard user={user} onLogout={handleLogout} />;
  }

  // Show landing page with navigation
  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-white to-slate-50 flex flex-col">
      {/* Header */}
      <header className="border-b border-slate-200 bg-white/90 backdrop-blur-sm sticky top-0 z-50 shadow-sm">
        <div className="max-w-7xl mx-auto px-4 py-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between">
            {/* Logo */}
            <div className="flex items-center gap-2">
              <div
                className="w-10 h-10 bg-gradient-to-br from-purple-primary to-cyan rounded-lg flex items-center justify-center cursor-pointer"
                onClick={() => handleNavigate('home')}
              >
                <Dna className="h-5 w-5 text-white" />
              </div>
              <span
                className="text-xl font-bold text-slate-900 cursor-pointer"
                onClick={() => handleNavigate('home')}
              >
                ATGC Flow
              </span>
            </div>

            {/* Desktop Navigation */}
            <nav className="hidden lg:flex items-center gap-6">
              <button
                onClick={() => handleNavigate('home')}
                className={`transition-colors ${currentPage === 'home' ? 'text-purple-primary font-semibold' : 'text-slate-600 hover:text-purple-primary'}`}
              >
                Home
              </button>
              <button
                onClick={() => handleNavigate('features')}
                className={`transition-colors ${currentPage === 'features' ? 'text-purple-primary font-semibold' : 'text-slate-600 hover:text-purple-primary'}`}
              >
                Features
              </button>
              <button
                onClick={() => handleNavigate('usecases')}
                className={`transition-colors ${currentPage === 'usecases' ? 'text-purple-primary font-semibold' : 'text-slate-600 hover:text-purple-primary'}`}
              >
                Use Cases
              </button>
              <button
                onClick={() => handleNavigate('publications')}
                className={`transition-colors ${currentPage === 'publications' ? 'text-purple-primary font-semibold' : 'text-slate-600 hover:text-purple-primary'}`}
              >
                Publications
              </button>
              <button
                onClick={() => handleNavigate('pricing')}
                className={`transition-colors ${currentPage === 'pricing' ? 'text-purple-primary font-semibold' : 'text-slate-600 hover:text-purple-primary'}`}
              >
                Pricing
              </button>
              <button
                onClick={() => handleNavigate('about')}
                className={`transition-colors ${currentPage === 'about' ? 'text-purple-primary font-semibold' : 'text-slate-600 hover:text-purple-primary'}`}
              >
                About
              </button>
              <button
                onClick={() => handleNavigate('contact')}
                className={`transition-colors ${currentPage === 'contact' ? 'text-purple-primary font-semibold' : 'text-slate-600 hover:text-purple-primary'}`}
              >
                Contact
              </button>
            </nav>

            {/* Desktop Sign In Button */}
            <div className="hidden lg:flex items-center gap-4">
              <Button
                className="bg-purple-primary hover:bg-purple-light text-white"
                onClick={() => handleNavigate('signin')}
              >
                Sign In
              </Button>
            </div>

            {/* Mobile Menu Button */}
            <button
              className="lg:hidden p-2 text-slate-600 hover:text-purple-primary"
              onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
            >
              {mobileMenuOpen ? <X className="h-6 w-6" /> : <Menu className="h-6 w-6" />}
            </button>
          </div>

          {/* Mobile Navigation */}
          {mobileMenuOpen && (
            <nav className="lg:hidden mt-4 pb-4 border-t border-slate-200 pt-4">
              <div className="flex flex-col gap-3">
                <button
                  onClick={() => handleNavigate('home')}
                  className={`text-left px-4 py-2 rounded-lg transition-colors ${currentPage === 'home' ? 'bg-purple-primary/10 text-purple-primary font-semibold' : 'text-slate-600 hover:bg-slate-100'}`}
                >
                  Home
                </button>
                <button
                  onClick={() => handleNavigate('features')}
                  className={`text-left px-4 py-2 rounded-lg transition-colors ${currentPage === 'features' ? 'bg-purple-primary/10 text-purple-primary font-semibold' : 'text-slate-600 hover:bg-slate-100'}`}
                >
                  Features
                </button>
                <button
                  onClick={() => handleNavigate('usecases')}
                  className={`text-left px-4 py-2 rounded-lg transition-colors ${currentPage === 'usecases' ? 'bg-purple-primary/10 text-purple-primary font-semibold' : 'text-slate-600 hover:bg-slate-100'}`}
                >
                  Use Cases
                </button>
                <button
                  onClick={() => handleNavigate('publications')}
                  className={`text-left px-4 py-2 rounded-lg transition-colors ${currentPage === 'publications' ? 'bg-purple-primary/10 text-purple-primary font-semibold' : 'text-slate-600 hover:bg-slate-100'}`}
                >
                  Publications
                </button>
                <button
                  onClick={() => handleNavigate('pricing')}
                  className={`text-left px-4 py-2 rounded-lg transition-colors ${currentPage === 'pricing' ? 'bg-purple-primary/10 text-purple-primary font-semibold' : 'text-slate-600 hover:bg-slate-100'}`}
                >
                  Pricing
                </button>
                <button
                  onClick={() => handleNavigate('about')}
                  className={`text-left px-4 py-2 rounded-lg transition-colors ${currentPage === 'about' ? 'bg-purple-primary/10 text-purple-primary font-semibold' : 'text-slate-600 hover:bg-slate-100'}`}
                >
                  About
                </button>
                <button
                  onClick={() => handleNavigate('contact')}
                  className={`text-left px-4 py-2 rounded-lg transition-colors ${currentPage === 'contact' ? 'bg-purple-primary/10 text-purple-primary font-semibold' : 'text-slate-600 hover:bg-slate-100'}`}
                >
                  Contact
                </button>
                <Button
                  className="bg-purple-primary hover:bg-purple-light text-white mt-2"
                  onClick={() => handleNavigate('signin')}
                >
                  Sign In
                </Button>
              </div>
            </nav>
          )}
        </div>
      </header>

      {/* Page Content */}
      <main className="flex-grow">
        {renderPage()}
      </main>

      {/* Footer */}
      <footer className="border-t border-purple-primary/20 bg-gradient-to-br from-slate-900 via-purple-primary to-purple-dark mt-12">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-12">
          <div className="grid md:grid-cols-4 gap-8 mb-8">
            {/* Brand */}
            <div>
              <div className="flex items-center gap-2 mb-4">
                <div className="w-8 h-8 bg-gradient-to-br from-cyan to-teal rounded-lg flex items-center justify-center">
                  <Dna className="h-4 w-4 text-white" />
                </div>
                <span className="font-bold text-white">ATGC Flow</span>
              </div>
              <p className="text-sm text-blue-100">
                Clinical-grade Whole Exome Sequencing analysis platform for researchers and clinicians worldwide
              </p>
            </div>

            {/* Product */}
            <div>
              <h4 className="font-semibold text-white mb-4">Product</h4>
              <ul className="space-y-2">
                <li><button onClick={() => handleNavigate('features')} className="text-sm text-blue-200 hover:text-white transition-colors">Features</button></li>
                <li><button onClick={() => handleNavigate('pricing')} className="text-sm text-blue-200 hover:text-white transition-colors">Pricing</button></li>
                <li><button onClick={() => handleNavigate('usecases')} className="text-sm text-blue-200 hover:text-white transition-colors">Use Cases</button></li>
                <li><button onClick={() => handleNavigate('research')} className="text-sm text-blue-200 hover:text-white transition-colors">Documentation</button></li>
              </ul>
            </div>

            {/* Company */}
            <div>
              <h4 className="font-semibold text-white mb-4">Company</h4>
              <ul className="space-y-2">
                <li><button onClick={() => handleNavigate('about')} className="text-sm text-blue-200 hover:text-white transition-colors">About</button></li>
                <li><button onClick={() => handleNavigate('publications')} className="text-sm text-blue-200 hover:text-white transition-colors">Publications</button></li>
                <li><button onClick={() => handleNavigate('contact')} className="text-sm text-blue-200 hover:text-white transition-colors">Contact</button></li>
              </ul>
            </div>

            {/* Legal */}
            <div>
              <h4 className="font-semibold text-white mb-4">Legal</h4>
              <ul className="space-y-2">
                <li><button onClick={() => handleNavigate('privacy')} className="text-sm text-blue-200 hover:text-white transition-colors">Privacy Policy</button></li>
                <li><button onClick={() => handleNavigate('terms')} className="text-sm text-blue-200 hover:text-white transition-colors">Terms of Service</button></li>
                <li><button onClick={() => handleNavigate('disclaimer')} className="text-sm text-blue-200 hover:text-white transition-colors">Disclaimer</button></li>
              </ul>
            </div>
          </div>

          {/* Bottom Bar */}
          <div className="border-t border-white/20 pt-8">
            <div className="text-center text-sm text-blue-100">
              © 2025 ATGC Flow | Research Project | Not for Clinical Use
            </div>
          </div>
        </div>
      </footer>
    </div>
  );
}
