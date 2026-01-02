import { initializeApp, getApps } from 'firebase/app';
import { getAuth } from 'firebase/auth';
import { getAnalytics } from 'firebase/analytics';

const firebaseConfig = {
  apiKey: "AIzaSyAcMhwJRa0diLwHw2bZaCLIr_akFZPAsWc",
  authDomain: "variant-ac1c6.firebaseapp.com",
  projectId: "variant-ac1c6",
  storageBucket: "variant-ac1c6.firebasestorage.app",
  messagingSenderId: "290831279354",
  appId: "1:290831279354:web:d1a3bb9be9a8345a455543",
  measurementId: "G-RYB1JVQGEE"
};

// Initialize Firebase only if it hasn't been initialized already
const app = getApps().length === 0 ? initializeApp(firebaseConfig) : getApps()[0];

// Initialize Firebase Authentication
export const auth = getAuth(app);

// Initialize Analytics (only in browser)
export const analytics = typeof window !== 'undefined' ? getAnalytics(app) : null;

export default app;
