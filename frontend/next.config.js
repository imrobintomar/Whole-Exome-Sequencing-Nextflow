/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,

  // Firebase Authentication URL rewrites
  async rewrites() {
    return [
      {
        source: '/__/auth/:path*',
        destination: 'https://variant-ac1c6.firebaseapp.com/__/auth/:path*',
      },
    ];
  },
}

module.exports = nextConfig
