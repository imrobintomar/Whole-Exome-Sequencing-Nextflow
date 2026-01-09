import type { Metadata } from 'next'
import { Inter } from 'next/font/google'
import './globals.css'
import { ThemeProvider } from '@/components/theme-provider'

const inter = Inter({ subsets: ['latin'] })

export const metadata: Metadata = {
  title: {
    default: 'ATGCFlow - Industry-Grade Whole Exome Sequencing Analysis Platform',
    template: '%s | ATGCFlow'
  },
  description: 'ATGCFlow is a comprehensive whole exome sequencing (WES) analysis platform featuring GATK best practices, ACMG variant classification, IGV browser integration, and real-time pipeline monitoring for genomic research.',
  keywords: [
    'whole exome sequencing',
    'WES analysis',
    'genomic sequencing',
    'variant analysis',
    'GATK pipeline',
    'ACMG classification',
    'bioinformatics',
    'genomics platform',
    'DNA sequencing',
    'variant calling',
    'clinical genomics',
    'NGS analysis',
    'exome analysis',
    'genomic variants',
    'precision medicine'
  ],
  authors: [{ name: 'ATGCFlow Team' }],
  creator: 'ATGCFlow',
  publisher: 'ATGCFlow',
  metadataBase: new URL(process.env.NEXT_PUBLIC_BASE_URL || 'https://atgcflow.com'),
  alternates: {
    canonical: '/',
  },
  openGraph: {
    type: 'website',
    locale: 'en_US',
    url: '/',
    title: 'ATGCFlow - Industry-Grade Whole Exome Sequencing Analysis Platform',
    description: 'Comprehensive whole exome sequencing analysis with GATK best practices, ACMG variant classification, and real-time monitoring for genomic research.',
    siteName: 'ATGCFlow',
    images: [
      {
        url: '/og-image.png',
        width: 1200,
        height: 630,
        alt: 'ATGCFlow - Whole Exome Sequencing Platform',
      },
    ],
  },
  twitter: {
    card: 'summary_large_image',
    title: 'ATGCFlow - Whole Exome Sequencing Analysis Platform',
    description: 'Industry-grade WES analysis with GATK, ACMG classification, and IGV browser integration.',
    images: ['/twitter-image.png'],
    creator: '@atgcflow',
  },
  robots: {
    index: true,
    follow: true,
    googleBot: {
      index: true,
      follow: true,
      'max-video-preview': -1,
      'max-image-preview': 'large',
      'max-snippet': -1,
    },
  },
  verification: {
    google: 'your-google-verification-code',
    // yandex: 'your-yandex-verification-code',
    // bing: 'your-bing-verification-code',
  },
  category: 'Bioinformatics',
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  const jsonLd = {
    '@context': 'https://schema.org',
    '@type': 'WebApplication',
    name: 'ATGCFlow',
    description: 'Industry-grade whole exome sequencing analysis platform with GATK best practices, ACMG variant classification, and IGV browser integration.',
    url: process.env.NEXT_PUBLIC_BASE_URL || 'https://atgcflow.com',
    applicationCategory: 'Bioinformatics',
    operatingSystem: 'Web Browser',
    offers: {
      '@type': 'Offer',
      price: '0',
      priceCurrency: 'USD',
    },
    featureList: [
      'Whole Exome Sequencing Analysis',
      'GATK Best Practices Pipeline',
      'ACMG Variant Classification',
      'IGV Browser Integration',
      'Real-time Pipeline Monitoring',
      'Variant Annotation with ANNOVAR',
      'Gene Panel Analysis',
      'Clinical Variant Filtering',
    ],
    provider: {
      '@type': 'Organization',
      name: 'ATGCFlow',
      url: process.env.NEXT_PUBLIC_BASE_URL || 'https://atgcflow.com',
    },
  }

  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        <script
          type="application/ld+json"
          dangerouslySetInnerHTML={{ __html: JSON.stringify(jsonLd) }}
        />
        <link rel="icon" href="/favicon.ico" sizes="any" />
        <link rel="apple-touch-icon" href="/apple-touch-icon.png" />
        <link rel="manifest" href="/manifest.json" />
        <meta name="theme-color" content="#2563eb" />
        <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=5" />
      </head>
      <body className={inter.className}>
        <ThemeProvider
          attribute="class"
          defaultTheme="system"
          enableSystem
          disableTransitionOnChange
        >
          {children}
        </ThemeProvider>
      </body>
    </html>
  )
}
