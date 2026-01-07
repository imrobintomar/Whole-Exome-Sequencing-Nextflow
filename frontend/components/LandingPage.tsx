'use client';

import { useEffect, useState } from 'react';
import { Button } from '@/components/ui/button';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import {
  ArrowRight,
  Dna,
  BarChart3,
  Database,
  Globe,
  Shield,
  Users,
  Clock,
  ChevronRight,
  Microscope,
  FlaskConical,
  Activity
} from 'lucide-react';

interface LandingPageProps {
  onNavigate: (page: string) => void;
  onSignIn: () => void;
}

export default function LandingPage({ onNavigate, onSignIn }: LandingPageProps) {
  const [jobCount, setJobCount] = useState(0);
  const [sampleCount, setSampleCount] = useState(0);
  const [variantCount, setVariantCount] = useState(0);

  useEffect(() => {
    const jobInterval = setInterval(() => {
      setJobCount(prev => prev < 156 ? prev + 2 : 156);
    }, 20);
    const sampleInterval = setInterval(() => {
      setSampleCount(prev => prev < 423 ? prev + 5 : 423);
    }, 20);
    const variantInterval = setInterval(() => {
      setVariantCount(prev => prev < 12547893 ? prev + 125478 : 12547893);
    }, 20);

    return () => {
      clearInterval(jobInterval);
      clearInterval(sampleInterval);
      clearInterval(variantInterval);
    };
  }, []);

  return (
    <>
      {/* Hero Section */}
      <section className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-16 sm:py-24">
        <div className="grid lg:grid-cols-2 gap-12 items-center">
          <div className="space-y-8">
            <div className="space-y-4">
              <Badge variant="outline" className="w-fit border-blue-200 bg-blue-50 text-blue-700">
                <span className="inline-block w-2 h-2 bg-blue-600 rounded-full mr-2"></span>
                Nextflow Pipeline Platform
              </Badge>

              <h1 className="text-4xl sm:text-5xl font-bold text-slate-900 leading-tight">
                Whole Exome Sequencing Analysis Pipeline
              </h1>

              <p className="text-lg text-slate-600 leading-relaxed">
                Production-grade WES pipeline powered by Nextflow. Process FASTQ files through BWA alignment, GATK variant calling, VEP annotation, and deliver comprehensive genomic insights with ACMG classification.
              </p>
            </div>

            <div className="flex flex-col sm:flex-row gap-4">
              <Button size="lg" className="bg-blue-600 hover:bg-blue-700 text-white" onClick={onSignIn}>
                Get Started
                <ArrowRight className="ml-2 h-4 w-4" />
              </Button>
              <Button size="lg" variant="outline" onClick={() => onNavigate('about')}>
                Learn More
              </Button>
            </div>

            {/* Stats */}
            <div className="grid grid-cols-3 gap-4 pt-4">
              <div>
                <div className="text-2xl font-bold text-slate-900">{jobCount.toLocaleString()}</div>
                <div className="text-sm text-slate-600">Jobs Processed</div>
              </div>
              <div>
                <div className="text-2xl font-bold text-slate-900">{sampleCount.toLocaleString()}</div>
                <div className="text-sm text-slate-600">Samples Analyzed</div>
              </div>
              <div>
                <div className="text-2xl font-bold text-slate-900">{(variantCount / 1000000).toFixed(1)}M</div>
                <div className="text-sm text-slate-600">Variants Called</div>
              </div>
            </div>
          </div>

          {/* Feature Highlight Card */}
          <div className="relative">
            <div className="absolute inset-0 bg-gradient-to-r from-blue-100 to-cyan-100 rounded-2xl blur-3xl opacity-30"></div>
            <Card className="relative border-slate-200 bg-white shadow-lg">
              <CardContent className="p-8">
                <div className="space-y-6">
                  <div className="flex items-start gap-4">
                    <div className="w-12 h-12 bg-blue-100 rounded-lg flex items-center justify-center flex-shrink-0">
                      <Database className="h-6 w-6 text-blue-600" />
                    </div>
                    <div>
                      <h3 className="font-semibold text-slate-900">Complete Pipeline</h3>
                      <p className="text-sm text-slate-600 mt-1">FastP QC → BWA alignment → GATK BQSR → Variant calling → VEP annotation</p>
                    </div>
                  </div>

                  <div className="flex items-start gap-4">
                    <div className="w-12 h-12 bg-cyan-100 rounded-lg flex items-center justify-center flex-shrink-0">
                      <FlaskConical className="h-6 w-6 text-cyan-600" />
                    </div>
                    <div>
                      <h3 className="font-semibold text-slate-900">ACMG Classification</h3>
                      <p className="text-sm text-slate-600 mt-1">Automated variant classification using ACMG/AMP 2015 guidelines with evidence scoring</p>
                    </div>
                  </div>

                  <div className="flex items-start gap-4">
                    <div className="w-12 h-12 bg-purple-100 rounded-lg flex items-center justify-center flex-shrink-0">
                      <Microscope className="h-6 w-6 text-purple-600" />
                    </div>
                    <div>
                      <h3 className="font-semibold text-slate-900">IGV Integration</h3>
                      <p className="text-sm text-slate-600 mt-1">Interactive genome browser for variant visualization and validation</p>
                    </div>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </section>

      {/* Pipeline Section */}
      <section className="bg-slate-900 text-white py-16 sm:py-24">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <h2 className="text-3xl sm:text-4xl font-bold mb-4">Multi-Step Analysis Pipeline</h2>
            <p className="text-lg text-slate-400">From raw FASTQ to annotated variants with clinical insights</p>
          </div>

          <div className="grid md:grid-cols-4 gap-6">
            {[
              {
                phase: 'Phase 1',
                title: 'Quality Control',
                duration: '~5-10 min',
                description: 'FastP quality filtering and adapter trimming with comprehensive QC metrics.',
                icon: Activity
              },
              {
                phase: 'Phase 2',
                title: 'Alignment & Processing',
                duration: '30-60 min',
                description: 'BWA-MEM alignment to GRCh38, sorting, duplicate marking, and base quality recalibration.',
                icon: Database
              },
              {
                phase: 'Phase 3',
                title: 'Variant Calling',
                duration: '45-90 min',
                description: 'GATK HaplotypeCaller for high-quality variant detection with confidence scoring.',
                icon: Dna
              },
              {
                phase: 'Phase 4',
                title: 'Annotation & Report',
                duration: '20-40 min',
                description: 'VEP annotation, ANNOVAR integration, 1000 Genomes frequencies, and ACMG classification.',
                icon: BarChart3
              }
            ].map((item, i) => {
              const IconComponent = item.icon;
              return (
                <div key={i} className="relative">
                  {i < 3 && (
                    <div className="hidden md:block absolute top-12 -right-3 text-slate-700 z-10">
                      <ChevronRight className="h-6 w-6" />
                    </div>
                  )}
                  <Card className="border-slate-700 bg-slate-800/50 h-full">
                    <CardContent className="p-6">
                      <div className="flex items-center gap-3 mb-4">
                        <div className="w-10 h-10 bg-blue-500/20 rounded-lg flex items-center justify-center">
                          <IconComponent className="h-5 w-5 text-blue-400" />
                        </div>
                        <span className="text-sm font-semibold text-blue-400">{item.phase}</span>
                      </div>
                      <h3 className="text-xl font-bold mb-2">{item.title}</h3>
                      <p className="text-slate-400 text-sm mb-4">{item.description}</p>
                      <div className="text-sm text-slate-500 flex items-center gap-2">
                        <Clock className="h-4 w-4" />
                        {item.duration}
                      </div>
                    </CardContent>
                  </Card>
                </div>
              );
            })}
          </div>
        </div>
      </section>

      {/* Features Grid */}
      <section className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-16 sm:py-24">
        <div className="text-center mb-16">
          <h2 className="text-3xl sm:text-4xl font-bold text-slate-900 mb-4">Platform Capabilities</h2>
          <p className="text-lg text-slate-600">Everything you need for production exome analysis</p>
        </div>

        <div className="grid md:grid-cols-2 gap-8">
          {[
            {
              icon: Globe,
              title: 'Gene Panel Filtering',
              description: 'PanelApp integration for disease-specific gene panels. Filter variants by ACMG secondary findings or custom gene lists.',
              color: 'cyan' as const
            },
            {
              icon: Shield,
              title: 'Enterprise Security',
              description: 'Firebase Authentication with secure session management. HTTPS encryption for all data transfers and API calls.',
              color: 'blue' as const
            },
            {
              icon: Users,
              title: 'Multi-User Platform',
              description: 'Per-user job tracking with SQLite database. Persistent sessions and job history for collaborative workflows.',
              color: 'purple' as const
            },
            {
              icon: Microscope,
              title: 'Interactive Visualization',
              description: 'IGV.js genome browser integration. Real-time variant visualization with BAM/VCF track support.',
              color: 'yellow' as const
            }
          ].map((feature, i) => {
            const IconComponent = feature.icon;
            const colorClasses: Record<'cyan' | 'blue' | 'purple' | 'yellow', string> = {
              cyan: 'bg-cyan-100 text-cyan-600',
              blue: 'bg-blue-100 text-blue-600',
              purple: 'bg-purple-100 text-purple-600',
              yellow: 'bg-yellow-100 text-yellow-600'
            };
            return (
              <Card key={i} className="border-slate-200">
                <CardContent className="p-8">
                  <div className={`w-12 h-12 rounded-lg flex items-center justify-center mb-4 ${colorClasses[feature.color]}`}>
                    <IconComponent className="h-6 w-6" />
                  </div>
                  <h3 className="text-xl font-bold text-slate-900 mb-2">{feature.title}</h3>
                  <p className="text-slate-600">{feature.description}</p>
                </CardContent>
              </Card>
            );
          })}
        </div>
      </section>

      {/* Metrics */}
      <section className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-16 sm:py-24">
        <div className="text-center mb-16">
          <h2 className="text-3xl sm:text-4xl font-bold text-slate-900 mb-4">Pipeline Metrics</h2>
          <p className="text-lg text-slate-600">Performance benchmarks on typical WES datasets</p>
        </div>

        <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-6">
          {[
            { label: 'Coverage', value: '100x', desc: 'Mean target coverage' },
            { label: 'Processing Time', value: '2-4h', desc: 'FASTQ to final VCF' },
            { label: 'Variant Calls', value: '~80K', desc: 'Per exome sample' },
            { label: 'Annotations', value: '20+', desc: 'Functional predictors' }
          ].map((metric, i) => (
            <Card key={i} className="border-slate-200 bg-gradient-to-br from-slate-50 to-white">
              <CardContent className="p-6">
                <div className="text-3xl font-bold text-blue-600 mb-2">{metric.value}</div>
                <div className="font-semibold text-slate-900 mb-1">{metric.label}</div>
                <div className="text-sm text-slate-600">{metric.desc}</div>
              </CardContent>
            </Card>
          ))}
        </div>
      </section>

      {/* CTA */}
      <section className="bg-gradient-to-r from-blue-600 to-cyan-600 text-white py-16 sm:py-24">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl sm:text-4xl font-bold mb-4">Ready to Analyze Your Exome Data?</h2>
          <p className="text-lg text-blue-100 mb-8 max-w-2xl mx-auto">
            Start processing your whole-exome sequencing data with our production-ready Nextflow pipeline and comprehensive variant annotation.
          </p>

          <Button size="lg" className="bg-white text-blue-600 hover:bg-blue-50 font-semibold" onClick={onSignIn}>
            Get Started Free
            <ArrowRight className="ml-2 h-4 w-4" />
          </Button>
        </div>
      </section>
    </>
  );
}
