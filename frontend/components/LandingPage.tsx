'use client';

import { Button } from '@/components/ui/button';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import {
  ArrowRight,
  Dna,
  BarChart3,
  Database,
  Shield,
  Users,
  Microscope,
  FlaskConical,
  Activity,
  Building2,
  Award,
  Globe,
  Zap,
  Target
} from 'lucide-react';

interface LandingPageProps {
  onNavigate: (page: 'home' | 'about' | 'research' | 'contact' | 'signin' | 'features' | 'usecases' | 'publications' | 'pricing') => void;
  onSignIn: () => void;
}

export default function LandingPage({ onNavigate, onSignIn }: LandingPageProps) {
  return (
    <>
      {/* Hero Section - Purple Gradient */}
      <section className="relative bg-gradient-to-br from-purple-primary via-purple-light to-purple-primary overflow-hidden py-24 sm:py-32">
        {/* Subtle floating elements */}
        <div className="absolute inset-0 opacity-10">
          <div className="absolute top-20 left-10 w-72 h-72 bg-cyan rounded-full blur-3xl"></div>
          <div className="absolute bottom-20 right-10 w-96 h-96 bg-accent-teal rounded-full blur-3xl"></div>
        </div>

        <div className="relative max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center space-y-8">
            <Badge variant="cyan" className="w-fit mx-auto text-sm px-4 py-2">
              <Dna className="inline-block w-4 h-4 mr-2" />
              Clinical-Grade Genomics Platform
            </Badge>

            <h1 className="text-4xl sm:text-5xl lg:text-6xl font-bold text-white leading-tight max-w-5xl mx-auto">
              Whole Exome Sequencing Analysis for Research & Clinical Applications
            </h1>

            <p className="text-xl text-blue-100 leading-relaxed max-w-3xl mx-auto">
              Production-grade WES pipeline delivering comprehensive genomic insights with ACMG classification.
              From FASTQ to clinically actionable variants in hours.
            </p>

            <div className="flex flex-col sm:flex-row gap-4 justify-center items-center pt-4">
              <Button
                size="lg"
                className="bg-cyan hover:bg-cyan-light text-white text-lg px-8 py-6 h-auto shadow-xl"
                onClick={onSignIn}
              >
                Start Free Trial
                <ArrowRight className="ml-2 h-5 w-5" />
              </Button>
              <Button
                size="lg"
                variant="outline"
                className="bg-white/10 hover:bg-white/20 text-white border-white/30 text-lg px-8 py-6 h-auto backdrop-blur-sm"
                onClick={() => onNavigate('features')}
              >
                Explore Features
              </Button>
            </div>

            {/* Trust Indicators */}
            <div className="grid grid-cols-3 gap-8 pt-12 max-w-3xl mx-auto">
              <div className="text-center">
                <div className="flex justify-center mb-3">
                  <div className="w-14 h-14 bg-white/10 backdrop-blur-sm rounded-xl flex items-center justify-center">
                    <Shield className="h-7 w-7 text-cyan" />
                  </div>
                </div>
                <div className="text-sm text-blue-100 font-medium">Clinical-Grade Security</div>
              </div>
              <div className="text-center">
                <div className="flex justify-center mb-3">
                  <div className="w-14 h-14 bg-white/10 backdrop-blur-sm rounded-xl flex items-center justify-center">
                    <Award className="h-7 w-7 text-cyan" />
                  </div>
                </div>
                <div className="text-sm text-blue-100 font-medium">Validated Pipeline</div>
              </div>
              <div className="text-center">
                <div className="flex justify-center mb-3">
                  <div className="w-14 h-14 bg-white/10 backdrop-blur-sm rounded-xl flex items-center justify-center">
                    <Globe className="h-7 w-7 text-cyan" />
                  </div>
                </div>
                <div className="text-sm text-blue-100 font-medium">Global Standards</div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Trusted By Section */}
      <section className="py-12 bg-slate-50 border-y border-slate-200">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <p className="text-center text-sm text-slate-600 mb-6 uppercase tracking-wide font-semibold">
            Trusted by Research Institutions Worldwide
          </p>
          <div className="flex flex-wrap justify-center items-center gap-12 opacity-60">
            {/* Placeholder for institution logos */}
            <div className="text-slate-400 text-sm font-medium">Research Labs</div>
            <div className="text-slate-400 text-sm font-medium">Clinical Institutions</div>
            <div className="text-slate-400 text-sm font-medium">Biotech Companies</div>
            <div className="text-slate-400 text-sm font-medium">Academic Centers</div>
          </div>
        </div>
      </section>

      {/* Pipeline Overview Section - Light Background */}
      <section className="py-24 bg-gradient-to-br from-slate-50 via-cyan-50/30 to-purple-50/20">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-purple-primary text-purple-primary">
              8-Stage Pipeline
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">Comprehensive Analysis Pipeline</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              From raw FASTQ to annotated variants with clinical insights
            </p>
          </div>

          <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-6">
            {[
              {
                phase: 'Phase 1',
                title: 'Quality Control',
                duration: '5-10 min',
                description: 'FastP quality filtering and adapter trimming with comprehensive QC metrics.',
                icon: Activity,
                color: 'cyan'
              },
              {
                phase: 'Phase 2',
                title: 'Alignment & Processing',
                duration: '30-60 min',
                description: 'BWA-MEM alignment to GRCh38, sorting, duplicate marking, and BQSR.',
                icon: Database,
                color: 'purple'
              },
              {
                phase: 'Phase 3',
                title: 'Variant Calling',
                duration: '45-90 min',
                description: 'GATK HaplotypeCaller for high-quality variant detection.',
                icon: Dna,
                color: 'teal'
              },
              {
                phase: 'Phase 4',
                title: 'Annotation & Classification',
                duration: '20-40 min',
                description: 'Multi-source annotation with ACMG classification and clinical insights.',
                icon: BarChart3,
                color: 'cyan'
              }
            ].map((item, i) => {
              const IconComponent = item.icon;
              const colorClasses = {
                cyan: 'bg-cyan/10 text-cyan',
                purple: 'bg-purple-primary/10 text-purple-primary',
                teal: 'bg-teal/10 text-teal'
              };
              return (
                <Card key={i} className="relative border-slate-200 bg-white hover:shadow-lg transition-all duration-300 hover:-translate-y-1">
                  <CardContent className="p-8">
                    <div className={`w-14 h-14 rounded-xl flex items-center justify-center mb-4 ${colorClasses[item.color as keyof typeof colorClasses]}`}>
                      <IconComponent className="h-7 w-7" />
                    </div>
                    <div className="text-xs font-semibold text-purple-primary mb-2">{item.phase}</div>
                    <h3 className="text-xl font-bold mb-3 text-slate-900">{item.title}</h3>
                    <p className="text-slate-600 text-sm mb-4 leading-relaxed">{item.description}</p>
                    <div className="text-xs text-slate-500 flex items-center gap-2">
                      <Activity className="h-3 w-3" />
                      {item.duration}
                    </div>
                  </CardContent>
                </Card>
              );
            })}
          </div>

          <div className="text-center mt-12">
            <Button
              variant="outline"
              className="border-purple-primary text-purple-primary hover:bg-purple-primary hover:text-white"
              onClick={() => onNavigate('features')}
            >
              View Full Pipeline Details
              <ArrowRight className="ml-2 h-4 w-4" />
            </Button>
          </div>
        </div>
      </section>

      {/* Features Grid - 2 Column */}
      <section className="py-24 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-cyan text-cyan">
              Platform Capabilities
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">Everything You Need for Clinical Genomics</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              Production-ready tools designed for research and clinical applications
            </p>
          </div>

          <div className="grid md:grid-cols-2 gap-12">
            {[
              {
                icon: FlaskConical,
                title: 'ACMG Classification Engine',
                description: 'Automated variant pathogenicity prediction using ACMG/AMP 2015 guidelines. Comprehensive evidence scoring with PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7 criteria.',
                color: 'cyan'
              },
              {
                icon: Target,
                title: 'Gene Panel Filtering',
                description: 'PanelApp integration for 400+ disease-specific gene panels. Filter variants by ACMG SF v3.2 secondary findings, cancer panels, or custom gene lists.',
                color: 'purple'
              },
              {
                icon: Microscope,
                title: 'IGV Genome Browser',
                description: 'Interactive visualization with IGV.js integration. Real-time variant inspection with BAM/VCF track support, read depth analysis, and allele fraction plots.',
                color: 'teal'
              },
              {
                icon: Shield,
                title: 'Enterprise Security',
                description: 'HIPAA-compliant architecture with Firebase authentication. End-to-end encryption, audit logging, and secure session management for clinical data protection.',
                color: 'cyan'
              }
            ].map((feature, i) => {
              const IconComponent = feature.icon;
              const colorClasses = {
                cyan: 'from-cyan/20 to-cyan/5 border-cyan/20',
                purple: 'from-purple-primary/20 to-purple-primary/5 border-purple-primary/20',
                teal: 'from-teal/20 to-teal/5 border-teal/20'
              };
              const iconColors = {
                cyan: 'text-cyan',
                purple: 'text-purple-primary',
                teal: 'text-teal'
              };
              return (
                <Card key={i} className={`border-2 bg-gradient-to-br ${colorClasses[feature.color as keyof typeof colorClasses]} hover:shadow-xl transition-all duration-300`}>
                  <CardContent className="p-10">
                    <div className={`w-16 h-16 rounded-2xl bg-white shadow-md flex items-center justify-center mb-6 ${iconColors[feature.color as keyof typeof iconColors]}`}>
                      <IconComponent className="h-8 w-8" />
                    </div>
                    <h3 className="text-2xl font-bold text-slate-900 mb-4">{feature.title}</h3>
                    <p className="text-slate-600 leading-relaxed">{feature.description}</p>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </div>
      </section>

      {/* Use Cases Section */}
      <section className="py-24 bg-gradient-to-br from-slate-50 to-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="teal" className="mb-4">
              Who We Serve
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">Built for Every Genomics Workflow</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              From research discovery to clinical diagnostics
            </p>
          </div>

          <div className="grid sm:grid-cols-2 lg:grid-cols-4 gap-8">
            {[
              {
                icon: Microscope,
                title: 'Research Laboratories',
                description: 'Discovery research, biomarker identification, and population studies with scalable analysis.',
                accent: 'cyan'
              },
              {
                icon: Building2,
                title: 'Clinical Diagnostics',
                description: 'Rare disease diagnosis, cancer genomics, and pharmacogenomics with validated pipelines.',
                accent: 'purple'
              },
              {
                icon: FlaskConical,
                title: 'Pharmaceutical R&D',
                description: 'Drug target discovery and clinical trial patient stratification with comprehensive annotation.',
                accent: 'teal'
              },
              {
                icon: Users,
                title: 'Individual Researchers',
                description: 'Academic projects, thesis research, and pilot studies with affordable pricing.',
                accent: 'cyan'
              }
            ].map((useCase, i) => {
              const IconComponent = useCase.icon;
              const accentColors = {
                cyan: 'hover:border-cyan/50 group-hover:text-cyan',
                purple: 'hover:border-purple-primary/50 group-hover:text-purple-primary',
                teal: 'hover:border-teal/50 group-hover:text-teal'
              };
              return (
                <Card
                  key={i}
                  className={`group border-slate-200 hover:shadow-lg transition-all duration-300 cursor-pointer ${accentColors[useCase.accent as keyof typeof accentColors]}`}
                  onClick={() => onNavigate('usecases')}
                >
                  <CardContent className="p-8 text-center">
                    <div className="w-16 h-16 bg-slate-100 rounded-2xl flex items-center justify-center mx-auto mb-6 group-hover:scale-110 transition-transform duration-300">
                      <IconComponent className={`h-8 w-8 text-slate-600 ${accentColors[useCase.accent as keyof typeof accentColors]}`} />
                    </div>
                    <h3 className="text-lg font-bold text-slate-900 mb-3">{useCase.title}</h3>
                    <p className="text-sm text-slate-600 leading-relaxed">{useCase.description}</p>
                  </CardContent>
                </Card>
              );
            })}
          </div>

          <div className="text-center mt-12">
            <Button
              onClick={() => onNavigate('usecases')}
              className="bg-purple-primary hover:bg-purple-light text-white"
            >
              Explore Use Cases
              <ArrowRight className="ml-2 h-4 w-4" />
            </Button>
          </div>
        </div>
      </section>

      {/* Metrics Section - Trust Indicators */}
      <section className="py-24 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <h2 className="text-4xl font-bold text-slate-900 mb-4">Platform Performance</h2>
            <p className="text-xl text-slate-600">Validated metrics on clinical-grade datasets</p>
          </div>

          <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-8">
            {[
              { label: 'Coverage Depth', value: '100x+', desc: 'Mean target coverage', icon: BarChart3 },
              { label: 'Processing Time', value: '2-4h', desc: 'FASTQ to final VCF', icon: Zap },
              { label: 'Variant Calls', value: '~80K', desc: 'Per exome sample', icon: Database },
              { label: 'Annotations', value: '20+', desc: 'Functional predictors', icon: Activity }
            ].map((metric, i) => {
              const IconComponent = metric.icon;
              return (
                <Card key={i} className="border-slate-200 bg-gradient-to-br from-white to-slate-50 hover:shadow-lg transition-all duration-300">
                  <CardContent className="p-8 text-center">
                    <IconComponent className="h-10 w-10 text-purple-primary mx-auto mb-4" />
                    <div className="text-4xl font-bold text-purple-primary mb-2">{metric.value}</div>
                    <div className="font-semibold text-slate-900 mb-1">{metric.label}</div>
                    <div className="text-sm text-slate-600">{metric.desc}</div>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </div>
      </section>

      {/* CTA Section - Purple Gradient */}
      <section className="relative py-24 bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark overflow-hidden">
        <div className="absolute inset-0 opacity-10">
          <div className="absolute top-10 right-20 w-96 h-96 bg-cyan rounded-full blur-3xl"></div>
        </div>

        <div className="relative max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-4xl sm:text-5xl font-bold text-white mb-6">
            Ready to Transform Your Genomic Research?
          </h2>
          <p className="text-xl text-blue-100 mb-10 max-w-2xl mx-auto leading-relaxed">
            Join leading research institutions using our clinical-grade WES pipeline.
            Start analyzing your exome data today with our free trial.
          </p>

          <div className="flex flex-col sm:flex-row gap-4 justify-center items-center">
            <Button
              size="lg"
              className="bg-cyan hover:bg-cyan-light text-white text-lg px-8 py-6 h-auto shadow-xl"
              onClick={onSignIn}
            >
              Start Free Trial
              <ArrowRight className="ml-2 h-5 w-5" />
            </Button>
            <Button
              size="lg"
              variant="outline"
              className="bg-white/10 hover:bg-white/20 text-white border-white/30 text-lg px-8 py-6 h-auto backdrop-blur-sm"
              onClick={() => onNavigate('pricing')}
            >
              View Pricing
            </Button>
          </div>

          <p className="text-sm text-blue-200 mt-8">
            No credit card required • 2 free jobs per month • Cancel anytime
          </p>
        </div>
      </section>
    </>
  );
}
