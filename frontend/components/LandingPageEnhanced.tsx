'use client';

import { useState, useEffect } from 'react';
import Image from 'next/image';
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
  Target,
  ChevronRight,
  CheckCircle2,
  Play,
  Clock,
  TrendingUp
} from 'lucide-react';

interface LandingPageProps {
  onNavigate: (page: 'home' | 'about' | 'research' | 'contact' | 'signin' | 'features' | 'usecases' | 'publications' | 'pricing') => void;
  onSignIn: () => void;
}

export default function LandingPageEnhanced({ onNavigate, onSignIn }: LandingPageProps) {
  const [activeSolution, setActiveSolution] = useState(0);
  const [isVisible, setIsVisible] = useState(false);

  // Fade-in animation on mount
  useEffect(() => {
    setIsVisible(true);
  }, []);

  // Auto-rotate solutions every 5 seconds
  useEffect(() => {
    const interval = setInterval(() => {
      setActiveSolution((prev) => (prev + 1) % 3);
    }, 5000);
    return () => clearInterval(interval);
  }, []);

  const solutions = [
    {
      title: 'Research & Discovery',
      subtitle: 'Accelerate Genomic Discovery',
      description: 'Comprehensive WES analysis for biomarker identification, population studies, and novel variant discovery. Scalable pipeline supporting batch processing of hundreds of samples with standardized quality metrics.',
      features: [
        'Batch processing up to 500 samples',
        'Population-scale variant analysis',
        'Custom gene panel creation',
        'Export ready for publication'
      ],
      icon: Microscope,
      image: '/AccelerateGenomicDiscovery.jpg',
      color: 'cyan'
    },
    {
      title: 'Clinical Diagnostics',
      subtitle: 'Clinical-Grade Variant Analysis',
      description: 'Production-ready pipeline for rare disease diagnosis, cancer genomics, and pharmacogenomics. ACMG classification engine with automated pathogenicity scoring and secondary findings reporting.',
      features: [
        'ACMG/AMP 2015 guidelines',
        'Automated clinical reports',
        'HIPAA-compliant infrastructure',
        '2-4 hour turnaround time'
      ],
      icon: Building2,
      image: '/Clinical-GradeVariantAnalysis.jpg',
      color: 'purple'
    },
    {
      title: 'Pharmaceutical R&D',
      subtitle: 'Drug Development Support',
      description: 'Streamlined variant analysis for clinical trials, patient stratification, and biomarker validation. API integration for seamless LIMS connectivity and multi-site collaboration.',
      features: [
        'Clinical trial patient selection',
        'Biomarker discovery pipeline',
        'RESTful API integration',
        'Multi-center data aggregation'
      ],
      icon: FlaskConical,
      image: '/DrugDevelopmentSupport.jpg',
      color: 'teal'
    }
  ];

  return (
    <>
      {/* ENHANCED Hero Section */}
      <section className="relative bg-gradient-to-br from-[#060140] via-[#0a0560] to-[#060140] overflow-hidden py-16 sm:py-24 lg:py-32">
        {/* Animated background elements */}
        <div className="absolute inset-0 opacity-20">
          <div className="absolute top-20 left-10 w-96 h-96 bg-cyan rounded-full blur-3xl animate-pulse"></div>
          <div className="absolute bottom-20 right-10 w-[32rem] h-[32rem] bg-purple-500 rounded-full blur-3xl animate-pulse" style={{ animationDelay: '1s' }}></div>
        </div>

        {/* Geometric pattern overlay */}
        <div className="absolute inset-0 opacity-5">
          <svg className="w-full h-full" xmlns="http://www.w3.org/2000/svg">
            <defs>
              <pattern id="grid" width="40" height="40" patternUnits="userSpaceOnUse">
                <path d="M 40 0 L 0 0 0 40" fill="none" stroke="white" strokeWidth="0.5"/>
              </pattern>
            </defs>
            <rect width="100%" height="100%" fill="url(#grid)" />
          </svg>
        </div>

        {/* DNA Helix Decorative Element */}
        <div className="absolute right-0 top-1/2 -translate-y-1/2 opacity-10 hidden lg:block">
          <Dna className="w-[30rem] h-[30rem] text-cyan animate-pulse" />
        </div>

        <div className={`relative max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 transition-all duration-1000 ${isVisible ? 'opacity-100 translate-y-0' : 'opacity-0 translate-y-10'}`}>
          <div className="grid lg:grid-cols-2 gap-12 items-center">
            {/* Left Column - Content */}
            <div className="space-y-8">
              {/* Badge with Social Proof */}
              <div className="inline-flex items-center gap-2 px-4 py-2 rounded-full bg-cyan/20 backdrop-blur-sm border border-cyan/30 animate-fade-in-down">
                <Award className="w-4 h-4 text-cyan" />
                <span className="text-sm font-medium text-white">Trusted by 500+ Researchers Worldwide</span>
              </div>

              {/* Main Headline - More Impactful */}
              <h1 className="text-4xl sm:text-5xl lg:text-6xl xl:text-7xl font-bold text-white leading-tight animate-fade-in-up">
                <span className="inline-block animate-word-reveal">Transform Raw</span>{' '}
                <span className="inline-block animate-word-reveal" style={{ animationDelay: '0.1s' }}>Genomic Data</span>
                <span className="block text-transparent bg-clip-text bg-gradient-to-r from-cyan via-teal to-cyan mt-2 animate-word-reveal" style={{ animationDelay: '0.2s' }}>
                  Into Clinical Insights
                </span>
                <span className="block text-2xl sm:text-3xl lg:text-4xl text-blue-200 mt-4 font-light animate-word-reveal" style={{ animationDelay: '0.3s' }}>
                  in Hours, Not Days
                </span>
              </h1>

              {/* Subheadline - Simplified & Benefit-Focused */}
              <p className="text-lg sm:text-xl text-blue-100 leading-relaxed max-w-2xl font-light animate-fade-in" style={{ animationDelay: '0.4s' }}>
                Production-grade WES pipeline with automated variant classification.
                <span className="block mt-2 text-cyan font-medium">From sequencing data to actionable reports—automatically.</span>
              </p>

              {/* Trust Indicators - Above Fold */}
              <div className="flex flex-wrap gap-6 text-sm text-blue-200 animate-fade-in" style={{ animationDelay: '0.5s' }}>
                <div className="flex items-center gap-2">
                  <CheckCircle2 className="w-5 h-5 text-cyan" />
                  <span>99.5% Sensitivity</span>
                </div>
                <div className="flex items-center gap-2">
                  <CheckCircle2 className="w-5 h-5 text-cyan" />
                  <span>ACMG Validated</span>
                </div>
                <div className="flex items-center gap-2">
                  <CheckCircle2 className="w-5 h-5 text-cyan" />
                  <span>HIPAA Compliant</span>
                </div>
              </div>

              {/* Enhanced CTAs */}
              <div className="flex flex-col sm:flex-row gap-4 pt-4 animate-fade-in" style={{ animationDelay: '0.6s' }}>
                <Button
                  size="lg"
                  className="bg-cyan hover:bg-cyan-light text-white text-lg px-10 py-7 h-auto shadow-2xl hover:shadow-cyan/50 transition-all duration-300 hover:scale-105 hover:-translate-y-1 gooey-button group"
                  onClick={onSignIn}
                >
                  Start Free Trial
                  <ArrowRight className="ml-2 h-5 w-5 group-hover:translate-x-1 transition-transform" />
                </Button>
                <Button
                  size="lg"
                  variant="outline"
                  className="bg-transparent hover:bg-white/10 text-white border-2 border-white/40 hover:border-white/60 text-lg px-10 py-7 h-auto backdrop-blur-sm transition-all duration-300 hover:scale-105 group"
                  onClick={() => onNavigate('features')}
                >
                  <Play className="mr-2 h-5 w-5 group-hover:scale-110 transition-transform" />
                  See Demo
                </Button>
              </div>

              {/* Micro-conversion - Low commitment */}
              <p className="text-sm text-blue-300 animate-fade-in" style={{ animationDelay: '0.7s' }}>
                🎯 No credit card required • 2 free jobs per month • Setup in 5 minutes
              </p>
            </div>

            {/* Right Column - Visual Element */}
            <div className="hidden lg:block relative animate-fade-in" style={{ animationDelay: '0.4s' }}>
              <div className="relative">
                {/* Floating Stats Cards */}
                <div className="absolute -top-6 -left-6 z-20">
                  <Card className="bg-white/95 backdrop-blur-sm border-2 border-cyan/20 shadow-2xl animate-float" style={{ animationDelay: '0s' }}>
                    <CardContent className="p-4">
                      <div className="flex items-center gap-3">
                        <div className="w-12 h-12 bg-cyan/10 rounded-xl flex items-center justify-center">
                          <Clock className="w-6 h-6 text-cyan" />
                        </div>
                        <div>
                          <div className="text-2xl font-bold text-[#060140]">2-4h</div>
                          <div className="text-xs text-slate-600">Turnaround Time</div>
                        </div>
                      </div>
                    </CardContent>
                  </Card>
                </div>

                <div className="absolute -bottom-6 -right-6 z-20">
                  <Card className="bg-white/95 backdrop-blur-sm border-2 border-teal/20 shadow-2xl animate-float" style={{ animationDelay: '0.5s' }}>
                    <CardContent className="p-4">
                      <div className="flex items-center gap-3">
                        <div className="w-12 h-12 bg-teal/10 rounded-xl flex items-center justify-center">
                          <TrendingUp className="w-6 h-6 text-teal" />
                        </div>
                        <div>
                          <div className="text-2xl font-bold text-[#060140]">99.5%</div>
                          <div className="text-xs text-slate-600">Sensitivity</div>
                        </div>
                      </div>
                    </CardContent>
                  </Card>
                </div>

                {/* Central Visual - Pipeline Diagram or Screenshot Placeholder */}
                <div className="aspect-square bg-gradient-to-br from-cyan/20 to-purple-500/20 rounded-3xl backdrop-blur-sm border-2 border-white/20 p-8 shadow-2xl">
                  <div className="w-full h-full bg-gradient-to-br from-white/10 to-white/5 rounded-2xl flex flex-col items-center justify-center gap-6 p-6">
                    {/* Pipeline Steps Visualization */}
                    {[
                      { icon: Database, label: 'Upload FASTQ', color: 'text-cyan' },
                      { icon: Activity, label: 'QC & Alignment', color: 'text-teal' },
                      { icon: Dna, label: 'Variant Calling', color: 'text-purple-300' },
                      { icon: BarChart3, label: 'ACMG Classification', color: 'text-cyan' }
                    ].map((step, idx) => {
                      const IconComponent = step.icon;
                      return (
                        <div key={idx} className="flex items-center gap-4 w-full">
                          <div className="flex items-center justify-center w-12 h-12 bg-white/20 rounded-xl backdrop-blur-sm">
                            <IconComponent className={`w-6 h-6 ${step.color}`} />
                          </div>
                          <div className="flex-1">
                            <div className="h-2 bg-white/10 rounded-full overflow-hidden">
                              <div
                                className="h-full bg-gradient-to-r from-cyan to-teal rounded-full animate-pulse"
                                style={{ width: `${100 - (idx * 15)}%` }}
                              />
                            </div>
                          </div>
                          <span className="text-xs font-medium text-white/80">{step.label}</span>
                        </div>
                      );
                    })}
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* Custom Animations */}
        <style jsx>{`
          @keyframes float {
            0%, 100% {
              transform: translateY(0px);
            }
            50% {
              transform: translateY(-10px);
            }
          }

          .animate-float {
            animation: float 3s ease-in-out infinite;
          }
        `}</style>
      </section>

      {/* What We Do Section with Slide-in Animation */}
      <section className="py-20 bg-white overflow-hidden">
        <div className="max-w-6xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl sm:text-4xl lg:text-5xl font-bold text-slate-900 mb-6 slide-in-left">
            What We Do
          </h2>
          <p className="text-xl text-slate-600 leading-relaxed max-w-4xl mx-auto slide-in-right">
            ATGC Flow provides a <strong className="text-[#060140]">production-grade Whole Exome Sequencing analysis platform</strong> that transforms raw genomic data into clinically actionable insights. Our automated pipeline combines best-practice bioinformatics tools with advanced annotation databases and machine learning-based variant classification, enabling researchers and clinicians to accelerate discovery and improve diagnostic accuracy.
          </p>
        </div>
      </section>

      {/* Solutions Section - Tab-Based with Smooth Transitions */}
      <section className="py-24 bg-gradient-to-br from-slate-50 via-white to-slate-50">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-[#060140] text-[#060140] text-sm px-4 py-2">
              Solutions
            </Badge>
            <h2 className="text-4xl sm:text-5xl font-bold text-slate-900 mb-4">
              Tailored for Your Genomics Workflow
            </h2>
            <p className="text-xl text-slate-600 max-w-3xl mx-auto">
              Whether you're advancing research, diagnosing rare diseases, or developing therapeutics,
              our platform adapts to your specific needs
            </p>
          </div>

          {/* Tab Navigation */}
          <div className="flex flex-col lg:flex-row gap-4 mb-8 justify-center">
            {solutions.map((solution, index) => {
              const IconComponent = solution.icon;
              return (
                <button
                  key={index}
                  onClick={() => setActiveSolution(index)}
                  className={`flex items-center gap-3 px-6 py-4 rounded-xl transition-all duration-500 ${
                    activeSolution === index
                      ? 'bg-[#060140] text-white shadow-2xl scale-105 -translate-y-1'
                      : 'bg-white text-slate-700 hover:bg-slate-50 border-2 border-slate-200 hover:border-[#060140] hover:shadow-lg'
                  }`}
                >
                  <IconComponent className="h-6 w-6" />
                  <span className="font-bold text-lg">{solution.title}</span>
                </button>
              );
            })}
          </div>

          {/* Tab Content with Fade & Slide Animation */}
          <div className="relative min-h-[500px]">
            {solutions.map((solution, index) => {
              const IconComponent = solution.icon;
              const colorClasses = {
                cyan: 'from-cyan/10 to-transparent',
                purple: 'from-purple-500/10 to-transparent',
                teal: 'from-teal/10 to-transparent'
              };
              return (
                <div
                  key={index}
                  className={`absolute inset-0 transition-all duration-700 ${
                    activeSolution === index
                      ? 'opacity-100 translate-x-0 z-10'
                      : activeSolution > index
                      ? 'opacity-0 -translate-x-full z-0'
                      : 'opacity-0 translate-x-full z-0'
                  }`}
                >
                  <Card className={`border-2 border-slate-200 bg-gradient-to-br ${colorClasses[solution.color as keyof typeof colorClasses]} hover:shadow-2xl transition-shadow duration-500`}>
                    <CardContent className="p-12">
                      <div className="grid lg:grid-cols-2 gap-12 items-center">
                        <div>
                          <div className="w-20 h-20 bg-[#060140] rounded-2xl flex items-center justify-center mb-6 hover:rotate-12 transition-transform duration-600">
                            <IconComponent className="h-10 w-10 text-cyan" />
                          </div>
                          <h3 className="text-3xl font-bold text-slate-900 mb-2">{solution.subtitle}</h3>
                          <p className="text-lg text-slate-600 leading-relaxed mb-8">
                            {solution.description}
                          </p>
                          <div className="space-y-3">
                            {solution.features.map((feature, idx) => (
                              <div
                                key={idx}
                                className="flex items-start gap-3 opacity-0 animate-slide-in-feature"
                                style={{ animationDelay: `${idx * 0.1}s`, animationFillMode: 'forwards' }}
                              >
                                <CheckCircle2 className="h-5 w-5 text-cyan flex-shrink-0 mt-1" />
                                <span className="text-slate-700">{feature}</span>
                              </div>
                            ))}
                          </div>
                          <Button
                            className="mt-8 bg-[#060140] hover:bg-[#0a0560] text-white transition-all duration-300 hover:shadow-lg hover:scale-105"
                            onClick={() => onNavigate('usecases')}
                          >
                            Learn More
                            <ChevronRight className="ml-2 h-4 w-4" />
                          </Button>
                        </div>
                        <div className="relative group">
                          <div className="aspect-square bg-gradient-to-br from-slate-100 to-slate-200 rounded-2xl overflow-hidden shadow-xl">
                            <div className="relative w-full h-full group-hover:rotate-6 transition-transform duration-600 ease-out">
                              <Image
                                src={solution.image}
                                alt={solution.subtitle}
                                fill
                                className="object-cover group-hover:scale-110 transition-transform duration-600"
                                priority={index === 0}
                                unoptimized
                              />
                              {/* Overlay gradient on hover */}
                              <div className="absolute inset-0 bg-gradient-to-t from-[#060140]/80 via-[#060140]/20 to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-600 flex items-end p-8">
                                <div className="text-white">
                                  <IconComponent className="h-12 w-12 mb-3" />
                                  <p className="text-sm font-semibold">{solution.title}</p>
                                </div>
                              </div>
                            </div>
                          </div>
                          {/* Decorative elements */}
                          <div className="absolute -top-4 -right-4 w-24 h-24 bg-cyan/20 rounded-full blur-2xl group-hover:scale-150 transition-transform duration-600"></div>
                          <div className="absolute -bottom-4 -left-4 w-32 h-32 bg-[#060140]/20 rounded-full blur-2xl group-hover:scale-150 transition-transform duration-600"></div>
                        </div>
                      </div>
                    </CardContent>
                  </Card>
                </div>
              );
            })}
          </div>

          {/* Pagination Dots */}
          <div className="flex justify-center gap-2 mt-8">
            {solutions.map((_, index) => (
              <button
                key={index}
                onClick={() => setActiveSolution(index)}
                className={`h-2 rounded-full transition-all duration-500 ${
                  activeSolution === index
                    ? 'w-8 bg-[#060140]'
                    : 'w-2 bg-slate-300 hover:bg-slate-400'
                }`}
                aria-label={`Go to solution ${index + 1}`}
              />
            ))}
          </div>
        </div>
      </section>

      {/* Company Values Section with Staggered Animation */}
      <section className="py-24 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <h2 className="text-4xl sm:text-5xl font-bold text-slate-900 mb-4">
              Our Core Values
            </h2>
            <p className="text-xl text-slate-600">What drives us to deliver excellence</p>
          </div>

          <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-8">
            {[
              {
                icon: Shield,
                title: 'Quality First',
                description: 'Clinical-grade quality is our core tenet. Every variant call undergoes rigorous validation against benchmark datasets.',
                shape: 'square'
              },
              {
                icon: Zap,
                title: 'Innovation',
                description: 'Constantly evolving our pipeline with cutting-edge algorithms and the latest genomic databases.',
                shape: 'circle'
              },
              {
                icon: Users,
                title: 'Collaboration',
                description: 'Built with feedback from leading researchers, clinicians, and bioinformaticians worldwide.',
                shape: 'triangle'
              },
              {
                icon: Globe,
                title: 'Global Impact',
                description: 'Democratizing access to clinical-grade genomics for institutions and researchers everywhere.',
                shape: 'hexagon'
              }
            ].map((value, index) => {
              const IconComponent = value.icon;
              return (
                <Card
                  key={index}
                  className="border-2 border-slate-200 hover:border-[#060140] hover:shadow-xl transition-all duration-500 group value-card"
                  style={{ animationDelay: `${index * 0.15}s` }}
                >
                  <CardContent className="p-8 text-center">
                    <div className="relative w-20 h-20 mx-auto mb-6 perspective-1000">
                      <div className="absolute inset-0 bg-gradient-to-br from-[#060140] to-cyan rounded-2xl group-hover:rotate-180 transition-transform duration-600 ease-out"></div>
                      <div className="absolute inset-0 flex items-center justify-center">
                        <IconComponent className="h-10 w-10 text-white relative z-10" />
                      </div>
                    </div>
                    <h3 className="text-xl font-bold text-slate-900 mb-3">{value.title}</h3>
                    <p className="text-slate-600 leading-relaxed text-sm">{value.description}</p>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </div>
      </section>

      {/* Trust Indicators with Fade-in */}
      <section className="py-16 bg-gradient-to-br from-slate-50 to-white border-y border-slate-200">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <p className="text-center text-sm text-slate-600 mb-8 uppercase tracking-wider font-semibold">
            Trusted by Leading Institutions
          </p>
          <div className="flex flex-wrap justify-center items-center gap-12 lg:gap-16">
            {[
              'Research Universities',
              'Clinical Labs',
              'Biotech Companies',
              'Pharmaceutical R&D',
              'Academic Medical Centers'
            ].map((name, i) => (
              <div
                key={i}
                className="text-slate-400 font-bold text-sm hover:text-[#060140] transition-all duration-500 cursor-pointer hover:scale-110 opacity-0 animate-fade-in-up"
                style={{ animationDelay: `${i * 0.1}s`, animationFillMode: 'forwards' }}
              >
                {name}
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Platform Metrics with Counter Animation */}
      <section className="py-24 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <h2 className="text-4xl sm:text-5xl font-bold text-slate-900 mb-4">
              Proven Performance
            </h2>
            <p className="text-xl text-slate-600">Validated on clinical and research datasets</p>
          </div>

          <div className="grid sm:grid-cols-2 lg:grid-cols-4 gap-8">
            {[
              {
                value: '99.5%',
                label: 'Sensitivity',
                desc: 'Validated on GIAB samples',
                icon: Target,
                color: 'cyan'
              },
              {
                value: '2-4h',
                label: 'Turnaround',
                desc: 'FASTQ to final report',
                icon: Zap,
                color: 'purple'
              },
              {
                value: '100x+',
                label: 'Coverage',
                desc: 'Mean target depth',
                icon: BarChart3,
                color: 'teal'
              },
              {
                value: '20+',
                label: 'Databases',
                desc: 'Annotation sources',
                icon: Database,
                color: 'cyan'
              }
            ].map((metric, i) => {
              const IconComponent = metric.icon;
              const colorClasses = {
                cyan: 'text-cyan',
                purple: 'text-[#060140]',
                teal: 'text-teal'
              };
              return (
                <div
                  key={i}
                  className="text-center group metric-card"
                  style={{ animationDelay: `${i * 0.1}s` }}
                >
                  <div className="w-16 h-16 bg-slate-100 rounded-2xl flex items-center justify-center mx-auto mb-4 group-hover:bg-[#060140] transition-all duration-500 group-hover:rotate-360">
                    <IconComponent className={`h-8 w-8 ${colorClasses[metric.color as keyof typeof colorClasses]} group-hover:text-white transition-colors duration-500`} />
                  </div>
                  <div className="text-5xl font-bold text-[#060140] mb-2 group-hover:scale-110 transition-transform duration-500">{metric.value}</div>
                  <div className="text-lg font-semibold text-slate-900 mb-1">{metric.label}</div>
                  <div className="text-sm text-slate-600">{metric.desc}</div>
                </div>
              );
            })}
          </div>
        </div>
      </section>

      {/* CTA Section with Pulse Effect */}
      <section className="relative py-32 bg-gradient-to-br from-[#060140] via-[#0a0560] to-[#030120] overflow-hidden">
        <div className="absolute inset-0 opacity-20">
          <div className="absolute top-10 right-20 w-[40rem] h-[40rem] bg-cyan rounded-full blur-3xl animate-pulse"></div>
          <div className="absolute bottom-10 left-20 w-[35rem] h-[35rem] bg-purple-500 rounded-full blur-3xl animate-pulse" style={{ animationDelay: '1.5s' }}></div>
        </div>

        <div className="relative max-w-5xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <Badge className="mb-6 bg-cyan/20 text-cyan border-cyan/30 px-6 py-2 text-sm animate-bounce-subtle">
            Get Started Today
          </Badge>
          <h2 className="text-4xl sm:text-5xl lg:text-6xl font-bold text-white mb-6 leading-tight">
            Ready to Transform Your
            <span className="block text-cyan mt-2">Genomic Research?</span>
          </h2>
          <p className="text-xl text-blue-100 mb-12 max-w-3xl mx-auto leading-relaxed">
            Join leading research institutions and clinical labs using our validated WES pipeline.
            Start analyzing your exome data today.
          </p>

          <div className="flex flex-col sm:flex-row gap-4 justify-center items-center mb-8">
            <Button
              size="lg"
              className="bg-cyan hover:bg-cyan-light text-white text-lg px-10 py-7 h-auto shadow-2xl hover:shadow-cyan/50 transition-all duration-300 hover:scale-105 hover:-translate-y-1 gooey-button"
              onClick={onSignIn}
            >
              Start Free Trial
              <ArrowRight className="ml-2 h-5 w-5" />
            </Button>
            <Button
              size="lg"
              variant="outline"
              className="bg-transparent hover:bg-white/10 text-white border-2 border-white/40 hover:border-white/60 text-lg px-10 py-7 h-auto backdrop-blur-sm transition-all duration-300 hover:scale-105"
              onClick={() => onNavigate('contact')}
            >
              Contact Sales
            </Button>
          </div>

          <p className="text-sm text-blue-200">
            No credit card required • 2 free jobs per month • Enterprise plans available
          </p>
        </div>
      </section>

      {/* Custom CSS for animations */}
      <style jsx>{`
        @keyframes fadeInDown {
          from {
            opacity: 0;
            transform: translateY(-20px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }

        @keyframes fadeInUp {
          from {
            opacity: 0;
            transform: translateY(20px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }

        @keyframes fadeIn {
          from {
            opacity: 0;
          }
          to {
            opacity: 1;
          }
        }

        @keyframes wordReveal {
          from {
            opacity: 0;
            transform: translateY(30px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }

        @keyframes slideInLeft {
          from {
            opacity: 0;
            transform: translateX(-50px);
          }
          to {
            opacity: 1;
            transform: translateX(0);
          }
        }

        @keyframes slideInRight {
          from {
            opacity: 0;
            transform: translateX(50px);
          }
          to {
            opacity: 1;
            transform: translateX(0);
          }
        }

        @keyframes slideInFeature {
          from {
            opacity: 0;
            transform: translateX(-20px);
          }
          to {
            opacity: 1;
            transform: translateX(0);
          }
        }

        @keyframes bounceSubtle {
          0%, 100% {
            transform: translateY(0);
          }
          50% {
            transform: translateY(-5px);
          }
        }

        @keyframes rotate360 {
          from {
            transform: rotate(0deg);
          }
          to {
            transform: rotate(360deg);
          }
        }

        @keyframes float {
          0%, 100% {
            transform: translateY(0px);
          }
          50% {
            transform: translateY(-10px);
          }
        }

        .animate-fade-in-down {
          animation: fadeInDown 0.8s ease-out;
        }

        .animate-fade-in-up {
          animation: fadeInUp 0.8s ease-out;
        }

        .animate-fade-in {
          animation: fadeIn 1s ease-out;
        }

        .animate-word-reveal {
          animation: wordReveal 0.8s ease-out;
          animation-fill-mode: backwards;
        }

        .slide-in-left {
          animation: slideInLeft 0.8s ease-out;
        }

        .slide-in-right {
          animation: slideInRight 0.8s ease-out;
        }

        .animate-slide-in-feature {
          animation: slideInFeature 0.6s ease-out;
        }

        .animate-bounce-subtle {
          animation: bounceSubtle 2s ease-in-out infinite;
        }

        .animate-float {
          animation: float 3s ease-in-out infinite;
        }

        .value-card {
          animation: fadeInUp 0.8s ease-out;
          animation-fill-mode: backwards;
        }

        .metric-card {
          animation: fadeInUp 0.8s ease-out;
          animation-fill-mode: backwards;
        }

        .group-hover\\:rotate-360:hover {
          animation: rotate360 0.6s ease-out;
        }

        .gooey-button {
          filter: url('#gooey');
        }

        .gooey-button:hover {
          filter: url('#gooey') brightness(1.1);
        }

        .perspective-1000 {
          perspective: 1000px;
        }
      `}</style>

      {/* SVG Filter for Gooey Effect */}
      <svg style={{ position: 'absolute', width: 0, height: 0 }}>
        <defs>
          <filter id="gooey">
            <feGaussianBlur in="SourceGraphic" stdDeviation="5" result="blur" />
            <feColorMatrix
              in="blur"
              mode="matrix"
              values="1 0 0 0 0  0 1 0 0 0  0 0 1 0 0  0 0 0 20 -10"
              result="gooey"
            />
            <feComposite in="SourceGraphic" in2="gooey" operator="atop" />
          </filter>
        </defs>
      </svg>
    </>
  );
}
