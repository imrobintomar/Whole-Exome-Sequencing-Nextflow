'use client';

import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Database, Zap, Globe, Target, Eye, ArrowRight, Dna, Activity, BarChart3 } from 'lucide-react';

interface AboutPageProps {
  onNavigate?: (page: 'features' | 'publications' | 'contact') => void;
}

export default function AboutPage({ onNavigate }: AboutPageProps) {
  return (
    <div className="min-h-screen">
      {/* Hero Section - Purple Gradient Banner */}
      <section className="relative bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark py-20 overflow-hidden">
        <div className="absolute inset-0 opacity-10">
          <div className="absolute top-10 right-20 w-96 h-96 bg-cyan rounded-full blur-3xl"></div>
        </div>

        <div className="relative max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <Badge variant="cyan" className="mb-6">
            About ATGC Flow
          </Badge>
          <h1 className="text-4xl sm:text-5xl lg:text-6xl font-bold text-white mb-6 leading-tight">
            Advancing Precision Medicine Through Genomic Innovation
          </h1>
          <p className="text-xl text-blue-100 leading-relaxed max-w-3xl mx-auto">
            We are developing a cutting-edge whole exome sequencing analysis platform designed to democratize
            access to clinical-grade variant interpretation for researchers and clinicians worldwide.
          </p>
        </div>
      </section>

      {/* Mission & Vision */}
      <section className="py-24 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid md:grid-cols-2 gap-12">
            <Card className="border-2 border-cyan/20 bg-gradient-to-br from-cyan/5 to-white hover:shadow-xl transition-all duration-300">
              <CardContent className="p-10">
                <div className="w-16 h-16 bg-cyan/10 rounded-2xl flex items-center justify-center mb-6">
                  <Target className="h-8 w-8 text-cyan" />
                </div>
                <h2 className="text-3xl font-bold text-slate-900 mb-4">Our Mission</h2>
                <p className="text-slate-600 leading-relaxed text-lg">
                  To accelerate precision medicine by providing researchers and clinicians with powerful,
                  accessible tools for analyzing exome sequencing data at scale. We believe genomic analysis
                  should be efficient, accurate, and available to all who need it.
                </p>
              </CardContent>
            </Card>

            <Card className="border-2 border-purple-primary/20 bg-gradient-to-br from-purple-primary/5 to-white hover:shadow-xl transition-all duration-300">
              <CardContent className="p-10">
                <div className="w-16 h-16 bg-purple-primary/10 rounded-2xl flex items-center justify-center mb-6">
                  <Eye className="h-8 w-8 text-purple-primary" />
                </div>
                <h2 className="text-3xl font-bold text-slate-900 mb-4">Our Vision</h2>
                <p className="text-slate-600 leading-relaxed text-lg">
                  To become the leading platform for whole exome variant analysis globally, enabling faster
                  diagnosis of genetic diseases and accelerating personalized medicine through accessible,
                  scalable technology and validated methodologies.
                </p>
              </CardContent>
            </Card>
          </div>
        </div>
      </section>

      {/* What Makes Us Different */}
      <section className="py-24 bg-gradient-to-br from-slate-50 to-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-purple-primary text-purple-primary">
              Our Differentiators
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">What Sets Us Apart</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              Clinical-grade analysis meets research innovation
            </p>
          </div>

          <div className="grid md:grid-cols-3 gap-8">
            {[
              {
                icon: Database,
                title: 'Production-Grade Pipeline',
                desc: 'Built on industry-standard tools including BWA-MEM, GATK 4.6, and ANNOVAR. Validated against benchmark datasets with reproducible results.',
                color: 'purple'
              },
              {
                icon: Zap,
                title: 'Comprehensive Annotation',
                desc: 'Multi-source annotation with 1000 Genomes, gnomAD, ClinVar, dbSNP. Automated ACMG/AMP 2015 classification with evidence scoring.',
                color: 'cyan'
              },
              {
                icon: Globe,
                title: 'Interactive Visualization',
                desc: 'IGV.js genome browser integration for real-time variant exploration. BAM/VCF track support with allele fraction analysis.',
                color: 'teal'
              }
            ].map((item, i) => {
              const IconComponent = item.icon;
              const colorClasses = {
                purple: 'bg-purple-primary/10 text-purple-primary border-purple-primary/20',
                cyan: 'bg-cyan/10 text-cyan border-cyan/20',
                teal: 'bg-teal/10 text-teal border-teal/20'
              };
              return (
                <Card key={i} className={`border-2 ${colorClasses[item.color as keyof typeof colorClasses].split(' ').slice(2).join(' ')} hover:shadow-xl transition-all duration-300`}>
                  <CardContent className="p-8">
                    <div className={`w-14 h-14 rounded-xl flex items-center justify-center mb-6 ${colorClasses[item.color as keyof typeof colorClasses].split(' ').slice(0, 2).join(' ')}`}>
                      <IconComponent className="h-7 w-7" />
                    </div>
                    <h3 className="text-xl font-bold text-slate-900 mb-3">{item.title}</h3>
                    <p className="text-slate-600 text-sm leading-relaxed">{item.desc}</p>
                  </CardContent>
                </Card>
              );
            })}
          </div>

          {onNavigate && (
            <div className="text-center mt-12">
              <Button
                onClick={() => onNavigate('features')}
                className="bg-purple-primary hover:bg-purple-light text-white"
              >
                Explore All Features
                <ArrowRight className="ml-2 h-4 w-4" />
              </Button>
            </div>
          )}
        </div>
      </section>

      {/* Research & Development */}
      <section className="py-24 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid lg:grid-cols-2 gap-12 items-center">
            <div>
              <Badge variant="teal" className="mb-4">
                Research & Innovation
              </Badge>
              <h2 className="text-4xl font-bold text-slate-900 mb-6">
                Advancing Genomic Science
              </h2>
              <p className="text-lg text-slate-600 leading-relaxed mb-6">
                We are committed to advancing genomic science through continued research and innovation.
                Our platform integrates best practices from the bioinformatics community and stays current
                with the latest advances in variant interpretation.
              </p>
              <div className="space-y-4">
                {[
                  'Continuous pipeline optimization and validation',
                  'Integration of latest ACMG guidelines and evidence criteria',
                  'Collaboration with clinical laboratories and research institutions',
                  'Open-source contributions to genomics community'
                ].map((item, i) => (
                  <div key={i} className="flex items-start gap-3">
                    <div className="w-6 h-6 bg-cyan/10 rounded-full flex items-center justify-center flex-shrink-0 mt-0.5">
                      <div className="w-2 h-2 bg-cyan rounded-full"></div>
                    </div>
                    <p className="text-slate-700">{item}</p>
                  </div>
                ))}
              </div>
              {onNavigate && (
                <div className="mt-8">
                  <Button
                    onClick={() => onNavigate('publications')}
                    variant="outline"
                    className="border-cyan text-cyan hover:bg-cyan hover:text-white"
                  >
                    View Publications
                    <ArrowRight className="ml-2 h-4 w-4" />
                  </Button>
                </div>
              )}
            </div>

            <div className="grid grid-cols-2 gap-6">
              {[
                { icon: Dna, label: 'Variant Calling', value: '99.5%', desc: 'Accuracy' },
                { icon: Activity, label: 'Pipeline Steps', value: '8', desc: 'Automated stages' },
                { icon: BarChart3, label: 'Annotations', value: '20+', desc: 'Data sources' },
                { icon: Database, label: 'Reference Data', value: '50+ GB', desc: 'Curated databases' }
              ].map((stat, i) => {
                const IconComponent = stat.icon;
                return (
                  <Card key={i} className="border-slate-200 bg-gradient-to-br from-white to-slate-50 hover:shadow-lg transition-all duration-300">
                    <CardContent className="p-6 text-center">
                      <IconComponent className="h-8 w-8 text-purple-primary mx-auto mb-3" />
                      <div className="text-3xl font-bold text-purple-primary mb-1">{stat.value}</div>
                      <div className="font-semibold text-slate-900 text-sm mb-1">{stat.label}</div>
                      <div className="text-xs text-slate-600">{stat.desc}</div>
                    </CardContent>
                  </Card>
                );
              })}
            </div>
          </div>
        </div>
      </section>

      {/* Technology Stack */}
      <section className="py-24 bg-gradient-to-br from-slate-50 via-cyan-50/20 to-purple-50/20">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-cyan text-cyan">
              Technology
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">Our Technology Stack</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              Built on proven, industry-standard bioinformatics tools
            </p>
          </div>

          <div className="grid md:grid-cols-3 gap-8">
            {[
              {
                category: 'Pipeline Tools',
                color: 'purple',
                tools: [
                  { name: 'Nextflow DSL2', desc: 'Workflow orchestration' },
                  { name: 'BWA-MEM', desc: 'Read alignment to GRCh38' },
                  { name: 'GATK 4.6', desc: 'Variant calling & BQSR' },
                  { name: 'FastP', desc: 'Quality control & filtering' }
                ]
              },
              {
                category: 'Annotation Tools',
                color: 'cyan',
                tools: [
                  { name: 'ANNOVAR', desc: 'Functional annotation' },
                  { name: 'SnpEff/SnpSift', desc: 'Variant effect prediction' },
                  { name: 'ClinVar', desc: 'Clinical significance' },
                  { name: 'gnomAD v4', desc: 'Population frequencies' }
                ]
              },
              {
                category: 'Platform & Visualization',
                color: 'teal',
                tools: [
                  { name: 'FastAPI', desc: 'REST API backend' },
                  { name: 'Next.js 14', desc: 'React frontend' },
                  { name: 'IGV.js', desc: 'Genome browser' },
                  { name: 'Firebase', desc: 'Authentication' }
                ]
              }
            ].map((stack, i) => {
              const colorClasses = {
                purple: 'border-purple-primary/30 bg-white',
                cyan: 'border-cyan/30 bg-white',
                teal: 'border-teal/30 bg-white'
              };
              const headerColors = {
                purple: 'text-purple-primary bg-purple-primary/5',
                cyan: 'text-cyan bg-cyan/5',
                teal: 'text-teal bg-teal/5'
              };
              return (
                <Card key={i} className={`border-2 ${colorClasses[stack.color as keyof typeof colorClasses]}`}>
                  <CardContent className="p-8">
                    <div className={`px-4 py-2 rounded-lg ${headerColors[stack.color as keyof typeof headerColors]} mb-6 inline-block`}>
                      <h3 className="font-bold text-sm">{stack.category}</h3>
                    </div>
                    <div className="space-y-4">
                      {stack.tools.map((tool, j) => (
                        <div key={j}>
                          <div className="font-semibold text-slate-900 text-sm">{tool.name}</div>
                          <div className="text-xs text-slate-600">{tool.desc}</div>
                        </div>
                      ))}
                    </div>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="py-20 bg-white">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl sm:text-4xl font-bold text-slate-900 mb-6">
            Ready to Learn More?
          </h2>
          <p className="text-xl text-slate-600 mb-10 max-w-2xl mx-auto">
            Explore our platform capabilities or get in touch with our team to discuss your genomics needs
          </p>
          {onNavigate && (
            <div className="flex flex-col sm:flex-row gap-4 justify-center">
              <Button
                size="lg"
                className="bg-purple-primary hover:bg-purple-light text-white"
                onClick={() => onNavigate('features')}
              >
                Explore Features
                <ArrowRight className="ml-2 h-5 w-5" />
              </Button>
              <Button
                size="lg"
                variant="outline"
                className="border-cyan text-cyan hover:bg-cyan hover:text-white"
                onClick={() => onNavigate('contact')}
              >
                Contact Us
              </Button>
            </div>
          )}
        </div>
      </section>
    </div>
  );
}
