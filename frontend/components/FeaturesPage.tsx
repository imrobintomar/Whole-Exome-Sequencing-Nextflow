'use client';

import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import {
  FileCheck,
  GitBranch,
  Database,
  Tags,
  Filter,
  Eye,
  FileBarChart,
  Download,
  Layers,
  Activity,
  ArrowRight,
  CheckCircle2,
  Play
} from 'lucide-react';

interface FeatureDetailProps {
  icon: any;
  title: string;
  description: string;
  features: string[];
  technical: string[];
  color: 'purple' | 'cyan' | 'teal';
  imageSide: 'left' | 'right';
}

export default function FeaturesPage() {
  const pipelineSteps = [
    {
      step: '01',
      icon: FileCheck,
      title: 'Quality Control',
      description: 'Pre-processing and quality assessment',
      color: 'purple'
    },
    {
      step: '02',
      icon: GitBranch,
      title: 'Read Alignment',
      description: 'Mapping reads to GRCh38 reference',
      color: 'cyan'
    },
    {
      step: '03',
      icon: Activity,
      title: 'Variant Calling',
      description: 'SNV and indel detection',
      color: 'teal'
    },
    {
      step: '04',
      icon: Database,
      title: 'Annotation',
      description: 'Multi-source functional annotation',
      color: 'purple'
    },
    {
      step: '05',
      icon: Tags,
      title: 'Classification',
      description: 'ACMG/AMP criteria application',
      color: 'cyan'
    },
    {
      step: '06',
      icon: Filter,
      title: 'Filtering',
      description: 'Gene panel and quality filters',
      color: 'teal'
    },
    {
      step: '07',
      icon: Eye,
      title: 'Visualization',
      description: 'Interactive genome browser',
      color: 'purple'
    },
    {
      step: '08',
      icon: FileBarChart,
      title: 'Reporting',
      description: 'Clinical report generation',
      color: 'cyan'
    }
  ];

  const featureDetails: FeatureDetailProps[] = [
    {
      icon: FileCheck,
      title: 'Quality Control & Preprocessing',
      description: 'Comprehensive quality assessment ensures only high-confidence data enters the analysis pipeline',
      features: [
        'Per-base quality score analysis',
        'Adapter sequence trimming',
        'Read length distribution',
        'GC content assessment',
        'Duplicate read identification',
        'Contamination screening'
      ],
      technical: [
        'FASTQ format validation',
        'Phred quality score thresholds',
        'Automatic quality filtering',
        'Detailed QC metrics reports'
      ],
      color: 'purple',
      imageSide: 'right'
    },
    {
      icon: GitBranch,
      title: 'Alignment & Variant Calling',
      description: 'State-of-the-art algorithms for accurate read mapping and variant detection',
      features: [
        'Alignment to GRCh38 reference genome',
        'Base quality score recalibration',
        'Indel realignment',
        'Duplicate marking',
        'SNV and indel calling',
        'Copy number variation detection'
      ],
      technical: [
        'Burrows-Wheeler Aligner (BWA)',
        'GATK HaplotypeCaller',
        'Minimum mapping quality: Q20',
        'Target coverage: 100x mean'
      ],
      color: 'cyan',
      imageSide: 'left'
    },
    {
      icon: Database,
      title: 'Multi-Source Annotation',
      description: 'Comprehensive variant annotation from 20+ curated databases',
      features: [
        'Population frequency (1000G, gnomAD v4)',
        'Clinical significance (ClinVar)',
        'Functional predictions (SIFT, PolyPhen)',
        'Protein domain annotations',
        'Gene constraint metrics',
        'Splice site predictions'
      ],
      technical: [
        'ANNOVAR annotation engine',
        'RefSeq and GENCODE transcripts',
        'dbSNP identifiers',
        'COSMIC cancer mutations'
      ],
      color: 'teal',
      imageSide: 'right'
    },
    {
      icon: Tags,
      title: 'ACMG Classification Engine',
      description: 'Automated variant classification following ACMG/AMP 2015 guidelines',
      features: [
        'Evidence-based pathogenicity scoring',
        'PVS, PS, PM, PP criteria evaluation',
        'BA, BS, BP criteria assessment',
        'Automated classification (P/LP/VUS/LB/B)',
        'Supporting evidence documentation',
        'Conflicting interpretation resolution'
      ],
      technical: [
        'InterVar algorithm integration',
        '28 evidence criteria evaluation',
        'ClinGen expert panel guidelines',
        'Manual override capability'
      ],
      color: 'purple',
      imageSide: 'left'
    },
    {
      icon: Filter,
      title: 'Gene Panel Filtering',
      description: 'Focus analysis on clinically relevant genes for your specific indication',
      features: [
        '500+ curated clinical gene panels',
        'Custom gene list upload',
        'Phenotype-driven filtering (HPO)',
        'Inheritance pattern filtering',
        'Frequency-based filtering',
        'Quality threshold filtering'
      ],
      technical: [
        'PanelApp gene panels',
        'OMIM disease-gene associations',
        'HPO phenotype ontology',
        'Flexible filter combinations'
      ],
      color: 'cyan',
      imageSide: 'right'
    },
    {
      icon: Eye,
      title: 'IGV Genome Browser',
      description: 'Interactive visualization for detailed variant inspection and validation',
      features: [
        'Real-time BAM/VCF viewing',
        'Coverage depth visualization',
        'Allele fraction analysis',
        'Splice junction display',
        'Multi-sample comparison',
        'Bookmark and annotation'
      ],
      technical: [
        'IGV.js integration',
        'WebGL-accelerated rendering',
        'Remote file streaming',
        'Custom track support'
      ],
      color: 'teal',
      imageSide: 'left'
    },
    {
      icon: FileBarChart,
      title: 'Interactive Visualizations',
      description: 'Publication-ready charts and graphs for data exploration',
      features: [
        'Variant distribution plots',
        'Coverage uniformity graphs',
        'Allele frequency histograms',
        'Quality metrics dashboards',
        'Comparison visualizations',
        'Exportable vector graphics'
      ],
      technical: [
        'D3.js and Plotly charts',
        'SVG/PNG export',
        'Responsive design',
        'Interactive tooltips'
      ],
      color: 'purple',
      imageSide: 'right'
    },
    {
      icon: Download,
      title: 'Report Generation',
      description: 'Professional clinical reports with comprehensive variant summaries',
      features: [
        'Structured clinical reports',
        'Variant summary tables',
        'ACMG classification evidence',
        'Gene-disease associations',
        'Literature references',
        'Quality metrics summary'
      ],
      technical: [
        'PDF/Excel export formats',
        'Customizable templates',
        'HGVS nomenclature',
        'CAP/CLIA compliance ready'
      ],
      color: 'cyan',
      imageSide: 'left'
    }
  ];

  return (
    <div className="min-h-screen">
      {/* Hero Section */}
      <section className="relative bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark py-20 overflow-hidden">
        <div className="absolute inset-0 opacity-10">
          <div className="absolute top-10 right-20 w-96 h-96 bg-cyan rounded-full blur-3xl"></div>
        </div>

        <div className="relative max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <Badge variant="cyan" className="mb-6">
            Platform Features
          </Badge>
          <h1 className="text-4xl sm:text-5xl lg:text-6xl font-bold text-white mb-6 leading-tight">
            Comprehensive WES Analysis Platform
          </h1>
          <p className="text-xl text-blue-100 leading-relaxed max-w-3xl mx-auto">
            End-to-end whole exome sequencing analysis with clinical-grade accuracy and automated interpretation
          </p>

          <div className="flex flex-wrap gap-4 justify-center mt-10">
            {[
              'Quality Control',
              'Variant Calling',
              'ACMG Classification',
              'Gene Panels',
              'IGV Browser',
              'Clinical Reports',
              'API Access',
              'Batch Processing'
            ].map((feature, i) => (
              <Badge key={i} variant="cyan" className="text-sm">
                {feature}
              </Badge>
            ))}
          </div>
        </div>
      </section>

      {/* Pipeline Visualization */}
      <section className="py-24 bg-gradient-to-br from-white via-slate-50 to-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-purple-primary text-purple-primary">
              Analysis Pipeline
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">8-Stage Analysis Pipeline</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              From raw FASTQ files to clinical reports in hours, not days
            </p>
          </div>

          <div className="relative">
            {/* Timeline Line */}
            <div className="absolute left-1/2 top-0 bottom-0 w-1 bg-gradient-to-b from-purple-primary via-cyan to-teal transform -translate-x-1/2 hidden lg:block"></div>

            {/* Pipeline Steps */}
            <div className="space-y-8">
              {pipelineSteps.map((step, index) => {
                const IconComponent = step.icon;
                const isEven = index % 2 === 0;
                const colorClasses = {
                  purple: 'bg-purple-primary text-white',
                  cyan: 'bg-cyan text-white',
                  teal: 'bg-teal text-white'
                };

                return (
                  <div key={index} className={`relative grid lg:grid-cols-2 gap-8 items-center ${isEven ? '' : 'lg:flex-row-reverse'}`}>
                    {/* Step Card */}
                    <div className={`${isEven ? 'lg:text-right lg:pr-12' : 'lg:pl-12 lg:col-start-2'}`}>
                      <Card className="border-2 border-slate-200 hover:shadow-xl transition-all duration-300">
                        <CardContent className="p-6">
                          <div className={`flex items-center gap-4 ${isEven ? 'lg:flex-row-reverse' : ''}`}>
                            <div className={`w-16 h-16 rounded-xl flex items-center justify-center ${colorClasses[step.color as keyof typeof colorClasses]}`}>
                              <IconComponent className="h-8 w-8" />
                            </div>
                            <div className={isEven ? 'lg:text-right' : ''}>
                              <div className="text-sm text-slate-500 font-semibold mb-1">STEP {step.step}</div>
                              <h3 className="text-xl font-bold text-slate-900 mb-1">{step.title}</h3>
                              <p className="text-sm text-slate-600">{step.description}</p>
                            </div>
                          </div>
                        </CardContent>
                      </Card>
                    </div>

                    {/* Timeline Dot */}
                    <div className="hidden lg:block absolute left-1/2 transform -translate-x-1/2">
                      <div className={`w-8 h-8 rounded-full ${colorClasses[step.color as keyof typeof colorClasses]} border-4 border-white shadow-lg flex items-center justify-center`}>
                        <CheckCircle2 className="h-4 w-4" />
                      </div>
                    </div>
                  </div>
                );
              })}
            </div>
          </div>

          <div className="text-center mt-16">
            <Card className="inline-block border-2 border-cyan/20 bg-gradient-to-br from-cyan/5 to-white">
              <CardContent className="p-6 flex items-center gap-4">
                <Play className="h-6 w-6 text-cyan" />
                <div className="text-left">
                  <div className="font-semibold text-slate-900">Average Turnaround Time</div>
                  <div className="text-2xl font-bold text-cyan">24-48 hours</div>
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </section>

      {/* Feature Deep Dives */}
      {featureDetails.map((feature, index) => {
        const IconComponent = feature.icon;
        const colorClasses = {
          purple: {
            icon: 'bg-purple-primary/10 text-purple-primary',
            border: 'border-purple-primary/20',
            badge: 'bg-purple-primary text-white'
          },
          cyan: {
            icon: 'bg-cyan/10 text-cyan',
            border: 'border-cyan/20',
            badge: 'bg-cyan text-white'
          },
          teal: {
            icon: 'bg-teal/10 text-teal',
            border: 'border-teal/20',
            badge: 'bg-teal text-white'
          }
        };

        return (
          <section
            key={index}
            className={`py-20 ${index % 2 === 0 ? 'bg-white' : 'bg-gradient-to-br from-slate-50 to-white'}`}
          >
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
              <div className={`grid lg:grid-cols-2 gap-12 items-center ${feature.imageSide === 'left' ? 'lg:flex-row-reverse' : ''}`}>
                {/* Content */}
                <div className={feature.imageSide === 'left' ? 'lg:order-2' : ''}>
                  <div className={`w-16 h-16 rounded-2xl flex items-center justify-center mb-6 ${colorClasses[feature.color].icon}`}>
                    <IconComponent className="h-8 w-8" />
                  </div>

                  <h2 className="text-3xl font-bold text-slate-900 mb-4">{feature.title}</h2>
                  <p className="text-lg text-slate-600 leading-relaxed mb-8">
                    {feature.description}
                  </p>

                  <div className="mb-8">
                    <h3 className="text-lg font-semibold text-slate-900 mb-4">Key Features</h3>
                    <div className="grid sm:grid-cols-2 gap-3">
                      {feature.features.map((item, i) => (
                        <div key={i} className="flex items-start gap-2">
                          <CheckCircle2 className={`h-5 w-5 mt-0.5 flex-shrink-0 ${colorClasses[feature.color].icon.split(' ')[1]}`} />
                          <span className="text-slate-700 text-sm">{item}</span>
                        </div>
                      ))}
                    </div>
                  </div>

                  <div>
                    <h3 className="text-lg font-semibold text-slate-900 mb-4">Technical Details</h3>
                    <div className="flex flex-wrap gap-2">
                      {feature.technical.map((tech, i) => (
                        <span key={i} className={`px-3 py-1 rounded-full text-xs font-medium ${colorClasses[feature.color].badge}`}>
                          {tech}
                        </span>
                      ))}
                    </div>
                  </div>
                </div>

                {/* Visual Placeholder */}
                <div className={feature.imageSide === 'left' ? 'lg:order-1' : ''}>
                  <Card className={`border-2 ${colorClasses[feature.color].border} bg-gradient-to-br from-slate-50 to-white`}>
                    <CardContent className="p-12">
                      <div className="aspect-video bg-gradient-to-br from-slate-100 to-slate-200 rounded-xl flex items-center justify-center">
                        <div className="text-center">
                          <Layers className="h-16 w-16 text-slate-400 mx-auto mb-4" />
                          <div className="text-slate-600 font-medium">Interactive Demo</div>
                          <div className="text-sm text-slate-500 mt-2">Visual representation</div>
                        </div>
                      </div>
                    </CardContent>
                  </Card>
                </div>
              </div>
            </div>
          </section>
        );
      })}

      {/* Technical Specifications */}
      <section className="py-24 bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="cyan" className="mb-6">
              Technical Specifications
            </Badge>
            <h2 className="text-4xl font-bold text-white mb-4">Platform Capabilities</h2>
            <p className="text-xl text-blue-100 max-w-2xl mx-auto">
              Production-ready infrastructure for high-throughput genomic analysis
            </p>
          </div>

          <div className="grid md:grid-cols-3 gap-8">
            {[
              {
                title: 'Input Formats',
                items: ['FASTQ (paired-end)', 'BAM/CRAM files', 'VCF files', 'BED gene panels']
              },
              {
                title: 'Output Formats',
                items: ['Annotated VCF', 'Excel reports', 'PDF clinical reports', 'JSON API data']
              },
              {
                title: 'Reference Genome',
                items: ['GRCh38/hg38', 'RefSeq transcripts', 'GENCODE v43', '50+ GB databases']
              },
              {
                title: 'Compute Requirements',
                items: ['16+ CPU cores', '64 GB RAM', '500 GB storage', 'GPU acceleration']
              },
              {
                title: 'Integration Options',
                items: ['RESTful API', 'CLI interface', 'Web dashboard', 'Batch processing']
              },
              {
                title: 'Quality Standards',
                items: ['CAP/CLIA ready', 'ISO 15189 aligned', 'HIPAA compliant', 'Audit logging']
              }
            ].map((spec, i) => (
              <Card key={i} className="border-2 border-white/20 bg-white/10 backdrop-blur-sm">
                <CardContent className="p-6">
                  <h3 className="text-lg font-bold text-white mb-4">{spec.title}</h3>
                  <ul className="space-y-2">
                    {spec.items.map((item, j) => (
                      <li key={j} className="flex items-start gap-2">
                        <div className="w-1.5 h-1.5 bg-cyan rounded-full mt-2 flex-shrink-0"></div>
                        <span className="text-blue-100 text-sm">{item}</span>
                      </li>
                    ))}
                  </ul>
                </CardContent>
              </Card>
            ))}
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="py-20 bg-white">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl sm:text-4xl font-bold text-slate-900 mb-6">
            Ready to Experience Our Platform?
          </h2>
          <p className="text-xl text-slate-600 mb-10 max-w-2xl mx-auto">
            Start analyzing your exome data with our comprehensive feature set
          </p>
          <div className="flex flex-col sm:flex-row gap-4 justify-center">
            <Button size="lg" className="bg-purple-primary hover:bg-purple-light text-white">
              Start Free Trial
              <ArrowRight className="ml-2 h-5 w-5" />
            </Button>
            <Button
              size="lg"
              variant="outline"
              className="border-cyan text-cyan hover:bg-cyan hover:text-white"
            >
              Schedule Demo
            </Button>
          </div>
        </div>
      </section>
    </div>
  );
}
