'use client';

import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import {
  Microscope,
  Building2,
  Pill,
  User,
  Clock,
  TrendingUp,
  CheckCircle2,
  ArrowRight,
  Download,
  Target,
  LucideIcon
} from 'lucide-react';

interface UseCaseProps {
  icon: LucideIcon;
  title: string;
  subtitle: string;
  description: string;
  challenges: string[];
  solutions: string[];
  outcomes: { metric: string; value: string }[];
  color: 'purple' | 'cyan' | 'teal';
}

interface CaseStudyProps {
  title: string;
  category: string;
  problem: string;
  solution: string;
  results: string[];
  color: 'purple' | 'cyan' | 'teal';
}

export default function UseCasesPage() {
  const useCases: UseCaseProps[] = [
    {
      icon: Microscope,
      title: 'Research Laboratories',
      subtitle: 'Discovery Research & Population Studies',
      description: 'Accelerate genomic discovery research with scalable analysis infrastructure designed for high-throughput variant detection and interpretation.',
      challenges: [
        'Managing large-scale sequencing datasets',
        'Reproducible bioinformatics workflows',
        'Variant interpretation complexity',
        'Limited computational resources',
        'Time-consuming manual analysis'
      ],
      solutions: [
        'Automated end-to-end pipeline',
        'Scalable cloud infrastructure',
        'Standardized annotation databases',
        'Batch processing capabilities',
        'Customizable gene panels',
        'API integration for automation'
      ],
      outcomes: [
        { metric: 'Analysis Time', value: '75% reduction' },
        { metric: 'Data Throughput', value: '10x increase' },
        { metric: 'Reproducibility', value: '100%' }
      ],
      color: 'purple'
    },
    {
      icon: Building2,
      title: 'Clinical Diagnostics',
      subtitle: 'Rare Disease & Genetic Disorder Diagnosis',
      description: 'Enable faster, more accurate genetic diagnoses with clinical-grade variant calling and ACMG classification for rare disease patients.',
      challenges: [
        'Diagnostic odyssey duration',
        'Variant interpretation accuracy',
        'Turnaround time pressure',
        'CAP/CLIA compliance requirements',
        'Manual curation burden'
      ],
      solutions: [
        'ACMG/AMP automated classification',
        'Curated disease gene panels',
        'ClinVar integration',
        'IGV visualization for validation',
        'Clinical report generation',
        'Quality assurance workflows'
      ],
      outcomes: [
        { metric: 'Diagnostic Yield', value: '35-40%' },
        { metric: 'Time to Diagnosis', value: '6 weeks' },
        { metric: 'Accuracy', value: '99.5%+' }
      ],
      color: 'cyan'
    },
    {
      icon: Pill,
      title: 'Pharmaceutical R&D',
      subtitle: 'Drug Target Discovery & Clinical Trials',
      description: 'Identify novel drug targets and stratify patients for clinical trials using comprehensive genomic analysis and variant databases.',
      challenges: [
        'Target identification efficiency',
        'Patient stratification complexity',
        'Variant-drug associations',
        'Pharmacogenomics interpretation',
        'Multi-site data harmonization'
      ],
      solutions: [
        'Comprehensive variant annotation',
        'Population frequency filtering',
        'Gene-disease associations',
        'Pharmacogenomics databases',
        'Batch cohort analysis',
        'API for LIMS integration'
      ],
      outcomes: [
        { metric: 'Target Discovery', value: '3x faster' },
        { metric: 'Trial Enrollment', value: '50% improvement' },
        { metric: 'Data Integration', value: 'Seamless' }
      ],
      color: 'teal'
    },
    {
      icon: User,
      title: 'Individual Researchers',
      subtitle: 'Academic Projects & Pilot Studies',
      description: 'Access enterprise-grade genomics analysis tools at affordable prices, perfect for thesis projects, pilot studies, and academic publications.',
      challenges: [
        'Limited computational resources',
        'Budget constraints',
        'Steep learning curve',
        'Lack of bioinformatics expertise',
        'Publication requirements'
      ],
      solutions: [
        'User-friendly web interface',
        'Flexible pricing (pay-per-use)',
        'Comprehensive documentation',
        'Pre-configured gene panels',
        'Publication-ready reports',
        'Educational resources'
      ],
      outcomes: [
        { metric: 'Cost Savings', value: '80% vs in-house' },
        { metric: 'Time to Results', value: '< 48 hours' },
        { metric: 'Publications', value: 'Citation-ready' }
      ],
      color: 'purple'
    }
  ];

  const caseStudies: CaseStudyProps[] = [
    {
      title: 'Rare Disease Diagnosis in Pediatric Neurology',
      category: 'Clinical Diagnostics',
      problem: 'A pediatric neurology department faced a 4-year average diagnostic odyssey for patients with suspected genetic disorders. Manual variant interpretation was time-consuming and required extensive molecular geneticist time.',
      solution: 'Implemented ATGC Flow with curated neurology gene panels and automated ACMG classification. Integrated IGV browser for rapid variant validation and generated structured clinical reports.',
      results: [
        'Diagnostic yield increased from 28% to 37%',
        'Average time to diagnosis reduced to 6 weeks',
        'Molecular geneticist time reduced by 60%',
        'Identified 42 novel variants in 12 months'
      ],
      color: 'cyan'
    },
    {
      title: 'Population Genomics Study in Southeast Asia',
      category: 'Research Laboratory',
      problem: 'A research consortium needed to analyze 5,000 exomes to identify population-specific variants and disease associations. Existing infrastructure could not scale, and inconsistent analysis methods across sites caused issues.',
      solution: 'Deployed ATGC Flow API with standardized pipeline across 8 research sites. Implemented batch processing with automated quality control and centralized variant database.',
      results: [
        'Analyzed 5,000 exomes in 3 months',
        '100% pipeline consistency across sites',
        'Identified 2,450 population-specific variants',
        'Published 3 high-impact papers'
      ],
      color: 'purple'
    },
    {
      title: 'Pharmacogenomics in Oncology Clinical Trial',
      category: 'Pharmaceutical R&D',
      problem: 'A phase II oncology trial required rapid identification of patients with specific biomarkers for targeted therapy enrollment. Manual screening was too slow, missing enrollment windows.',
      solution: 'Integrated ATGC Flow API with trial enrollment system. Created custom gene panels for biomarker detection with automated flagging of eligible patients.',
      results: [
        'Enrollment rate increased by 45%',
        'Screening turnaround time: 24 hours',
        'Successfully stratified 180 patients',
        'Trial completed 3 months ahead of schedule'
      ],
      color: 'teal'
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
            Use Cases
          </Badge>
          <h1 className="text-4xl sm:text-5xl lg:text-6xl font-bold text-white mb-6 leading-tight">
            Empowering Genomic Discovery Across Disciplines
          </h1>
          <p className="text-xl text-blue-100 leading-relaxed max-w-3xl mx-auto">
            From research labs to clinical diagnostics, see how ATGC Flow accelerates genomic analysis for diverse applications
          </p>
        </div>
      </section>

      {/* Use Cases - Detailed Sections */}
      {useCases.map((useCase, index) => {
        const IconComponent = useCase.icon;
        const colorClasses = {
          purple: {
            bg: 'bg-purple-primary',
            bgLight: 'bg-purple-primary/10',
            text: 'text-purple-primary',
            border: 'border-purple-primary/20',
            gradient: 'from-purple-primary/5 to-white'
          },
          cyan: {
            bg: 'bg-cyan',
            bgLight: 'bg-cyan/10',
            text: 'text-cyan',
            border: 'border-cyan/20',
            gradient: 'from-cyan/5 to-white'
          },
          teal: {
            bg: 'bg-teal',
            bgLight: 'bg-teal/10',
            text: 'text-teal',
            border: 'border-teal/20',
            gradient: 'from-teal/5 to-white'
          }
        };

        return (
          <section key={index} className={`py-24 ${index % 2 === 0 ? 'bg-white' : 'bg-gradient-to-br from-slate-50 to-white'}`}>
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
              {/* Header */}
              <div className="max-w-3xl mb-12">
                <div className={`w-16 h-16 ${colorClasses[useCase.color].bgLight} rounded-2xl flex items-center justify-center mb-6`}>
                  <IconComponent className={`h-8 w-8 ${colorClasses[useCase.color].text}`} />
                </div>
                <h2 className="text-4xl font-bold text-slate-900 mb-3">{useCase.title}</h2>
                <p className={`text-xl font-semibold ${colorClasses[useCase.color].text} mb-4`}>
                  {useCase.subtitle}
                </p>
                <p className="text-lg text-slate-600 leading-relaxed">
                  {useCase.description}
                </p>
              </div>

              {/* Content Grid */}
              <div className="grid lg:grid-cols-3 gap-8">
                {/* Challenges */}
                <Card className="border-2 border-slate-200">
                  <CardContent className="p-8">
                    <div className="flex items-center gap-3 mb-6">
                      <Target className="h-6 w-6 text-slate-400" />
                      <h3 className="text-xl font-bold text-slate-900">Challenges</h3>
                    </div>
                    <ul className="space-y-3">
                      {useCase.challenges.map((challenge, i) => (
                        <li key={i} className="flex items-start gap-3 text-slate-600 text-sm">
                          <div className="w-1.5 h-1.5 bg-slate-400 rounded-full mt-2 flex-shrink-0"></div>
                          <span>{challenge}</span>
                        </li>
                      ))}
                    </ul>
                  </CardContent>
                </Card>

                {/* Solutions */}
                <Card className={`border-2 ${colorClasses[useCase.color].border} bg-gradient-to-br ${colorClasses[useCase.color].gradient}`}>
                  <CardContent className="p-8">
                    <div className="flex items-center gap-3 mb-6">
                      <CheckCircle2 className={`h-6 w-6 ${colorClasses[useCase.color].text}`} />
                      <h3 className="text-xl font-bold text-slate-900">Our Solutions</h3>
                    </div>
                    <ul className="space-y-3">
                      {useCase.solutions.map((solution, i) => (
                        <li key={i} className="flex items-start gap-3 text-slate-700 text-sm">
                          <CheckCircle2 className={`h-4 w-4 ${colorClasses[useCase.color].text} flex-shrink-0 mt-0.5`} />
                          <span className="font-medium">{solution}</span>
                        </li>
                      ))}
                    </ul>
                  </CardContent>
                </Card>

                {/* Outcomes */}
                <Card className={`border-2 ${colorClasses[useCase.color].border} ${colorClasses[useCase.color].bg} text-white`}>
                  <CardContent className="p-8">
                    <div className="flex items-center gap-3 mb-6">
                      <TrendingUp className="h-6 w-6" />
                      <h3 className="text-xl font-bold">Success Metrics</h3>
                    </div>
                    <div className="space-y-6">
                      {useCase.outcomes.map((outcome, i) => (
                        <div key={i}>
                          <div className="text-sm font-medium mb-1 opacity-90">{outcome.metric}</div>
                          <div className="text-3xl font-bold">{outcome.value}</div>
                        </div>
                      ))}
                    </div>
                  </CardContent>
                </Card>
              </div>

              {/* CTA */}
              <div className="mt-10 text-center">
                <Button className={`${colorClasses[useCase.color].bg} hover:opacity-90 text-white`}>
                  Learn More About This Use Case
                  <ArrowRight className="ml-2 h-4 w-4" />
                </Button>
              </div>
            </div>
          </section>
        );
      })}

      {/* Case Studies Section */}
      <section className="py-24 bg-gradient-to-br from-slate-900 via-purple-primary to-purple-light">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="cyan" className="mb-6">
              Real-World Impact
            </Badge>
            <h2 className="text-4xl font-bold text-white mb-4">Case Studies</h2>
            <p className="text-xl text-blue-100 max-w-2xl mx-auto">
              See how organizations are using ATGC Flow to transform their genomic workflows
            </p>
          </div>

          <div className="grid lg:grid-cols-3 gap-8">
            {caseStudies.map((study, index) => {
              const colorClasses = {
                purple: 'border-purple-light/50',
                cyan: 'border-cyan/50',
                teal: 'border-teal/50'
              };

              return (
                <Card key={index} className={`border-2 ${colorClasses[study.color]} bg-white/10 backdrop-blur-sm hover:bg-white/15 transition-all duration-300`}>
                  <CardContent className="p-8">
                    <Badge variant={study.color} className="mb-4">
                      {study.category}
                    </Badge>

                    <h3 className="text-xl font-bold text-white mb-6 leading-tight">
                      {study.title}
                    </h3>

                    <div className="space-y-6">
                      <div>
                        <div className="text-sm font-semibold text-cyan mb-2">Problem</div>
                        <p className="text-sm text-blue-100 leading-relaxed">{study.problem}</p>
                      </div>

                      <div>
                        <div className="text-sm font-semibold text-cyan mb-2">Solution</div>
                        <p className="text-sm text-blue-100 leading-relaxed">{study.solution}</p>
                      </div>

                      <div>
                        <div className="text-sm font-semibold text-cyan mb-3">Results</div>
                        <ul className="space-y-2">
                          {study.results.map((result, i) => (
                            <li key={i} className="flex items-start gap-2 text-sm text-blue-100">
                              <CheckCircle2 className="h-4 w-4 text-cyan flex-shrink-0 mt-0.5" />
                              <span>{result}</span>
                            </li>
                          ))}
                        </ul>
                      </div>
                    </div>

                    <Button
                      variant="outline"
                      className="w-full mt-6 border-white/30 text-white hover:bg-white/20"
                    >
                      <Download className="mr-2 h-4 w-4" />
                      Download Case Study
                    </Button>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </div>
      </section>

      {/* Common Benefits Section */}
      <section className="py-24 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-purple-primary text-purple-primary">
              Universal Benefits
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">Why Organizations Choose ATGC Flow</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              Key advantages that benefit all use cases
            </p>
          </div>

          <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-8">
            {[
              {
                icon: Clock,
                title: 'Faster Results',
                description: 'From days to hours with automated pipeline',
                color: 'purple'
              },
              {
                icon: CheckCircle2,
                title: 'Higher Accuracy',
                description: '99.5%+ sensitivity with validated methods',
                color: 'cyan'
              },
              {
                icon: TrendingUp,
                title: 'Scalable',
                description: 'From single samples to thousands of exomes',
                color: 'teal'
              },
              {
                icon: Target,
                title: 'Cost Effective',
                description: 'Reduce infrastructure and personnel costs',
                color: 'purple'
              }
            ].map((benefit, i) => {
              const IconComponent = benefit.icon;
              const colorClasses = {
                purple: 'bg-purple-primary/10 text-purple-primary',
                cyan: 'bg-cyan/10 text-cyan',
                teal: 'bg-teal/10 text-teal'
              };

              return (
                <Card key={i} className="border-2 border-slate-200 text-center hover:shadow-xl transition-all duration-300">
                  <CardContent className="p-8">
                    <div className={`w-16 h-16 ${colorClasses[benefit.color as keyof typeof colorClasses]} rounded-2xl flex items-center justify-center mx-auto mb-4`}>
                      <IconComponent className="h-8 w-8" />
                    </div>
                    <h3 className="text-lg font-bold text-slate-900 mb-2">{benefit.title}</h3>
                    <p className="text-sm text-slate-600">{benefit.description}</p>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </div>
      </section>

      {/* Final CTA */}
      <section className="py-20 bg-gradient-to-br from-slate-50 to-white">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl sm:text-4xl font-bold text-slate-900 mb-6">
            Ready to Transform Your Genomic Workflow?
          </h2>
          <p className="text-xl text-slate-600 mb-10 max-w-2xl mx-auto">
            Join organizations worldwide using ATGC Flow for their exome sequencing analysis
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
              Schedule Consultation
            </Button>
          </div>
        </div>
      </section>
    </div>
  );
}
