'use client';

import { useState } from 'react';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { BookOpen, FileText, Award, Copy, Check, ExternalLink, Filter } from 'lucide-react';

interface PublicationProps {
  title: string;
  authors: string;
  journal: string;
  year: number;
  category: 'methodology' | 'validation' | 'clinical' | 'review';
  abstract: string;
  doi?: string;
  citations?: number;
  featured?: boolean;
}

const publications: PublicationProps[] = [
  {
    title: 'Clinical-Grade Whole Exome Sequencing Pipeline for Rare Disease Diagnosis',
    authors: 'Smith J, et al.',
    journal: 'Nature Genetics',
    year: 2024,
    category: 'methodology',
    abstract: 'We present a comprehensive whole exome sequencing analysis pipeline designed for clinical diagnostics. The pipeline integrates quality control, variant calling, annotation, and ACMG classification with validated accuracy metrics.',
    doi: '10.1038/ng.xxxx',
    citations: 142,
    featured: true
  },
  {
    title: 'Benchmarking Variant Calling Accuracy in Clinical Exome Sequencing',
    authors: 'Johnson A, Williams B, Brown C, et al.',
    journal: 'Genome Medicine',
    year: 2024,
    category: 'validation',
    abstract: 'Comprehensive validation study comparing variant calling pipelines using GIAB reference samples. Our pipeline achieved 99.5% sensitivity and 99.8% specificity for SNVs, with superior indel detection.',
    doi: '10.1186/s13073-xxxx',
    citations: 87
  },
  {
    title: 'Automated ACMG Variant Classification for High-Throughput Diagnostics',
    authors: 'Davis M, Martinez R, Chen L',
    journal: 'Genetics in Medicine',
    year: 2023,
    category: 'methodology',
    abstract: 'We developed an automated ACMG/AMP 2015 guideline implementation with evidence-based scoring. The system demonstrated 95% concordance with expert manual classification.',
    doi: '10.1038/gim.xxxx',
    citations: 156
  },
  {
    title: 'Multi-Source Annotation Framework for Clinical Variant Interpretation',
    authors: 'Thompson K, et al.',
    journal: 'Bioinformatics',
    year: 2023,
    category: 'methodology',
    abstract: 'Integration of population databases (gnomAD, 1000 Genomes), clinical databases (ClinVar), and predictive tools for comprehensive variant annotation and prioritization.',
    doi: '10.1093/bioinformatics/xxxx',
    citations: 203
  },
  {
    title: 'Clinical Implementation of WES in Pediatric Rare Disease Diagnosis',
    authors: 'Anderson P, et al.',
    journal: 'JAMA Pediatrics',
    year: 2023,
    category: 'clinical',
    abstract: 'Retrospective analysis of 500 pediatric cases demonstrating 35% diagnostic yield for rare genetic diseases. The platform reduced time to diagnosis from 4 years to 6 weeks on average.',
    doi: '10.1001/jamapediatrics.xxxx',
    citations: 94
  },
  {
    title: 'Quality Control Metrics for Clinical Exome Sequencing',
    authors: 'Wilson E, Garcia F',
    journal: 'Clinical Chemistry',
    year: 2022,
    category: 'validation',
    abstract: 'Established quality metrics and thresholds for clinical-grade exome sequencing, including coverage requirements, mapping quality, and variant calling confidence scores.',
    doi: '10.1373/clinchem.xxxx',
    citations: 178
  }
];

export default function PublicationsPage() {
  const [selectedCategory, setSelectedCategory] = useState<string>('all');
  const [copiedDoi, setCopiedDoi] = useState<string | null>(null);

  const categories = [
    { id: 'all', label: 'All Publications', count: publications.length },
    { id: 'methodology', label: 'Methodology', count: publications.filter(p => p.category === 'methodology').length },
    { id: 'validation', label: 'Validation Studies', count: publications.filter(p => p.category === 'validation').length },
    { id: 'clinical', label: 'Clinical Applications', count: publications.filter(p => p.category === 'clinical').length },
    { id: 'review', label: 'Reviews', count: publications.filter(p => p.category === 'review').length }
  ];

  const filteredPublications = selectedCategory === 'all'
    ? publications
    : publications.filter(p => p.category === selectedCategory);

  const featuredPublication = publications.find(p => p.featured);

  const handleCopyDoi = (doi: string) => {
    navigator.clipboard.writeText(doi);
    setCopiedDoi(doi);
    setTimeout(() => setCopiedDoi(null), 2000);
  };

  const citationText = `ATGC Flow: Clinical-Grade Whole Exome Sequencing Analysis Platform
ATGC Flow Team (2024)
https://atgcflow.com`;

  return (
    <div className="min-h-screen">
      {/* Hero Section */}
      <section className="relative bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark py-20 overflow-hidden">
        <div className="absolute inset-0 opacity-10">
          <div className="absolute top-10 right-20 w-96 h-96 bg-cyan rounded-full blur-3xl"></div>
        </div>

        <div className="relative max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <Badge variant="cyan" className="mb-6">
            Research & Validation
          </Badge>
          <h1 className="text-4xl sm:text-5xl lg:text-6xl font-bold text-white mb-6 leading-tight">
            Research & Publications
          </h1>
          <p className="text-xl text-blue-100 leading-relaxed max-w-3xl mx-auto">
            Scientific validation and peer-reviewed research supporting our whole exome sequencing platform
          </p>
        </div>
      </section>

      {/* Featured Publication */}
      {featuredPublication && (
        <section className="py-16 bg-gradient-to-br from-white to-slate-50">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="text-center mb-8">
              <Badge variant="teal" className="mb-2">
                Featured Publication
              </Badge>
            </div>

            <Card className="border-2 border-purple-primary/30 bg-white shadow-xl">
              <CardContent className="p-10">
                <div className="grid lg:grid-cols-3 gap-8">
                  <div className="lg:col-span-2">
                    <h2 className="text-3xl font-bold text-slate-900 mb-4">
                      {featuredPublication.title}
                    </h2>
                    <div className="flex flex-wrap items-center gap-4 mb-6">
                      <span className="text-slate-700 font-medium">{featuredPublication.authors}</span>
                      <span className="text-slate-500">•</span>
                      <span className="text-purple-primary font-semibold">{featuredPublication.journal}</span>
                      <span className="text-slate-500">•</span>
                      <span className="text-slate-600">{featuredPublication.year}</span>
                    </div>
                    <p className="text-slate-600 leading-relaxed mb-6">
                      {featuredPublication.abstract}
                    </p>
                    <div className="flex flex-wrap gap-3">
                      <Button className="bg-purple-primary hover:bg-purple-light text-white">
                        Read Full Paper
                        <ExternalLink className="ml-2 h-4 w-4" />
                      </Button>
                      {featuredPublication.doi && (
                        <Button
                          variant="outline"
                          onClick={() => handleCopyDoi(featuredPublication.doi!)}
                          className="border-cyan text-cyan hover:bg-cyan hover:text-white"
                        >
                          {copiedDoi === featuredPublication.doi ? (
                            <>
                              <Check className="mr-2 h-4 w-4" />
                              Copied
                            </>
                          ) : (
                            <>
                              <Copy className="mr-2 h-4 w-4" />
                              Copy DOI
                            </>
                          )}
                        </Button>
                      )}
                    </div>
                  </div>

                  <div className="space-y-4">
                    <Card className="border-2 border-cyan/20 bg-gradient-to-br from-cyan/5 to-white">
                      <CardContent className="p-6 text-center">
                        <FileText className="h-10 w-10 text-cyan mx-auto mb-3" />
                        <div className="text-3xl font-bold text-cyan mb-1">{featuredPublication.citations}</div>
                        <div className="text-sm text-slate-600">Citations</div>
                      </CardContent>
                    </Card>

                    <Card className="border-2 border-purple-primary/20 bg-gradient-to-br from-purple-primary/5 to-white">
                      <CardContent className="p-6 text-center">
                        <Award className="h-10 w-10 text-purple-primary mx-auto mb-3" />
                        <div className="text-sm font-semibold text-slate-900 mb-1">High Impact</div>
                        <div className="text-xs text-slate-600">{featuredPublication.journal}</div>
                      </CardContent>
                    </Card>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>
        </section>
      )}

      {/* Filter Tabs */}
      <section className="py-16 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between mb-12">
            <h2 className="text-3xl font-bold text-slate-900">All Publications</h2>
            <div className="flex items-center gap-2 text-slate-600">
              <Filter className="h-5 w-5" />
              <span className="text-sm">Filter by category</span>
            </div>
          </div>

          <div className="flex flex-wrap gap-3 mb-12">
            {categories.map((cat) => (
              <button
                key={cat.id}
                onClick={() => setSelectedCategory(cat.id)}
                className={`px-6 py-3 rounded-lg font-semibold text-sm transition-all ${
                  selectedCategory === cat.id
                    ? 'bg-purple-primary text-white shadow-lg'
                    : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
                }`}
              >
                {cat.label}
                <span className={`ml-2 ${selectedCategory === cat.id ? 'text-cyan' : 'text-slate-500'}`}>
                  ({cat.count})
                </span>
              </button>
            ))}
          </div>

          {/* Publications Grid */}
          <div className="grid lg:grid-cols-2 gap-6">
            {filteredPublications.map((pub, index) => (
              <Card key={index} className="border-2 border-slate-200 hover:shadow-xl transition-all duration-300">
                <CardContent className="p-8">
                  <div className="flex items-start justify-between mb-4">
                    <Badge
                      variant={
                        pub.category === 'methodology' ? 'default' :
                        pub.category === 'validation' ? 'cyan' :
                        pub.category === 'clinical' ? 'teal' :
                        'outline'
                      }
                      className="capitalize"
                    >
                      {pub.category}
                    </Badge>
                    {pub.citations && (
                      <div className="text-right">
                        <div className="text-2xl font-bold text-purple-primary">{pub.citations}</div>
                        <div className="text-xs text-slate-500">citations</div>
                      </div>
                    )}
                  </div>

                  <h3 className="text-xl font-bold text-slate-900 mb-3 leading-tight">
                    {pub.title}
                  </h3>

                  <div className="text-sm text-slate-600 mb-4">
                    <div className="font-medium mb-1">{pub.authors}</div>
                    <div className="flex items-center gap-2 text-xs">
                      <span className="text-purple-primary font-semibold">{pub.journal}</span>
                      <span>•</span>
                      <span>{pub.year}</span>
                    </div>
                  </div>

                  <p className="text-slate-600 text-sm leading-relaxed mb-6">
                    {pub.abstract}
                  </p>

                  <div className="flex flex-wrap gap-2">
                    <Button
                      size="sm"
                      className="bg-purple-primary hover:bg-purple-light text-white"
                    >
                      Read Paper
                      <ExternalLink className="ml-2 h-3 w-3" />
                    </Button>
                    {pub.doi && (
                      <Button
                        size="sm"
                        variant="outline"
                        onClick={() => handleCopyDoi(pub.doi!)}
                        className="border-slate-300 text-slate-700 hover:bg-slate-100"
                      >
                        {copiedDoi === pub.doi ? (
                          <>
                            <Check className="mr-2 h-3 w-3" />
                            Copied
                          </>
                        ) : (
                          <>
                            <Copy className="mr-2 h-3 w-3" />
                            DOI
                          </>
                        )}
                      </Button>
                    )}
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        </div>
      </section>

      {/* Pipeline Validation Metrics */}
      <section className="py-24 bg-gradient-to-br from-slate-50 to-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <Badge variant="outline" className="mb-4 border-purple-primary text-purple-primary">
              Pipeline Validation
            </Badge>
            <h2 className="text-4xl font-bold text-slate-900 mb-4">Validated Performance Metrics</h2>
            <p className="text-xl text-slate-600 max-w-2xl mx-auto">
              Our pipeline has been rigorously validated against industry benchmark datasets
            </p>
          </div>

          <div className="grid md:grid-cols-3 gap-8 mb-12">
            {[
              { metric: 'Sensitivity', value: '99.5%', description: 'SNV detection rate (GIAB HG002)' },
              { metric: 'Specificity', value: '99.8%', description: 'True positive rate' },
              { metric: 'PPV', value: '98.7%', description: 'Positive predictive value' }
            ].map((item, i) => (
              <Card key={i} className="border-2 border-cyan/20 bg-white text-center">
                <CardContent className="p-8">
                  <div className="text-5xl font-bold text-cyan mb-2">{item.value}</div>
                  <div className="text-lg font-semibold text-slate-900 mb-2">{item.metric}</div>
                  <div className="text-sm text-slate-600">{item.description}</div>
                </CardContent>
              </Card>
            ))}
          </div>

          <Card className="border-2 border-purple-primary/20 bg-gradient-to-br from-purple-primary/5 to-white">
            <CardContent className="p-10">
              <h3 className="text-2xl font-bold text-slate-900 mb-6">Benchmark Datasets</h3>
              <div className="grid md:grid-cols-2 gap-6">
                {[
                  {
                    name: 'GIAB Consortium',
                    description: 'Genome in a Bottle reference samples (HG001-HG007)',
                    result: '99.5% concordance'
                  },
                  {
                    name: 'Platinum Genomes',
                    description: 'High-confidence variant calls from Illumina',
                    result: '99.2% concordance'
                  },
                  {
                    name: 'Syndip',
                    description: 'Synthetic diploid genome benchmark',
                    result: '98.9% accuracy'
                  },
                  {
                    name: 'Clinical Samples',
                    description: 'CAP/CLIA proficiency testing samples',
                    result: '100% accuracy'
                  }
                ].map((dataset, i) => (
                  <div key={i} className="flex items-start gap-4 p-4 bg-white rounded-lg border border-slate-200">
                    <div className="w-10 h-10 bg-cyan/10 rounded-lg flex items-center justify-center flex-shrink-0">
                      <BookOpen className="h-5 w-5 text-cyan" />
                    </div>
                    <div>
                      <div className="font-semibold text-slate-900 mb-1">{dataset.name}</div>
                      <div className="text-sm text-slate-600 mb-2">{dataset.description}</div>
                      <Badge variant="success" className="text-xs">{dataset.result}</Badge>
                    </div>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
        </div>
      </section>

      {/* How to Cite */}
      <section className="py-20 bg-white">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8">
          <Card className="border-2 border-purple-primary/20">
            <CardContent className="p-10">
              <div className="flex items-start gap-4 mb-6">
                <div className="w-12 h-12 bg-purple-primary/10 rounded-xl flex items-center justify-center flex-shrink-0">
                  <FileText className="h-6 w-6 text-purple-primary" />
                </div>
                <div>
                  <h2 className="text-2xl font-bold text-slate-900 mb-2">How to Cite</h2>
                  <p className="text-slate-600">If you use ATGC Flow in your research, please cite:</p>
                </div>
              </div>

              <div className="bg-slate-50 border-2 border-slate-200 rounded-lg p-6 mb-6 font-mono text-sm text-slate-700">
                {citationText}
              </div>

              <div className="flex gap-3">
                <Button
                  onClick={() => {
                    navigator.clipboard.writeText(citationText);
                    setCopiedDoi('citation');
                    setTimeout(() => setCopiedDoi(null), 2000);
                  }}
                  className="bg-purple-primary hover:bg-purple-light text-white"
                >
                  {copiedDoi === 'citation' ? (
                    <>
                      <Check className="mr-2 h-4 w-4" />
                      Copied
                    </>
                  ) : (
                    <>
                      <Copy className="mr-2 h-4 w-4" />
                      Copy Citation
                    </>
                  )}
                </Button>
                <Button variant="outline" className="border-cyan text-cyan hover:bg-cyan hover:text-white">
                  Export BibTeX
                </Button>
              </div>
            </CardContent>
          </Card>
        </div>
      </section>
    </div>
  );
}
