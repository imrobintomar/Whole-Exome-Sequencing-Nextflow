'use client';

import { Card, CardContent } from '@/components/ui/card';
import { Database, Zap, Globe } from 'lucide-react';

export default function AboutPage() {
  return (
    <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 py-16 sm:py-24">
      <div className="space-y-12">
        <div>
          <h1 className="text-4xl sm:text-5xl font-bold text-slate-900 mb-6">About Us</h1>
          <p className="text-lg text-slate-600 leading-relaxed">
            We are developing a cutting-edge whole exome sequencing analysis platform designed to democratize access to production-grade variant interpretation. Our mission is to accelerate precision medicine by providing researchers and clinicians with powerful tools to analyze exome sequencing data at scale.
          </p>
        </div>

        <div className="grid md:grid-cols-2 gap-12 ">
          <div>
            <h2 className="text-2xl font-bold text-slate-900 mb-4">Our Mission</h2>
            <p className="text-slate-600 leading-relaxed">
              We believe genomic analysis should be accessible, efficient, and accurate. Our platform combines advanced bioinformatics pipelines with comprehensive variant annotation to help researchers unlock insights hidden in genetic data.
            </p>
          </div>
          <div>
            <h2 className="text-2xl font-bold text-slate-900 mb-4">Our Vision</h2>
            <p className="text-slate-600 leading-relaxed">
              To become a leading platform for whole exome variant analysis, enabling faster diagnosis of genetic diseases and accelerating personalized medicine through accessible, scalable technology.
            </p>
          </div>
        </div>

        <div>
          <h2 className="text-2xl font-bold text-slate-900 mb-6">What Makes Us Different</h2>
          <div className="grid md:grid-cols-3 gap-6">
            {[
              {
                icon: Database,
                title: 'Production-Grade Pipeline',
                desc: 'Built on Nextflow with industry-standard tools: BWA, GATK, VEP, and ANNOVAR.',
                color: 'blue'
              },
              {
                icon: Zap,
                title: 'Comprehensive Annotation',
                desc: 'Multi-phase processing with VEP annotation, 1000 Genomes frequencies, and ACMG classification.',
                color: 'cyan'
              },
              {
                icon: Globe,
                title: 'Interactive Visualization',
                desc: 'IGV.js browser integration for real-time variant exploration and validation.',
                color: 'purple'
              }
            ].map((item, i) => {
              const IconComponent = item.icon;
              const colorClasses: Record<string, string> = {
                blue: 'bg-blue-100 text-blue-600',
                cyan: 'bg-cyan-100 text-cyan-600',
                purple: 'bg-purple-100 text-purple-600'
              };
              return (
                <Card key={i} className="border-[#02042e] bg-[#02042e]">
                  <CardContent className="p-6">
                    <div className={`w-12 h-12 rounded-lg flex items-center justify-center mb-4 ${colorClasses[item.color]}`}>
                      <IconComponent className="h-6 w-6" />
                    </div>
                    <h3 className="font-bold text-white mb-2">{item.title}</h3>
                    <p className="text-white text-sm">{item.desc}</p>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </div>

        <div className="bg-[#02042e] border border-blue-200 rounded-lg p-8">
          <h3 className="text-lg font-bold text-white mb-2">Research & Development</h3>
          <p className="text-white">
            We are committed to advancing genomic science through continued research and innovation. Our platform integrates best practices from the bioinformatics community and stays current with the latest advances in variant interpretation.
          </p>
        </div>

        <div>
          <h2 className="text-2xl font-bold text-slate-900 mb-6">Our Technology Stack</h2>
          <div className="grid md:grid-cols-2 gap-6">
            <div className="space-y-3">
              <h3 className="font-semibold text-slate-900">Pipeline Tools</h3>
              <ul className="text-slate-600 space-y-2">
                <li>• Nextflow - Workflow orchestration</li>
                <li>• BWA-MEM - Read alignment</li>
                <li>• GATK 4.6 - Variant calling & BQSR</li>
                <li>• VEP - Variant effect prediction</li>
                <li>• ANNOVAR - Functional annotation</li>
              </ul>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
