'use client';

import { Card, CardContent } from '@/components/ui/card';

export default function ResearchPage() {
  return (
    <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 py-16 sm:py-24">
      <div className="space-y-12">
        <div>
          <h1 className="text-4xl sm:text-5xl font-bold text-slate-900 mb-6">Research Project</h1>
          <div className="bg-yellow-50 border border-yellow-300 rounded-lg p-6 mb-6">
            <p className="text-yellow-900 font-semibold">
              ⚠️ Disclaimer: This is a research project and should NOT be used for clinical purposes. For clinical variant interpretation, please consult with qualified genetic counselors and clinical laboratories.
            </p>
          </div>
        </div>

        <div>
          <h2 className="text-2xl font-bold text-slate-900 mb-4">Project Overview</h2>
          <p className="text-slate-600 leading-relaxed mb-4">
            This platform is an experimental whole exome sequencing analysis system developed to demonstrate scalable approaches to variant interpretation. The project integrates industry-standard bioinformatics tools (BWA, GATK, VEP, ANNOVAR) within a Nextflow pipeline to process FASTQ files through alignment, variant calling, and comprehensive annotation.
          </p>
        </div>

        <div>
          <h2 className="text-2xl font-bold text-slate-900 mb-4">Technical Approach</h2>
          <div className="space-y-4">
            {[
              {
                title: 'Nextflow Pipeline',
                desc: 'Scalable workflow management using Nextflow DSL2 for reproducible bioinformatics analysis with automatic parallelization and resource optimization.'
              },
              {
                title: 'Multi-Step Processing',
                desc: 'FastP quality control → BWA-MEM alignment → GATK base quality recalibration → HaplotypeCaller variant calling → Multi-source annotation.'
              },
              {
                title: 'Comprehensive Annotation',
                desc: 'VEP v113 for variant effect prediction, ANNOVAR for functional annotation, 1000 Genomes population frequencies, and dbNSFP pathogenicity scores.'
              },
              {
                title: 'ACMG Classification',
                desc: 'Automated variant classification following ACMG/AMP 2015 guidelines with evidence-based scoring and clinical interpretation support.'
              },
              {
                title: 'Web Platform',
                desc: 'FastAPI backend with Next.js frontend, Firebase authentication, SQLite job management, and IGV.js genome browser integration.'
              }
            ].map((item, i) => (
              <Card key={i} className="border-slate-200">
                <CardContent className="p-6">
                  <h3 className="font-bold text-slate-900 mb-2">{item.title}</h3>
                  <p className="text-slate-600 text-sm">{item.desc}</p>
                </CardContent>
              </Card>
            ))}
          </div>
        </div>

        <div>
          <h2 className="text-2xl font-bold text-slate-900 mb-4">Reference Genome & Databases</h2>
          <div className="bg-slate-50 border border-slate-200 rounded-lg p-6 space-y-3">
            <p className="text-slate-900">
              <strong>Reference Genome:</strong> GRCh38/hg38 (Genome Reference Consortium Human Build 38)
            </p>
            <p className="text-slate-900">
              <strong>Known Variants:</strong> Homo sapiens assembly38 known indels, Mills & 1000G gold standard indels
            </p>
            <p className="text-slate-900">
              <strong>Annotation Databases:</strong> VEP cache v113, ANNOVAR humandb, 1000 Genomes Project, dbNSFP v4.7, ClinVar
            </p>
            <p className="text-slate-900">
              <strong>Gene Panels:</strong> PanelApp integration, ACMG SF v3.2 secondary findings genes
            </p>
          </div>
        </div>

        <div>
          <h2 className="text-2xl font-bold text-slate-900 mb-4">Pipeline Specifications</h2>
          <div className="grid md:grid-cols-2 gap-6">
            <Card className="border-slate-200">
              <CardContent className="p-6">
                <h3 className="font-bold text-slate-900 mb-3">Input Requirements</h3>
                <ul className="text-slate-600 space-y-2 text-sm">
                  <li>• Paired-end FASTQ files (.fastq.gz)</li>
                  <li>• Whole exome sequencing data</li>
                  <li>• Minimum 50x mean coverage recommended</li>
                  <li>• Illumina platform supported</li>
                </ul>
              </CardContent>
            </Card>
            <Card className="border-slate-200">
              <CardContent className="p-6">
                <h3 className="font-bold text-slate-900 mb-3">Output Files</h3>
                <ul className="text-slate-600 space-y-2 text-sm">
                  <li>• BAM file with BQSR (indexed)</li>
                  <li>• Raw VCF (gzipped & indexed)</li>
                  <li>• Annotated VCF with VEP & ANNOVAR</li>
                  <li>• Filtered TSV with functional annotations</li>
                </ul>
              </CardContent>
            </Card>
          </div>
        </div>

        <div>
          <h2 className="text-2xl font-bold text-slate-900 mb-4">Limitations & Considerations</h2>
          <div className="grid md:grid-cols-2 gap-8">
            <div>
              <h3 className="font-bold text-slate-900 mb-3">Current Limitations</h3>
              <ul className="text-slate-600 space-y-2 text-sm">
                <li>• Research use only - not clinically validated</li>
                <li>• ACMG classifications are automated (require expert review)</li>
                <li>• Limited to WES data (not optimized for WGS or targeted panels)</li>
                <li>• Requires bioinformatics expertise for result interpretation</li>
                <li>• Processing time varies with system resources</li>
              </ul>
            </div>
            <div>
              <h3 className="font-bold text-slate-900 mb-3">Future Enhancements</h3>
              <ul className="text-slate-600 space-y-2 text-sm">
                <li>• Clinical validation studies</li>
                <li>• Enhanced variant filtering options</li>
                <li>• Structural variant detection</li>
                <li>• Pharmacogenomics annotations</li>
                <li>• Export to standard formats (HGVS, VCF 4.3)</li>
                <li>• Batch processing for cohort analysis</li>
              </ul>
            </div>
          </div>
        </div>

        <div className="bg-slate-900 text-white rounded-lg p-8">
          <h3 className="text-lg font-bold mb-4">Pipeline Performance</h3>
          <div className="grid md:grid-cols-3 gap-6">
            <div>
              <div className="text-2xl font-bold text-cyan-400 mb-1">2-4 hours</div>
              <div className="text-slate-300 text-sm">Total processing time</div>
            </div>
            <div>
              <div className="text-2xl font-bold text-cyan-400 mb-1">~80,000</div>
              <div className="text-slate-300 text-sm">Variants per exome</div>
            </div>
            <div>
              <div className="text-2xl font-bold text-cyan-400 mb-1">20+</div>
              <div className="text-slate-300 text-sm">Annotation sources</div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
