'use client';

import { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Badge } from './ui/badge';
import { CheckCircle, XCircle, Info } from 'lucide-react';
import { cn } from '@/lib/utils';
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
} from './ui/dialog';

interface ACMGEvidenceMatrixProps {
  evidence: Record<string, {
    met: boolean;
    strength: string;
    description: string;
    auto_applied: boolean;
  }>;
  evidenceSummary: {
    PVS: number;
    PS: number;
    PM: number;
    PP: number;
    BA: number;
    BS: number;
    BP: number;
  };
  className?: string;
}

// ACMG criteria definitions
const ACMG_CRITERIA = {
  pathogenic: {
    PVS1: 'Null variant in LOF-intolerant gene',
    PS1: 'Same amino acid change as pathogenic',
    PS2: 'De novo (confirmed paternity/maternity)',
    PS3: 'Functional studies support damaging',
    PS4: 'Prevalence in affected > controls',
    PM1: 'In protein functional domain',
    PM2: 'Absent/rare in population databases',
    PM3: 'Detected in trans with pathogenic',
    PM4: 'Protein length change',
    PM5: 'Novel missense at pathogenic site',
    PM6: 'Assumed de novo (no paternity)',
    PP1: 'Cosegregation with disease',
    PP2: 'Missense in constrained gene',
    PP3: 'Multiple computational evidence',
    PP4: 'Phenotype/family history specific',
    PP5: 'Reputable source pathogenic',
  },
  benign: {
    BA1: 'Allele frequency > 5% in population',
    BS1: 'Allele frequency > expected',
    BS2: 'Healthy adult with recessive condition',
    BS3: 'Functional studies show no damaging effect',
    BS4: 'Lack of segregation',
    BP1: 'Missense in LOF-only gene',
    BP2: 'In trans with pathogenic (dominant)',
    BP3: 'In-frame indel in non-repeat',
    BP4: 'Multiple computational evidence benign',
    BP5: 'Variant in case with alternate cause',
    BP6: 'Reputable source benign',
    BP7: 'Synonymous, no splice impact',
  }
};

export default function ACMGEvidenceMatrix({
  evidence,
  evidenceSummary,
  className
}: ACMGEvidenceMatrixProps) {
  const [selectedCriterion, setSelectedCriterion] = useState<string | null>(null);

  const getStrengthColor = (strength: string) => {
    switch (strength.toLowerCase()) {
      case 'very_strong':
      case 'stand_alone':
        return 'bg-purple-100 text-purple-800 border-purple-300';
      case 'strong':
        return 'bg-red-100 text-red-800 border-red-300';
      case 'moderate':
        return 'bg-orange-100 text-orange-800 border-orange-300';
      case 'supporting':
        return 'bg-yellow-100 text-yellow-800 border-yellow-300';
      default:
        return 'bg-gray-100 text-gray-800 border-gray-300';
    }
  };

  const renderCriterionCell = (code: string, description: string, category: 'pathogenic' | 'benign') => {
    const ev = evidence[code];
    const met = ev?.met || false;
    const hasEvidence = ev !== undefined;

    return (
      <div
        key={code}
        className={cn(
          'relative p-3 rounded-lg border-2 transition-all cursor-pointer hover:shadow-md',
          met ? 'border-primary bg-primary/10' : 'border-gray-200 bg-gray-50',
          hasEvidence && !met && 'opacity-60'
        )}
        onClick={() => hasEvidence && setSelectedCriterion(code)}
      >
        <div className="flex items-center justify-between mb-1">
          <span className="font-semibold text-sm">{code}</span>
          {met ? (
            <CheckCircle className="h-4 w-4 text-green-600" />
          ) : hasEvidence ? (
            <XCircle className="h-4 w-4 text-gray-400" />
          ) : (
            <Info className="h-4 w-4 text-gray-300" />
          )}
        </div>
        <p className="text-xs text-muted-foreground line-clamp-2">{description}</p>
        {ev && (
          <Badge className={cn('mt-2 text-[10px]', getStrengthColor(ev.strength))}>
            {ev.strength.replace('_', ' ')}
          </Badge>
        )}
      </div>
    );
  };

  const selectedEvidence = selectedCriterion ? evidence[selectedCriterion] : null;
  const selectedDescription = selectedCriterion
    ? ACMG_CRITERIA.pathogenic[selectedCriterion as keyof typeof ACMG_CRITERIA.pathogenic] ||
      ACMG_CRITERIA.benign[selectedCriterion as keyof typeof ACMG_CRITERIA.benign]
    : '';

  return (
    <>
      <Card className={className}>
        <CardHeader>
          <CardTitle>ACMG Evidence Matrix</CardTitle>
          <CardDescription>
            Click on criteria for details. Green = Met, Gray = Evaluated but not met
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-6">
          {/* Evidence Summary */}
          <div className="grid grid-cols-7 gap-2 p-4 bg-muted rounded-lg">
            <div className="text-center">
              <div className="text-2xl font-bold text-purple-600">{evidenceSummary.PVS}</div>
              <div className="text-xs text-muted-foreground">PVS</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-red-600">{evidenceSummary.PS}</div>
              <div className="text-xs text-muted-foreground">PS</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-orange-600">{evidenceSummary.PM}</div>
              <div className="text-xs text-muted-foreground">PM</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-yellow-600">{evidenceSummary.PP}</div>
              <div className="text-xs text-muted-foreground">PP</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-blue-600">{evidenceSummary.BA}</div>
              <div className="text-xs text-muted-foreground">BA</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-green-600">{evidenceSummary.BS}</div>
              <div className="text-xs text-muted-foreground">BS</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-teal-600">{evidenceSummary.BP}</div>
              <div className="text-xs text-muted-foreground">BP</div>
            </div>
          </div>

          {/* Pathogenic Criteria */}
          <div>
            <h3 className="text-lg font-semibold mb-3 text-red-700">Pathogenic Evidence</h3>

            {/* PVS */}
            <div className="mb-4">
              <h4 className="text-sm font-medium mb-2 text-purple-700">Very Strong (PVS)</h4>
              <div className="grid grid-cols-1 gap-2">
                {renderCriterionCell('PVS1', ACMG_CRITERIA.pathogenic.PVS1, 'pathogenic')}
              </div>
            </div>

            {/* PS */}
            <div className="mb-4">
              <h4 className="text-sm font-medium mb-2 text-red-700">Strong (PS)</h4>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
                {['PS1', 'PS2', 'PS3', 'PS4'].map(code =>
                  renderCriterionCell(
                    code,
                    ACMG_CRITERIA.pathogenic[code as keyof typeof ACMG_CRITERIA.pathogenic],
                    'pathogenic'
                  )
                )}
              </div>
            </div>

            {/* PM */}
            <div className="mb-4">
              <h4 className="text-sm font-medium mb-2 text-orange-700">Moderate (PM)</h4>
              <div className="grid grid-cols-2 md:grid-cols-6 gap-2">
                {['PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6'].map(code =>
                  renderCriterionCell(
                    code,
                    ACMG_CRITERIA.pathogenic[code as keyof typeof ACMG_CRITERIA.pathogenic],
                    'pathogenic'
                  )
                )}
              </div>
            </div>

            {/* PP */}
            <div className="mb-4">
              <h4 className="text-sm font-medium mb-2 text-yellow-700">Supporting (PP)</h4>
              <div className="grid grid-cols-2 md:grid-cols-5 gap-2">
                {['PP1', 'PP2', 'PP3', 'PP4', 'PP5'].map(code =>
                  renderCriterionCell(
                    code,
                    ACMG_CRITERIA.pathogenic[code as keyof typeof ACMG_CRITERIA.pathogenic],
                    'pathogenic'
                  )
                )}
              </div>
            </div>
          </div>

          {/* Benign Criteria */}
          <div>
            <h3 className="text-lg font-semibold mb-3 text-green-700">Benign Evidence</h3>

            {/* BA */}
            <div className="mb-4">
              <h4 className="text-sm font-medium mb-2 text-blue-700">Stand-alone (BA)</h4>
              <div className="grid grid-cols-1 gap-2">
                {renderCriterionCell('BA1', ACMG_CRITERIA.benign.BA1, 'benign')}
              </div>
            </div>

            {/* BS */}
            <div className="mb-4">
              <h4 className="text-sm font-medium mb-2 text-green-700">Strong (BS)</h4>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
                {['BS1', 'BS2', 'BS3', 'BS4'].map(code =>
                  renderCriterionCell(
                    code,
                    ACMG_CRITERIA.benign[code as keyof typeof ACMG_CRITERIA.benign],
                    'benign'
                  )
                )}
              </div>
            </div>

            {/* BP */}
            <div className="mb-4">
              <h4 className="text-sm font-medium mb-2 text-teal-700">Supporting (BP)</h4>
              <div className="grid grid-cols-2 md:grid-cols-7 gap-2">
                {['BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'].map(code =>
                  renderCriterionCell(
                    code,
                    ACMG_CRITERIA.benign[code as keyof typeof ACMG_CRITERIA.benign],
                    'benign'
                  )
                )}
              </div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Evidence Detail Dialog */}
      <Dialog open={selectedCriterion !== null} onOpenChange={() => setSelectedCriterion(null)}>
        <DialogContent>
          <DialogHeader>
            <DialogTitle className="flex items-center gap-2">
              <span>{selectedCriterion}</span>
              {selectedEvidence && (
                <Badge className={getStrengthColor(selectedEvidence.strength)}>
                  {selectedEvidence.strength.replace('_', ' ')}
                </Badge>
              )}
            </DialogTitle>
            <DialogDescription>{selectedDescription}</DialogDescription>
          </DialogHeader>
          {selectedEvidence && (
            <div className="space-y-4">
              <div>
                <h4 className="font-medium mb-2">Status</h4>
                <div className="flex items-center gap-2">
                  {selectedEvidence.met ? (
                    <>
                      <CheckCircle className="h-5 w-5 text-green-600" />
                      <span className="text-green-600 font-medium">Met</span>
                    </>
                  ) : (
                    <>
                      <XCircle className="h-5 w-5 text-gray-400" />
                      <span className="text-gray-600">Not Met</span>
                    </>
                  )}
                </div>
              </div>

              {selectedEvidence.description && (
                <div>
                  <h4 className="font-medium mb-2">Evidence</h4>
                  <p className="text-sm text-muted-foreground bg-muted p-3 rounded">
                    {selectedEvidence.description}
                  </p>
                </div>
              )}

              <div>
                <h4 className="font-medium mb-2">Application</h4>
                <p className="text-sm text-muted-foreground">
                  {selectedEvidence.auto_applied
                    ? 'Automatically applied by classification algorithm'
                    : 'Manually curated'}
                </p>
              </div>
            </div>
          )}
        </DialogContent>
      </Dialog>
    </>
  );
}
