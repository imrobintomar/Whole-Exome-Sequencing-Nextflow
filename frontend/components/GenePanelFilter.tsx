'use client';

import { useState } from 'react';
import { panelApi, GenePanel } from '@/lib/api';
import { Search, Filter, ChevronDown, ChevronRight, Star } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Input } from './ui/input';
import { Badge } from './ui/badge';

interface GenePanelFilterProps {
  onPanelSelect?: (panelId: number, genes: string[]) => void;
  onAcmgSelect?: (genes: string[]) => void;
  className?: string;
}

export default function GenePanelFilter({ onPanelSelect, onAcmgSelect, className }: GenePanelFilterProps) {
  const [searchQuery, setSearchQuery] = useState('');
  const [searching, setSearching] = useState(false);
  const [panels, setPanels] = useState<GenePanel[]>([]);
  const [selectedPanelId, setSelectedPanelId] = useState<number | null>(null);
  const [selectedGenes, setSelectedGenes] = useState<string[]>([]);
  const [expandedPanels, setExpandedPanels] = useState<Set<number>>(new Set());
  const [loadingGenes, setLoadingGenes] = useState<number | null>(null);

  const handleSearch = async () => {
    if (!searchQuery.trim()) return;

    setSearching(true);
    try {
      const results = await panelApi.searchPanels(searchQuery);
      setPanels(results);
    } catch (error) {
      console.error('Failed to search panels:', error);
      alert('Failed to search gene panels. Please try again.');
    } finally {
      setSearching(false);
    }
  };

  const togglePanelExpand = async (panelId: number) => {
    const newExpanded = new Set(expandedPanels);

    if (expandedPanels.has(panelId)) {
      newExpanded.delete(panelId);
      setExpandedPanels(newExpanded);
    } else {
      // Load genes if not already loaded
      if (selectedPanelId !== panelId) {
        setLoadingGenes(panelId);
        try {
          const result = await panelApi.getPanelGenes(panelId);
          setSelectedGenes(result.genes);
          setSelectedPanelId(panelId);
        } catch (error) {
          console.error('Failed to load genes:', error);
          alert('Failed to load panel genes. Please try again.');
          setLoadingGenes(null);
          return;
        } finally {
          setLoadingGenes(null);
        }
      }
      newExpanded.add(panelId);
      setExpandedPanels(newExpanded);
    }
  };

  const handleSelectPanel = (panelId: number) => {
    if (selectedPanelId === panelId && selectedGenes.length > 0) {
      onPanelSelect?.(panelId, selectedGenes);
    }
  };

  const handleLoadAcmgSF = async () => {
    setLoadingGenes(-1);
    try {
      const result = await panelApi.getAcmgSecondaryFindings();
      setSelectedGenes(result.genes);
      setSelectedPanelId(-1);
      onAcmgSelect?.(result.genes);
    } catch (error) {
      console.error('Failed to load ACMG SF genes:', error);
      alert('Failed to load ACMG Secondary Findings. Please try again.');
    } finally {
      setLoadingGenes(null);
    }
  };

  return (
    <Card className={className}>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Filter className="h-5 w-5" />
          Gene Panel Filtering
        </CardTitle>
        <CardDescription>
          Filter variants to clinically relevant gene panels
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* ACMG Secondary Findings Quick Access */}
        <div className="border rounded-lg p-4 bg-amber-50 dark:bg-amber-950/20">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <Star className="h-4 w-4 text-amber-600" />
              <div>
                <p className="font-medium text-sm">ACMG Secondary Findings v3.2</p>
                <p className="text-xs text-muted-foreground">81 medically actionable genes</p>
              </div>
            </div>
            <Button
              size="sm"
              variant="outline"
              onClick={handleLoadAcmgSF}
              disabled={loadingGenes === -1}
            >
              {loadingGenes === -1 ? 'Loading...' : 'Apply Filter'}
            </Button>
          </div>
        </div>

        {/* Search Panel */}
        <div className="space-y-2">
          <div className="flex gap-2">
            <Input
              placeholder="Search gene panels (e.g., epilepsy, cardiac, cancer)"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              onKeyPress={(e) => e.key === 'Enter' && handleSearch()}
            />
            <Button onClick={handleSearch} disabled={searching || !searchQuery.trim()}>
              <Search className="h-4 w-4 mr-2" />
              {searching ? 'Searching...' : 'Search'}
            </Button>
          </div>
        </div>

        {/* Results */}
        {panels.length > 0 && (
          <div className="space-y-2 max-h-96 overflow-y-auto">
            <p className="text-sm text-muted-foreground">
              Found {panels.length} panel{panels.length !== 1 ? 's' : ''}
            </p>
            {panels.map((panel) => (
              <div key={panel.id} className="border rounded-lg">
                <div
                  className="p-3 cursor-pointer hover:bg-accent/50 flex items-center justify-between"
                  onClick={() => togglePanelExpand(panel.id)}
                >
                  <div className="flex items-center gap-2 flex-1">
                    {expandedPanels.has(panel.id) ? (
                      <ChevronDown className="h-4 w-4 flex-shrink-0" />
                    ) : (
                      <ChevronRight className="h-4 w-4 flex-shrink-0" />
                    )}
                    <div className="flex-1 min-w-0">
                      <p className="font-medium text-sm truncate">{panel.name}</p>
                      <p className="text-xs text-muted-foreground truncate">
                        {panel.disease_group || panel.disease_sub_group}
                      </p>
                    </div>
                    <Badge variant="secondary" className="flex-shrink-0">
                      v{panel.version}
                    </Badge>
                  </div>
                </div>

                {expandedPanels.has(panel.id) && (
                  <div className="p-3 pt-0 border-t bg-muted/30">
                    {loadingGenes === panel.id ? (
                      <p className="text-sm text-muted-foreground">Loading genes...</p>
                    ) : selectedPanelId === panel.id ? (
                      <>
                        <div className="mb-2">
                          <p className="text-sm font-medium">
                            {selectedGenes.length} genes (Green/High confidence)
                          </p>
                          <div className="flex flex-wrap gap-1 mt-2 max-h-32 overflow-y-auto">
                            {selectedGenes.slice(0, 50).map((gene) => (
                              <Badge key={gene} variant="outline" className="text-xs">
                                {gene}
                              </Badge>
                            ))}
                            {selectedGenes.length > 50 && (
                              <Badge variant="secondary" className="text-xs">
                                +{selectedGenes.length - 50} more
                              </Badge>
                            )}
                          </div>
                        </div>
                        <Button
                          size="sm"
                          className="w-full mt-2"
                          onClick={() => handleSelectPanel(panel.id)}
                        >
                          Apply Filter ({selectedGenes.length} genes)
                        </Button>
                      </>
                    ) : (
                      <p className="text-sm text-muted-foreground">Click to load genes</p>
                    )}
                  </div>
                )}
              </div>
            ))}
          </div>
        )}

        {panels.length === 0 && searchQuery && !searching && (
          <p className="text-sm text-muted-foreground text-center py-4">
            No panels found. Try searching for disease names like "epilepsy", "cardiac", or "cancer"
          </p>
        )}
      </CardContent>
    </Card>
  );
}
