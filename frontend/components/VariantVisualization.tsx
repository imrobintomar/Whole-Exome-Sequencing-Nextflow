"use client"

import { useState, useEffect } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "./ui/card"
import { Button } from "./ui/button"
import { ArrowLeft, Loader2 } from "lucide-react"
import { variantApi, type VariantMetrics } from "@/lib/api"
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, PieChart, Pie, Cell } from "recharts"

interface VariantVisualizationProps {
  jobId: string
  onBack: () => void
}

export function VariantVisualization({ jobId, onBack }: VariantVisualizationProps) {
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [metrics, setMetrics] = useState<VariantMetrics | null>(null)
  const [normalizedChromosomes, setNormalizedChromosomes] = useState(false)
  const [proteinAlteringOnly, setProteinAlteringOnly] = useState(false)

  const loadMetrics = async () => {
    try {
      setLoading(true)
      setError(null)
      const data = await variantApi.getVariantMetrics(jobId)
      setMetrics(data)
    } catch (err: any) {
      setError(err.response?.data?.detail || "Failed to load variant metrics")
    } finally {
      setLoading(false)
    }
  }

  useEffect(() => {
    loadMetrics()
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobId])

  if (loading) {
    return (
      <div className="flex flex-col items-center justify-center h-96 space-y-4">
        <Loader2 className="h-8 w-8 animate-spin text-primary" />
        <p className="text-sm text-muted-foreground">Loading variant analysis...</p>
      </div>
    )
  }

  if (error) {
    return (
      <div className="flex flex-col items-center justify-center h-96 space-y-4">
        <p className="text-sm text-destructive">{error}</p>
        <Button onClick={onBack} variant="outline">
          <ArrowLeft className="mr-2 h-4 w-4" />
          Back to Jobs
        </Button>
      </div>
    )
  }

  if (!metrics || !metrics.metrics) return null

  const { headline, chromosome_distribution, chromosome_distribution_normalized, gene_distribution, gene_distribution_protein_altering, functional_impact } = metrics.metrics

  // Prepare chromosome data
  const chromosomeData = normalizedChromosomes ? chromosome_distribution_normalized : chromosome_distribution

  // Prepare gene data
  const geneData = proteinAlteringOnly ? gene_distribution_protein_altering : gene_distribution

  // Prepare functional impact data for pie chart
  const functionalCategoriesData = functional_impact?.categories
    ? Object.entries(functional_impact.categories)
        .filter(([_, count]) => count > 0)
        .map(([category, count]) => ({
          name: category.charAt(0).toUpperCase() + category.slice(1),
          value: count as number
        }))
    : []

  const exonicSubcategoriesData = functional_impact?.exonic_subcategories
    ? Object.entries(functional_impact.exonic_subcategories)
        .filter(([_, count]) => count > 0)
        .map(([category, count]) => ({
          name: category.charAt(0).toUpperCase() + category.slice(1),
          value: count as number
        }))
    : []

  // Colors for charts
  const CHROMOSOME_COLORS = {
    autosome: "#3b82f6", // blue
    sex: "#f59e0b" // amber
  }

  const PIE_COLORS = ["#3b82f6", "#8b5cf6", "#ec4899", "#f59e0b", "#10b981", "#6366f1", "#f43f5e"]

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4">
        <div>
          <h2 className="text-2xl sm:text-3xl font-bold">Variant Analysis</h2>
          <p className="text-sm text-muted-foreground mt-1">
            Sample: {metrics.sample_name}
          </p>
        </div>
        <Button onClick={onBack} variant="outline">
          <ArrowLeft className="mr-2 h-4 w-4" />
          Back to Jobs
        </Button>
      </div>

      {/* 1. Headline Metrics - KPI Cards */}
      <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-5 gap-4">
        <Card>
          <CardHeader className="pb-2">
            <CardDescription className="text-xs">Total Variants</CardDescription>
          </CardHeader>
          <CardContent>
            <p className="text-2xl font-bold">{headline.total_variants.toLocaleString()}</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2">
            <CardDescription className="text-xs">Total Genes</CardDescription>
          </CardHeader>
          <CardContent>
            <p className="text-2xl font-bold">{headline.total_genes.toLocaleString()}</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2">
            <CardDescription className="text-xs">SNVs</CardDescription>
          </CardHeader>
          <CardContent>
            <p className="text-2xl font-bold">{headline.snvs.toLocaleString()}</p>
            <p className="text-xs text-muted-foreground mt-1">
              {headline.total_variants > 0
                ? `${((headline.snvs / headline.total_variants) * 100).toFixed(1)}%`
                : "0%"}
            </p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2">
            <CardDescription className="text-xs">INDELs</CardDescription>
          </CardHeader>
          <CardContent>
            <p className="text-2xl font-bold">{headline.indels.toLocaleString()}</p>
            <p className="text-xs text-muted-foreground mt-1">
              {headline.total_variants > 0
                ? `${((headline.indels / headline.total_variants) * 100).toFixed(1)}%`
                : "0%"}
            </p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2">
            <CardDescription className="text-xs">High Confidence (PASS)</CardDescription>
          </CardHeader>
          <CardContent>
            <p className="text-2xl font-bold">{headline.high_confidence.toLocaleString()}</p>
            <p className="text-xs text-muted-foreground mt-1">
              {headline.total_variants > 0
                ? `${((headline.high_confidence / headline.total_variants) * 100).toFixed(1)}%`
                : "0%"}
            </p>
          </CardContent>
        </Card>
      </div>

      {/* 2. Chromosome-wise Distribution */}
      <Card>
        <CardHeader>
          <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
            <div>
              <CardTitle>Chromosome-wise Variant Distribution</CardTitle>
              <CardDescription className="mt-1">
                {normalizedChromosomes
                  ? "Variants per megabase (Mb) for each chromosome"
                  : "Total variant counts per chromosome (PASS variants only)"}
              </CardDescription>
            </div>
            <Button
              variant="outline"
              size="sm"
              onClick={() => setNormalizedChromosomes(!normalizedChromosomes)}
            >
              {normalizedChromosomes ? "Show Absolute Counts" : "Normalize by Length"}
            </Button>
          </div>
        </CardHeader>
        <CardContent>
          <ResponsiveContainer width="100%" height={400}>
            <BarChart data={chromosomeData} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis
                dataKey="chromosome"
                angle={-45}
                textAnchor="end"
                height={80}
                tick={{ fontSize: 12 }}
              />
              <YAxis
                label={{
                  value: normalizedChromosomes ? "Variants per Mb" : "Variant Count",
                  angle: -90,
                  position: "insideLeft"
                }}
                tick={{ fontSize: 12 }}
              />
              <Tooltip
                content={({ active, payload }) => {
                  if (active && payload && payload.length) {
                    const data = payload[0].payload
                    return (
                      <div className="bg-background border rounded p-2 shadow-lg">
                        <p className="font-semibold">{data.chromosome}</p>
                        <p className="text-sm">
                          {normalizedChromosomes
                            ? `${data.count.toFixed(2)} variants/Mb`
                            : `${data.count} variants`}
                        </p>
                        <p className="text-xs text-muted-foreground capitalize">{data.type}</p>
                      </div>
                    )
                  }
                  return null
                }}
              />
              <Legend
                wrapperStyle={{ paddingTop: "20px" }}
                content={() => (
                  <div className="flex justify-center gap-4 text-sm">
                    <div className="flex items-center gap-2">
                      <div className="w-3 h-3" style={{ backgroundColor: CHROMOSOME_COLORS.autosome }} />
                      <span>Autosome</span>
                    </div>
                    <div className="flex items-center gap-2">
                      <div className="w-3 h-3" style={{ backgroundColor: CHROMOSOME_COLORS.sex }} />
                      <span>Sex Chromosome</span>
                    </div>
                  </div>
                )}
              />
              <Bar dataKey="count" fill="#3b82f6">
                {chromosomeData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={CHROMOSOME_COLORS[entry.type]} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </CardContent>
      </Card>

      {/* 3. Gene-wise Distribution */}
      <Card>
        <CardHeader>
          <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
            <div>
              <CardTitle>Gene-wise Variant Distribution</CardTitle>
              <CardDescription className="mt-1">
                {proteinAlteringOnly
                  ? "Top genes with protein-altering variants (missense, nonsense, frameshift, splicing)"
                  : "Top genes by total variant count"}
              </CardDescription>
            </div>
            <Button
              variant="outline"
              size="sm"
              onClick={() => setProteinAlteringOnly(!proteinAlteringOnly)}
            >
              {proteinAlteringOnly ? "Show All Variants" : "Protein-Altering Only"}
            </Button>
          </div>
        </CardHeader>
        <CardContent>
          <ResponsiveContainer width="100%" height={400}>
            <BarChart data={geneData} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis
                dataKey="gene"
                angle={-45}
                textAnchor="end"
                height={80}
                tick={{ fontSize: 12 }}
              />
              <YAxis
                label={{ value: "Variant Count", angle: -90, position: "insideLeft" }}
                tick={{ fontSize: 12 }}
              />
              <Tooltip
                content={({ active, payload }) => {
                  if (active && payload && payload.length) {
                    const data = payload[0].payload
                    return (
                      <div className="bg-background border rounded p-2 shadow-lg">
                        <p className="font-semibold">{data.gene}</p>
                        <p className="text-sm">{data.count} variants</p>
                      </div>
                    )
                  }
                  return null
                }}
              />
              <Bar dataKey="count" fill="#8b5cf6" />
            </BarChart>
          </ResponsiveContainer>
        </CardContent>
      </Card>

      {/* 4. Functional Impact Snapshot */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Main Functional Categories */}
        <Card>
          <CardHeader>
            <CardTitle>Functional Categories</CardTitle>
            <CardDescription>Distribution across genomic regions</CardDescription>
          </CardHeader>
          <CardContent>
            <ResponsiveContainer width="100%" height={350}>
              <PieChart>
                <Pie
                  data={functionalCategoriesData}
                  cx="50%"
                  cy="50%"
                  labelLine={false}
                  label={({ name, percent }) => `${name}: ${((percent ?? 0) * 100).toFixed(1)}%`}
                  outerRadius={100}
                  fill="#8884d8"
                  dataKey="value"
                >
                  {functionalCategoriesData.map((entry, index) => (
                    <Cell key={`cell-${index}`} fill={PIE_COLORS[index % PIE_COLORS.length]} />
                  ))}
                </Pie>
                <Tooltip
                  content={({ active, payload }) => {
                    if (active && payload && payload.length) {
                      const data = payload[0]
                      return (
                        <div className="bg-background border rounded p-2 shadow-lg">
                          <p className="font-semibold">{data.name}</p>
                          <p className="text-sm">{data.value} variants</p>
                        </div>
                      )
                    }
                    return null
                  }}
                />
              </PieChart>
            </ResponsiveContainer>
          </CardContent>
        </Card>

        {/* Exonic Subcategories */}
        <Card>
          <CardHeader>
            <CardTitle>Exonic Variant Types</CardTitle>
            <CardDescription>Breakdown of coding sequence variants</CardDescription>
          </CardHeader>
          <CardContent>
            <ResponsiveContainer width="100%" height={350}>
              <PieChart>
                <Pie
                  data={exonicSubcategoriesData}
                  cx="50%"
                  cy="50%"
                  labelLine={false}
                  label={({ name, percent }) => `${name}: ${((percent ?? 0) * 100).toFixed(1)}%`}
                  outerRadius={100}
                  fill="#8884d8"
                  dataKey="value"
                >
                  {exonicSubcategoriesData.map((entry, index) => (
                    <Cell key={`cell-${index}`} fill={PIE_COLORS[index % PIE_COLORS.length]} />
                  ))}
                </Pie>
                <Tooltip
                  content={({ active, payload }) => {
                    if (active && payload && payload.length) {
                      const data = payload[0]
                      return (
                        <div className="bg-background border rounded p-2 shadow-lg">
                          <p className="font-semibold">{data.name}</p>
                          <p className="text-sm">{data.value} variants</p>
                        </div>
                      )
                    }
                    return null
                  }}
                />
              </PieChart>
            </ResponsiveContainer>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
