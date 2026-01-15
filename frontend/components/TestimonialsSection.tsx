'use client';

import { useState } from 'react';
import { Card, CardContent } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { ChevronLeft, ChevronRight, Star, Quote } from 'lucide-react';

interface Testimonial {
  quote: string;
  author: string;
  title: string;
  institution: string;
  avatar?: string;
  rating: number;
  category: 'research' | 'clinical' | 'pharma' | 'academic';
}

const testimonials: Testimonial[] = [
  {
    quote: "ATGC Flow has transformed our rare disease diagnostic workflow. The automated ACMG classification saves our team 15+ hours per week, and the accuracy is outstanding. We've increased our diagnostic yield by 12% since implementation.",
    author: "Dr. Sarah Chen",
    title: "Director of Molecular Diagnostics",
    institution: "Children's Research Hospital",
    rating: 5,
    category: 'clinical'
  },
  {
    quote: "As a PhD student, I needed enterprise-grade tools on a student budget. ATGC Flow's flexible pricing let me analyze 200 exomes for my thesis. The IGV integration and automated reports made my publications so much easier.",
    author: "Michael Rodriguez",
    title: "PhD Candidate",
    institution: "Stanford University",
    rating: 5,
    category: 'academic'
  },
  {
    quote: "The pipeline validation metrics are exceptional. We compared ATGC Flow against our in-house system using GIAB samples and it outperformed in both sensitivity and specificity. The team's technical support has been invaluable.",
    author: "Dr. Jennifer Park",
    title: "Senior Bioinformatics Scientist",
    institution: "Genome Research Institute",
    rating: 5,
    category: 'research'
  },
  {
    quote: "We needed to stratify 500 patients for a phase II trial in under 2 months. ATGC Flow's API integration with our LIMS and batch processing capabilities made this possible. The trial enrolled ahead of schedule.",
    author: "Dr. Robert Thompson",
    title: "Head of Translational Genomics",
    institution: "BioPharma Inc.",
    rating: 5,
    category: 'pharma'
  },
  {
    quote: "The gene panel filtering is incredibly powerful. We can quickly focus on relevant variants for specific phenotypes using HPO terms. This has cut our analysis time in half and improved our diagnostic accuracy.",
    author: "Dr. Amanda Foster",
    title: "Clinical Geneticist",
    institution: "Metropolitan Medical Center",
    rating: 5,
    category: 'clinical'
  },
  {
    quote: "Running a population genomics study across 8 sites was a logistical nightmare until we standardized on ATGC Flow. The API allowed us to maintain pipeline consistency and we analyzed 5,000 exomes in record time.",
    author: "Prof. David Liu",
    title: "Principal Investigator",
    institution: "National Genomics Consortium",
    rating: 5,
    category: 'research'
  }
];

interface TestimonialsSectionProps {
  title?: string;
  subtitle?: string;
  showBadge?: boolean;
  badgeText?: string;
  maxVisible?: number;
  autoScroll?: boolean;
  backgroundColor?: 'white' | 'slate' | 'gradient';
}

export default function TestimonialsSection({
  title = "What Our Users Say",
  subtitle = "Trusted by researchers and clinicians worldwide",
  showBadge = true,
  badgeText = "Testimonials",
  maxVisible = 3,
  autoScroll = false,
  backgroundColor = 'white'
}: TestimonialsSectionProps) {
  const [currentIndex, setCurrentIndex] = useState(0);

  const nextSlide = () => {
    setCurrentIndex((prev) => (prev + 1) % Math.ceil(testimonials.length / maxVisible));
  };

  const prevSlide = () => {
    setCurrentIndex((prev) =>
      prev === 0 ? Math.ceil(testimonials.length / maxVisible) - 1 : prev - 1
    );
  };

  const visibleTestimonials = testimonials.slice(
    currentIndex * maxVisible,
    currentIndex * maxVisible + maxVisible
  );

  const bgClasses = {
    white: 'bg-white',
    slate: 'bg-gradient-to-br from-slate-50 to-white',
    gradient: 'bg-gradient-to-br from-purple-primary/5 via-cyan/5 to-white'
  };

  const categoryColors = {
    research: 'purple',
    clinical: 'cyan',
    pharma: 'teal',
    academic: 'purple'
  };

  return (
    <section className={`py-24 ${bgClasses[backgroundColor]}`}>
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        {/* Header */}
        <div className="text-center mb-16">
          {showBadge && (
            <Badge variant="outline" className="mb-4 border-purple-primary text-purple-primary">
              {badgeText}
            </Badge>
          )}
          <h2 className="text-4xl font-bold text-slate-900 mb-4">{title}</h2>
          <p className="text-xl text-slate-600 max-w-2xl mx-auto">{subtitle}</p>
        </div>

        {/* Testimonials Carousel */}
        <div className="relative">
          <div className={`grid ${maxVisible === 1 ? 'grid-cols-1' : maxVisible === 2 ? 'md:grid-cols-2' : 'lg:grid-cols-3'} gap-8`}>
            {visibleTestimonials.map((testimonial, index) => (
              <Card key={index} className="border-2 border-slate-200 hover:shadow-xl transition-all duration-300 relative">
                <CardContent className="p-8">
                  {/* Quote Icon */}
                  <div className="absolute top-6 right-6">
                    <Quote className="h-12 w-12 text-purple-primary/10" />
                  </div>

                  {/* Rating */}
                  <div className="flex gap-1 mb-4">
                    {[...Array(testimonial.rating)].map((_, i) => (
                      <Star key={i} className="h-5 w-5 fill-purple-primary text-purple-primary" />
                    ))}
                  </div>

                  {/* Quote */}
                  <p className="text-slate-700 leading-relaxed mb-6 italic">
                    "{testimonial.quote}"
                  </p>

                  {/* Author Info */}
                  <div className="border-t border-slate-200 pt-6">
                    <div className="flex items-start gap-4">
                      {/* Avatar */}
                      <div className="w-12 h-12 bg-gradient-to-br from-purple-primary to-cyan rounded-full flex items-center justify-center text-white font-bold flex-shrink-0">
                        {testimonial.author.split(' ').map(n => n[0]).join('')}
                      </div>

                      {/* Details */}
                      <div className="flex-1">
                        <div className="font-bold text-slate-900 mb-1">{testimonial.author}</div>
                        <div className="text-sm text-slate-600 mb-2">{testimonial.title}</div>
                        <div className="flex items-center gap-2 flex-wrap">
                          <span className="text-xs text-slate-500">{testimonial.institution}</span>
                          <Badge
                            variant={categoryColors[testimonial.category] as any}
                            className="text-xs"
                          >
                            {testimonial.category === 'research' ? 'Research' :
                             testimonial.category === 'clinical' ? 'Clinical' :
                             testimonial.category === 'pharma' ? 'Pharma' :
                             'Academic'}
                          </Badge>
                        </div>
                      </div>
                    </div>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>

          {/* Navigation Buttons */}
          {testimonials.length > maxVisible && (
            <div className="flex items-center justify-center gap-4 mt-12">
              <button
                onClick={prevSlide}
                className="w-12 h-12 rounded-full bg-purple-primary text-white hover:bg-purple-light transition-colors flex items-center justify-center"
                aria-label="Previous testimonials"
              >
                <ChevronLeft className="h-6 w-6" />
              </button>

              {/* Dots Indicator */}
              <div className="flex gap-2">
                {[...Array(Math.ceil(testimonials.length / maxVisible))].map((_, index) => (
                  <button
                    key={index}
                    onClick={() => setCurrentIndex(index)}
                    className={`h-2 rounded-full transition-all ${
                      index === currentIndex
                        ? 'w-8 bg-purple-primary'
                        : 'w-2 bg-slate-300 hover:bg-slate-400'
                    }`}
                    aria-label={`Go to testimonial set ${index + 1}`}
                  />
                ))}
              </div>

              <button
                onClick={nextSlide}
                className="w-12 h-12 rounded-full bg-purple-primary text-white hover:bg-purple-light transition-colors flex items-center justify-center"
                aria-label="Next testimonials"
              >
                <ChevronRight className="h-6 w-6" />
              </button>
            </div>
          )}
        </div>

        {/* Stats Summary */}
        <div className="mt-16 grid grid-cols-2 md:grid-cols-4 gap-8">
          {[
            { label: 'Average Rating', value: '5.0', icon: Star },
            { label: 'Happy Users', value: '500+', icon: null },
            { label: 'Research Papers', value: '200+', icon: null },
            { label: 'Clinical Labs', value: '50+', icon: null }
          ].map((stat, i) => {
            const IconComponent = stat.icon;
            return (
              <div key={i} className="text-center">
                <div className="text-4xl font-bold text-purple-primary mb-2 flex items-center justify-center gap-2">
                  {stat.value}
                  {IconComponent && <IconComponent className="h-6 w-6 fill-purple-primary" />}
                </div>
                <div className="text-sm text-slate-600 font-medium">{stat.label}</div>
              </div>
            );
          })}
        </div>
      </div>
    </section>
  );
}
