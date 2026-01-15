'use client';

import { useState } from 'react';
import { Card, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Mail, Users, Microscope, Building2, Clock, CheckCircle2 } from 'lucide-react';

export default function ContactPage() {
  const [formData, setFormData] = useState({
    name: '',
    email: '',
    inquiryType: 'General Inquiry',
    message: ''
  });
  const [submitted, setSubmitted] = useState(false);

  const handleSubmit = () => {
    if (formData.name && formData.email && formData.message) {
      setSubmitted(true);
      setTimeout(() => setSubmitted(false), 4000);
      setFormData({ name: '', email: '', inquiryType: 'General Inquiry', message: '' });
    }
  };

  return (
    <div className="min-h-screen">
      {/* Hero Section */}
      <section className="bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark py-16 overflow-hidden">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <Badge variant="cyan" className="mb-6">
            Get in Touch
          </Badge>
          <h1 className="text-4xl sm:text-5xl lg:text-5xl font-bold text-white mb-6 leading-tight">
            Get in Touch with Our Team
          </h1>
          <p className="text-xl text-blue-100 leading-relaxed max-w-2xl mx-auto">
            We're here to support your genomic research and answer any questions about our platform
          </p>
        </div>
      </section>

      {/* Main Content */}
      <section className="py-24 bg-white">
        <div className="max-w-6xl mx-auto px-4 sm:px-6 lg:px-8">
          {/* Contact Methods Grid */}
          <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-6 mb-16">
            {[
              {
                icon: Mail,
                title: 'General Inquiries',
                content: 'aiimsgenomics@gmail.com',
                color: 'purple'
              },
              {
                icon: Users,
                title: 'Technical Support',
                content: 'aiimsgenomics@gmail.com',
                color: 'cyan'
              },
              {
                icon: Microscope,
                title: 'Research Collaboration',
                content: 'drprabudhgoel@gmail.com',
                color: 'teal'
              },
              {
                icon: Building2,
                title: 'Platform Access',
                content: 'aiimsgenomics@gmail.com',
                color: 'purple'
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
                  <CardContent className="p-6 text-center">
                    <div className="flex justify-center mb-4">
                      <div className={`w-14 h-14 rounded-xl flex items-center justify-center ${colorClasses[item.color as keyof typeof colorClasses].split(' ').slice(0, 2).join(' ')}`}>
                        <IconComponent className="h-7 w-7" />
                      </div>
                    </div>
                    <h3 className="font-bold text-slate-900 mb-2 text-sm">{item.title}</h3>
                    <p className="text-slate-600 text-xs break-words">{item.content}</p>
                  </CardContent>
                </Card>
              );
            })}
          </div>

          {/* Contact Form */}
          <div className="grid lg:grid-cols-3 gap-12">
            <div className="lg:col-span-2">
              <Card className="border-2 border-slate-200">
                <CardContent className="p-10">
                  <h2 className="text-2xl font-bold text-slate-900 mb-2">Send Us a Message</h2>
                  <p className="text-slate-600 mb-8">Fill out the form below and we'll get back to you as soon as possible</p>

                  {submitted && (
                    <div className="mb-8 p-5 bg-cyan/5 border-2 border-cyan/20 rounded-xl flex items-start gap-3">
                      <CheckCircle2 className="h-6 w-6 text-cyan flex-shrink-0 mt-0.5" />
                      <div>
                        <p className="text-cyan font-semibold text-lg">Message sent successfully!</p>
                        <p className="text-slate-600 text-sm mt-1">We'll get back to you within 24-48 hours.</p>
                      </div>
                    </div>
                  )}

                  <div className="space-y-6">
                    <div className="grid md:grid-cols-2 gap-6">
                      <div>
                        <label className="block text-sm font-semibold text-slate-700 mb-2">Name *</label>
                        <input
                          type="text"
                          value={formData.name}
                          onChange={(e) => setFormData({ ...formData, name: e.target.value })}
                          className="w-full px-4 py-3 border-2 border-slate-200 rounded-lg focus:outline-none focus:border-cyan focus:ring-2 focus:ring-cyan/20 transition-all"
                          placeholder="Your full name"
                        />
                      </div>
                      <div>
                        <label className="block text-sm font-semibold text-slate-700 mb-2">Email *</label>
                        <input
                          type="email"
                          value={formData.email}
                          onChange={(e) => setFormData({ ...formData, email: e.target.value })}
                          className="w-full px-4 py-3 border-2 border-slate-200 rounded-lg focus:outline-none focus:border-cyan focus:ring-2 focus:ring-cyan/20 transition-all"
                          placeholder="your@email.com"
                        />
                      </div>
                    </div>

                    <div>
                      <label className="block text-sm font-semibold text-slate-700 mb-2">Inquiry Type *</label>
                      <select
                        value={formData.inquiryType}
                        onChange={(e) => setFormData({ ...formData, inquiryType: e.target.value })}
                        className="w-full px-4 py-3 border-2 border-slate-200 rounded-lg focus:outline-none focus:border-cyan focus:ring-2 focus:ring-cyan/20 transition-all bg-white"
                      >
                        <option>General Inquiry</option>
                        <option>Technical Support</option>
                        <option>Sales/Pricing</option>
                        <option>Research Collaboration</option>
                        <option>Bug Report</option>
                        <option>Feature Request</option>
                      </select>
                    </div>

                    <div>
                      <label className="block text-sm font-semibold text-slate-700 mb-2">Message *</label>
                      <textarea
                        rows={8}
                        value={formData.message}
                        onChange={(e) => setFormData({ ...formData, message: e.target.value })}
                        className="w-full px-4 py-3 border-2 border-slate-200 rounded-lg focus:outline-none focus:border-cyan focus:ring-2 focus:ring-cyan/20 transition-all resize-none"
                        placeholder="Tell us more about your inquiry..."
                      />
                    </div>

                    <Button
                      onClick={handleSubmit}
                      className="w-full bg-purple-primary hover:bg-purple-light text-white py-6 text-base font-semibold"
                      disabled={!formData.name || !formData.email || !formData.message}
                    >
                      Send Message
                    </Button>
                  </div>
                </CardContent>
              </Card>
            </div>

            {/* Sidebar */}
            <div className="space-y-6">
              {/* Response Time */}
              <Card className="border-2 border-cyan/20 bg-gradient-to-br from-cyan/5 to-white">
                <CardContent className="p-8">
                  <div className="w-12 h-12 bg-cyan/10 rounded-xl flex items-center justify-center mb-4">
                    <Clock className="h-6 w-6 text-cyan" />
                  </div>
                  <h3 className="text-lg font-bold text-slate-900 mb-3">Response Time</h3>
                  <div className="space-y-3 text-sm">
                    <div className="flex justify-between">
                      <span className="text-slate-600">General Inquiries</span>
                      <span className="font-semibold text-slate-900">24-48h</span>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-slate-600">Technical Support</span>
                      <span className="font-semibold text-slate-900">12-24h</span>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-slate-600">Urgent Issues</span>
                      <span className="font-semibold text-slate-900">4-8h</span>
                    </div>
                  </div>
                  <p className="text-xs text-slate-500 mt-4 pt-4 border-t border-slate-200">
                    Response times are for business days (Mon-Fri)
                  </p>
                </CardContent>
              </Card>

              {/* Research Collaboration */}
              <Card className="border-2 border-purple-primary/20 bg-gradient-to-br from-purple-primary/5 to-white">
                <CardContent className="p-8">
                  <div className="w-12 h-12 bg-purple-primary/10 rounded-xl flex items-center justify-center mb-4">
                    <Microscope className="h-6 w-6 text-purple-primary" />
                  </div>
                  <h3 className="text-lg font-bold text-slate-900 mb-3">Research Collaborations</h3>
                  <p className="text-slate-600 text-sm leading-relaxed mb-4">
                    We welcome collaborations with research institutions and clinical laboratories.
                  </p>
                  <ul className="space-y-2 text-sm text-slate-600">
                    <li className="flex items-start gap-2">
                      <div className="w-1.5 h-1.5 bg-purple-primary rounded-full mt-1.5 flex-shrink-0"></div>
                      <span>Joint research projects</span>
                    </li>
                    <li className="flex items-start gap-2">
                      <div className="w-1.5 h-1.5 bg-purple-primary rounded-full mt-1.5 flex-shrink-0"></div>
                      <span>Platform validation studies</span>
                    </li>
                    <li className="flex items-start gap-2">
                      <div className="w-1.5 h-1.5 bg-purple-primary rounded-full mt-1.5 flex-shrink-0"></div>
                      <span>Custom feature development</span>
                    </li>
                  </ul>
                </CardContent>
              </Card>

              {/* Quick Tips */}
              <Card className="border-2 border-teal/20 bg-gradient-to-br from-teal/5 to-white">
                <CardContent className="p-6">
                  <h3 className="text-sm font-bold text-slate-900 mb-3">Quick Tips</h3>
                  <ul className="space-y-2 text-xs text-slate-600">
                    <li className="flex items-start gap-2">
                      <span className="text-teal">•</span>
                      <span>Include your institution name for research inquiries</span>
                    </li>
                    <li className="flex items-start gap-2">
                      <span className="text-teal">•</span>
                      <span>Attach screenshots for technical issues</span>
                    </li>
                    <li className="flex items-start gap-2">
                      <span className="text-teal">•</span>
                      <span>Check our documentation before submitting</span>
                    </li>
                  </ul>
                </CardContent>
              </Card>
            </div>
          </div>
        </div>
      </section>
    </div>
  );
}
