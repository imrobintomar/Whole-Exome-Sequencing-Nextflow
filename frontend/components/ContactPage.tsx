'use client';

import { useState } from 'react';
import { Card, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Mail, Users, Microscope, Building2, Clock, CheckCircle2, AlertCircle } from 'lucide-react';

interface FormErrors {
  name?: string;
  email?: string;
  message?: string;
}

// Email configuration - can be moved to environment variables
const CONTACT_EMAILS = {
  general: process.env.NEXT_PUBLIC_CONTACT_EMAIL || 'aiimsgenomics@gmail.com',
  support: process.env.NEXT_PUBLIC_SUPPORT_EMAIL || 'aiimsgenomics@gmail.com',
  collaboration: process.env.NEXT_PUBLIC_COLLABORATION_EMAIL || 'drprabudhgoel@gmail.com',
  access: process.env.NEXT_PUBLIC_ACCESS_EMAIL || 'aiimsgenomics@gmail.com'
};

export default function ContactPage() {
  const [formData, setFormData] = useState({
    name: '',
    email: '',
    inquiryType: 'General Inquiry',
    message: ''
  });
  const [submitted, setSubmitted] = useState(false);
  const [errors, setErrors] = useState<FormErrors>({});
  const [touched, setTouched] = useState<{ [key: string]: boolean }>({});

  // Email validation function
  const isValidEmail = (email: string): boolean => {
    const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
    return emailRegex.test(email);
  };

  // Validate single field
  const validateField = (name: string, value: string): string | undefined => {
    switch (name) {
      case 'name':
        if (!value.trim()) return 'Name is required';
        if (value.trim().length < 2) return 'Name must be at least 2 characters';
        return undefined;
      case 'email':
        if (!value.trim()) return 'Email is required';
        if (!isValidEmail(value)) return 'Please enter a valid email address';
        return undefined;
      case 'message':
        if (!value.trim()) return 'Message is required';
        if (value.trim().length < 10) return 'Message must be at least 10 characters';
        return undefined;
      default:
        return undefined;
    }
  };

  // Validate all fields
  const validateForm = (): boolean => {
    const newErrors: FormErrors = {
      name: validateField('name', formData.name),
      email: validateField('email', formData.email),
      message: validateField('message', formData.message)
    };

    setErrors(newErrors);
    return !newErrors.name && !newErrors.email && !newErrors.message;
  };

  const handleBlur = (field: string) => {
    setTouched({ ...touched, [field]: true });
    const error = validateField(field, formData[field as keyof typeof formData]);
    setErrors({ ...errors, [field]: error });
  };

  const handleChange = (field: string, value: string) => {
    setFormData({ ...formData, [field]: value });
    // Clear error when user starts typing
    if (touched[field]) {
      const error = validateField(field, value);
      setErrors({ ...errors, [field]: error });
    }
  };

  const handleSubmit = () => {
    // Mark all fields as touched
    setTouched({ name: true, email: true, message: true });

    if (validateForm()) {
      setSubmitted(true);
      setTimeout(() => setSubmitted(false), 4000);
      setFormData({ name: '', email: '', inquiryType: 'General Inquiry', message: '' });
      setErrors({});
      setTouched({});
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
                content: CONTACT_EMAILS.general,
                color: 'purple'
              },
              {
                icon: Users,
                title: 'Technical Support',
                content: CONTACT_EMAILS.support,
                color: 'cyan'
              },
              {
                icon: Microscope,
                title: 'Research Collaboration',
                content: CONTACT_EMAILS.collaboration,
                color: 'teal'
              },
              {
                icon: Building2,
                title: 'Platform Access',
                content: CONTACT_EMAILS.access,
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
                    <div className="mb-8 p-5 bg-cyan/5 border-2 border-cyan/20 rounded-xl flex items-start gap-3" role="alert" aria-live="polite">
                      <CheckCircle2 className="h-6 w-6 text-cyan flex-shrink-0 mt-0.5" />
                      <div>
                        <p className="text-cyan font-semibold text-lg">Message sent successfully!</p>
                        <p className="text-slate-600 text-sm mt-1">We'll get back to you within 24-48 hours.</p>
                      </div>
                    </div>
                  )}

                  <div className="space-y-6">
                    <div className="grid md:grid-cols-2 gap-6">
                      {/* Name Field */}
                      <div>
                        <label htmlFor="name" className="block text-sm font-semibold text-slate-700 mb-2">
                          Name <span className="text-red-500">*</span>
                        </label>
                        <input
                          id="name"
                          type="text"
                          value={formData.name}
                          onChange={(e) => handleChange('name', e.target.value)}
                          onBlur={() => handleBlur('name')}
                          className={`w-full px-4 py-3 border-2 rounded-lg focus:outline-none focus:ring-2 focus:ring-cyan/20 transition-all ${
                            errors.name && touched.name
                              ? 'border-red-500 focus:border-red-500'
                              : 'border-slate-200 focus:border-cyan'
                          }`}
                          placeholder="Your full name"
                          aria-required="true"
                          aria-invalid={!!(errors.name && touched.name)}
                          aria-describedby={errors.name && touched.name ? 'name-error' : undefined}
                        />
                        {errors.name && touched.name && (
                          <div id="name-error" className="flex items-center gap-2 mt-2 text-red-600 text-sm" role="alert">
                            <AlertCircle className="h-4 w-4 flex-shrink-0" />
                            <span>{errors.name}</span>
                          </div>
                        )}
                      </div>

                      {/* Email Field */}
                      <div>
                        <label htmlFor="email" className="block text-sm font-semibold text-slate-700 mb-2">
                          Email <span className="text-red-500">*</span>
                        </label>
                        <input
                          id="email"
                          type="email"
                          value={formData.email}
                          onChange={(e) => handleChange('email', e.target.value)}
                          onBlur={() => handleBlur('email')}
                          className={`w-full px-4 py-3 border-2 rounded-lg focus:outline-none focus:ring-2 focus:ring-cyan/20 transition-all ${
                            errors.email && touched.email
                              ? 'border-red-500 focus:border-red-500'
                              : 'border-slate-200 focus:border-cyan'
                          }`}
                          placeholder="your@email.com"
                          aria-required="true"
                          aria-invalid={!!(errors.email && touched.email)}
                          aria-describedby={errors.email && touched.email ? 'email-error' : undefined}
                        />
                        {errors.email && touched.email && (
                          <div id="email-error" className="flex items-center gap-2 mt-2 text-red-600 text-sm" role="alert">
                            <AlertCircle className="h-4 w-4 flex-shrink-0" />
                            <span>{errors.email}</span>
                          </div>
                        )}
                      </div>
                    </div>

                    {/* Inquiry Type Field */}
                    <div>
                      <label htmlFor="inquiryType" className="block text-sm font-semibold text-slate-700 mb-2">
                        Inquiry Type <span className="text-red-500">*</span>
                      </label>
                      <select
                        id="inquiryType"
                        value={formData.inquiryType}
                        onChange={(e) => handleChange('inquiryType', e.target.value)}
                        className="w-full px-4 py-3 border-2 border-slate-200 rounded-lg focus:outline-none focus:border-cyan focus:ring-2 focus:ring-cyan/20 transition-all bg-white"
                        aria-required="true"
                      >
                        <option>General Inquiry</option>
                        <option>Technical Support</option>
                        <option>Sales/Pricing</option>
                        <option>Research Collaboration</option>
                        <option>Bug Report</option>
                        <option>Feature Request</option>
                      </select>
                    </div>

                    {/* Message Field */}
                    <div>
                      <label htmlFor="message" className="block text-sm font-semibold text-slate-700 mb-2">
                        Message <span className="text-red-500">*</span>
                      </label>
                      <textarea
                        id="message"
                        rows={8}
                        value={formData.message}
                        onChange={(e) => handleChange('message', e.target.value)}
                        onBlur={() => handleBlur('message')}
                        className={`w-full px-4 py-3 border-2 rounded-lg focus:outline-none focus:ring-2 focus:ring-cyan/20 transition-all resize-none ${
                          errors.message && touched.message
                            ? 'border-red-500 focus:border-red-500'
                            : 'border-slate-200 focus:border-cyan'
                        }`}
                        placeholder="Tell us more about your inquiry..."
                        aria-required="true"
                        aria-invalid={!!(errors.message && touched.message)}
                        aria-describedby={errors.message && touched.message ? 'message-error' : undefined}
                      />
                      {errors.message && touched.message && (
                        <div id="message-error" className="flex items-center gap-2 mt-2 text-red-600 text-sm" role="alert">
                          <AlertCircle className="h-4 w-4 flex-shrink-0" />
                          <span>{errors.message}</span>
                        </div>
                      )}
                    </div>

                    <Button
                      onClick={handleSubmit}
                      className="w-full bg-purple-primary hover:bg-purple-light text-white py-6 text-base font-semibold"
                      type="submit"
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
