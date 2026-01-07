'use client';

import { useState } from 'react';
import { Card, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Globe, Users, Microscope, Mail } from 'lucide-react';

export default function ContactPage() {
  const [formData, setFormData] = useState({ name: '', email: '', subject: '', message: '' });
  const [submitted, setSubmitted] = useState(false);

  const handleSubmit = () => {
    if (formData.name && formData.email && formData.subject && formData.message) {
      setSubmitted(true);
      setTimeout(() => setSubmitted(false), 3000);
      setFormData({ name: '', email: '', subject: '', message: '' });
    }
  };

  return (
    <div className="max-w-2xl mx-auto px-4 sm:px-6 lg:px-8 py-16 sm:py-24">
      <div className="space-y-12">
        <div>
          <h1 className="text-4xl sm:text-5xl font-bold text-slate-900 mb-4">Contact Us</h1>
          <p className="text-lg text-slate-600 leading-relaxed">
            Have questions about our whole exome sequencing platform? We'd love to hear from you. Get in touch with our team.
          </p>
        </div>

        <div className="grid md:grid-cols-2 gap-8 mb-12">
          {[
            {
              icon: Mail,
              title: 'General Inquiries',
              content: 'aiimsgenomics@gmail.com'
            },
            {
              icon: Users,
              title: 'Technical Support',
              content: 'aiimsgenomics@gmail.com'
            },
            {
              icon: Microscope,
              title: 'Research Collaboration',
              content: 'aiimsgenomics@gmail.com'
            },
            {
              icon: Globe,
              title: 'Platform Access',
              content: 'aiimsgenomics@gmail.com'
            }
          ].map((item, i) => {
            const IconComponent = item.icon;
            return (
              <Card key={i} className="border-slate-200">
                <CardContent className="p-6 text-center">
                  <div className="flex justify-center mb-4">
                    <div className="w-12 h-12 bg-blue-100 rounded-lg flex items-center justify-center">
                      <IconComponent className="h-6 w-6 text-blue-600" />
                    </div>
                  </div>
                  <h3 className="font-bold text-slate-900 mb-2">{item.title}</h3>
                  <p className="text-slate-600 text-sm">{item.content}</p>
                </CardContent>
              </Card>
            );
          })}
        </div>

        <Card className="border-slate-200">
          <CardContent className="p-8">
            {submitted && (
              <div className="mb-6 p-4 bg-green-50 border border-green-200 rounded-lg">
                <p className="text-green-800 font-semibold">✓ Message sent successfully! We'll get back to you soon.</p>
              </div>
            )}

            <div className="space-y-6">
              <div className="grid md:grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-2">Name</label>
                  <input
                    type="text"
                    value={formData.name}
                    onChange={(e) => setFormData({ ...formData, name: e.target.value })}
                    className="w-full px-4 py-2 border border-slate-300 rounded-lg focus:outline-none focus:border-blue-500"
                    placeholder="Your name"
                  />
                </div>
                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-2">Email</label>
                  <input
                    type="email"
                    value={formData.email}
                    onChange={(e) => setFormData({ ...formData, email: e.target.value })}
                    className="w-full px-4 py-2 border border-slate-300 rounded-lg focus:outline-none focus:border-blue-500"
                    placeholder="your@email.com"
                  />
                </div>
              </div>

              <div>
                <label className="block text-sm font-medium text-slate-700 mb-2">Subject</label>
                <input
                  type="text"
                  value={formData.subject}
                  onChange={(e) => setFormData({ ...formData, subject: e.target.value })}
                  className="w-full px-4 py-2 border border-slate-300 rounded-lg focus:outline-none focus:border-blue-500"
                  placeholder="How can we help?"
                />
              </div>

              <div>
                <label className="block text-sm font-medium text-slate-700 mb-2">Message</label>
                <textarea
                  rows={6}
                  value={formData.message}
                  onChange={(e) => setFormData({ ...formData, message: e.target.value })}
                  className="w-full px-4 py-2 border border-slate-300 rounded-lg focus:outline-none focus:border-blue-500"
                  placeholder="Your message..."
                />
              </div>

              <Button
                onClick={handleSubmit}
                className="w-full bg-blue-600 hover:bg-blue-700 text-white py-2"
              >
                Send Message
              </Button>
            </div>
          </CardContent>
        </Card>

        <div className="bg-slate-50 border border-slate-200 rounded-lg p-8">
          <h2 className="text-lg font-bold text-slate-900 mb-4">Response Time</h2>
          <p className="text-slate-600">
            We typically respond to inquiries within 24-48 hours during business days. For urgent technical support requests, please mark your subject as "Urgent Support Request".
          </p>
        </div>

        <div className="bg-blue-50 border border-blue-200 rounded-lg p-8">
          <h2 className="text-lg font-bold text-blue-900 mb-4">Research Collaborations</h2>
          <p className="text-blue-800">
            We welcome collaborations with research institutions and clinical laboratories. If you're interested in using our platform for your genomics research or have ideas for improving variant analysis, please reach out to discuss potential partnerships.
          </p>
        </div>
      </div>
    </div>
  );
}
