'use client';

import { useState, useEffect, useRef } from 'react';
import { chatApi } from '../lib/api';

interface AdminConversation {
  id: number;
  user_email: string;
  subject: string;
  status: string;
  last_message_at: string;
  unread_count?: number;
}

interface AdminMessage {
  id: number;
  message: string;
  sender_role: string;
  created_at: string;
}

export default function AdminChatPanel() {
  const [conversations, setConversations] = useState<AdminConversation[]>([]);
  const [selectedConversation, setSelectedConversation] = useState<number | null>(null);
  const [messages, setMessages] = useState<AdminMessage[]>([]);
  const [newMessage, setNewMessage] = useState('');
  const [loading, setLoading] = useState(true);
  const [sending, setSending] = useState(false);
  const [statusFilter, setStatusFilter] = useState<'all' | 'open' | 'resolved' | 'closed'>('open');
  const messagesEndRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    loadConversations();
    const interval = setInterval(loadConversations, 15000);
    return () => clearInterval(interval);
  }, [statusFilter]);

  useEffect(() => {
    if (selectedConversation) {
      loadMessages(selectedConversation);
      const interval = setInterval(() => loadMessages(selectedConversation), 10000);
      return () => clearInterval(interval);
    }
  }, [selectedConversation]);

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  };

  const loadConversations = async () => {
    try {
      const params = statusFilter !== 'all' ? { status: statusFilter } : undefined;
      const data = await chatApi.adminGetConversations(params);
      setConversations(data.conversations);
    } catch (error) {
      console.error('Failed to load conversations:', error);
    } finally {
      setLoading(false);
    }
  };

  const loadMessages = async (conversationId: number) => {
    try {
      const data = await chatApi.adminGetConversationMessages(conversationId);
      setMessages(data.messages);
    } catch (error) {
      console.error('Failed to load messages:', error);
    }
  };

  const handleSendMessage = async () => {
    if (!newMessage.trim() || !selectedConversation) return;

    setSending(true);
    try {
      await chatApi.adminSendMessage(selectedConversation, newMessage);
      setNewMessage('');
      await loadMessages(selectedConversation);
    } catch (error) {
      alert('Failed to send message. Please try again.');
    } finally {
      setSending(false);
    }
  };

  const handleUpdateStatus = async (conversationId: number, status: 'open' | 'resolved' | 'closed') => {
    try {
      await chatApi.adminUpdateStatus(conversationId, status);
      await loadConversations();
      if (selectedConversation === conversationId) {
        await loadMessages(conversationId);
      }
    } catch (error) {
      alert('Failed to update status. Please try again.');
    }
  };

  const selectedConv = conversations.find((c) => c.id === selectedConversation);

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-[400px]">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-indigo-600"></div>
      </div>
    );
  }

  return (
    <div className="grid md:grid-cols-3 gap-6 h-[calc(100vh-300px)]">
      {/* Conversations List */}
      <div className="bg-white rounded-xl shadow-lg overflow-hidden flex flex-col border border-gray-100">
        <div className="p-5 bg-gradient-to-r from-indigo-600 to-purple-600">
          <h2 className="font-bold text-white text-lg mb-3">Support Tickets</h2>
          <div className="flex gap-2 flex-wrap">
            {['all', 'open', 'resolved', 'closed'].map((filter) => (
              <button
                key={filter}
                onClick={() => setStatusFilter(filter as any)}
                className={`px-3 py-1.5 text-xs font-semibold rounded-lg transition-all ${
                  statusFilter === filter
                    ? 'bg-white text-indigo-700 shadow-md'
                    : 'bg-white/20 text-white hover:bg-white/30'
                }`}
              >
                {filter.charAt(0).toUpperCase() + filter.slice(1)}
              </button>
            ))}
          </div>
        </div>

        <div className="flex-1 overflow-y-auto">
          {conversations.length === 0 ? (
            <div className="p-6 text-center">
              <div className="bg-gray-50 rounded-lg p-8">
                <svg className="h-12 w-12 mx-auto mb-3 text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M20 13V6a2 2 0 00-2-2H6a2 2 0 00-2 2v7m16 0v5a2 2 0 01-2 2H6a2 2 0 01-2-2v-5m16 0h-2.586a1 1 0 00-.707.293l-2.414 2.414a1 1 0 01-.707.293h-3.172a1 1 0 01-.707-.293l-2.414-2.414A1 1 0 006.586 13H4" />
                </svg>
                <p className="text-gray-600 text-sm font-medium">No {statusFilter !== 'all' && statusFilter} tickets</p>
              </div>
            </div>
          ) : (
            conversations.map((conv) => (
              <button
                key={conv.id}
                onClick={() => setSelectedConversation(conv.id)}
                className={`w-full p-4 text-left border-b border-gray-100 hover:bg-gradient-to-r hover:from-indigo-50 hover:to-purple-50 transition-all ${
                  selectedConversation === conv.id
                    ? 'bg-gradient-to-r from-indigo-50 to-purple-50 border-l-4 border-l-indigo-600'
                    : ''
                }`}
              >
                <div className="flex items-start justify-between mb-2">
                  <div className="flex-1 min-w-0">
                    <h3 className="font-semibold text-gray-900 text-sm truncate">{conv.subject}</h3>
                    <div className="flex items-center gap-2 mt-1">
                      <div className="flex items-center gap-1">
                        <svg className="h-3 w-3 text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M16 7a4 4 0 11-8 0 4 4 0 018 0zM12 14a7 7 0 00-7 7h14a7 7 0 00-7-7z" />
                        </svg>
                        <p className="text-xs text-gray-600 truncate">{conv.user_email}</p>
                      </div>
                    </div>
                  </div>
                  <span
                    className={`px-2.5 py-1 text-xs font-bold rounded-full flex-shrink-0 ml-2 ${
                      conv.status === 'open'
                        ? 'bg-green-100 text-green-700'
                        : conv.status === 'resolved'
                        ? 'bg-blue-100 text-blue-700'
                        : 'bg-gray-100 text-gray-700'
                    }`}
                  >
                    {conv.status}
                  </span>
                </div>
                <p className="text-xs text-gray-500 flex items-center gap-1">
                  <svg className="h-3 w-3" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
                  </svg>
                  {new Date(conv.last_message_at).toLocaleString()}
                </p>
              </button>
            ))
          )}
        </div>
      </div>

      {/* Messages Area */}
      <div className="md:col-span-2 bg-white rounded-xl shadow-lg flex flex-col border border-gray-100 overflow-hidden">
        {selectedConversation && selectedConv ? (
          <>
            {/* Header */}
            <div className="p-5 bg-gradient-to-r from-indigo-600 to-purple-600 text-white">
              <div className="flex items-center justify-between mb-3">
                <div className="flex-1">
                  <h2 className="font-bold text-lg">{selectedConv.subject}</h2>
                  <div className="flex items-center gap-2 mt-1">
                    <svg className="h-4 w-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M16 7a4 4 0 11-8 0 4 4 0 018 0zM12 14a7 7 0 00-7 7h14a7 7 0 00-7-7z" />
                    </svg>
                    <p className="text-sm opacity-90">{selectedConv.user_email}</p>
                  </div>
                </div>
                <div className="flex gap-2">
                  {selectedConv.status !== 'resolved' && (
                    <button
                      onClick={() => handleUpdateStatus(selectedConversation, 'resolved')}
                      className="px-4 py-2 bg-white/20 hover:bg-white/30 rounded-lg text-sm font-medium transition-all backdrop-blur-sm"
                    >
                      ✓ Resolve
                    </button>
                  )}
                  {selectedConv.status !== 'closed' && (
                    <button
                      onClick={() => handleUpdateStatus(selectedConversation, 'closed')}
                      className="px-4 py-2 bg-white/20 hover:bg-white/30 rounded-lg text-sm font-medium transition-all backdrop-blur-sm"
                    >
                      Close
                    </button>
                  )}
                  {selectedConv.status === 'closed' && (
                    <button
                      onClick={() => handleUpdateStatus(selectedConversation, 'open')}
                      className="px-4 py-2 bg-white/20 hover:bg-white/30 rounded-lg text-sm font-medium transition-all backdrop-blur-sm"
                    >
                      ↻ Reopen
                    </button>
                  )}
                </div>
              </div>
              <span
                className={`inline-flex items-center gap-1.5 px-3 py-1 text-xs font-bold rounded-full ${
                  selectedConv.status === 'open'
                    ? 'bg-green-100 text-green-700'
                    : selectedConv.status === 'resolved'
                    ? 'bg-blue-100 text-blue-700'
                    : 'bg-gray-100 text-gray-700'
                }`}
              >
                <span className={`w-2 h-2 rounded-full ${
                  selectedConv.status === 'open' ? 'bg-green-500 animate-pulse' : 'bg-gray-400'
                }`}></span>
                {selectedConv.status}
              </span>
            </div>

            {/* Messages */}
            <div className="flex-1 overflow-y-auto p-6 space-y-4 bg-gradient-to-b from-gray-50 to-white">
              {messages.map((msg) => (
                <div
                  key={msg.id}
                  className={`flex items-end gap-2 ${msg.sender_role === 'admin' ? 'justify-end' : 'justify-start'}`}
                >
                  {msg.sender_role === 'user' && (
                    <div className="w-8 h-8 rounded-full bg-gradient-to-br from-gray-400 to-gray-600 flex items-center justify-center text-white font-bold text-xs flex-shrink-0">
                      {selectedConv.user_email[0].toUpperCase()}
                    </div>
                  )}
                  <div
                    className={`max-w-[70%] rounded-2xl px-4 py-3 shadow-md ${
                      msg.sender_role === 'admin'
                        ? 'bg-gradient-to-r from-indigo-600 to-purple-600 text-white rounded-br-sm'
                        : 'bg-white text-gray-900 border border-gray-200 rounded-bl-sm'
                    }`}
                  >
                    <p className="text-sm whitespace-pre-wrap leading-relaxed">{msg.message}</p>
                    <p
                      className={`text-xs mt-2 flex items-center gap-1 ${
                        msg.sender_role === 'admin' ? 'text-indigo-100' : 'text-gray-500'
                      }`}
                    >
                      {msg.sender_role === 'user' && <span className="font-medium">Customer •</span>}
                      {msg.sender_role === 'admin' && <span className="font-medium">You •</span>}
                      <span>{new Date(msg.created_at).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}</span>
                    </p>
                  </div>
                  {msg.sender_role === 'admin' && (
                    <div className="w-8 h-8 rounded-full bg-gradient-to-br from-indigo-600 to-purple-600 flex items-center justify-center text-white font-bold text-xs flex-shrink-0">
                      A
                    </div>
                  )}
                </div>
              ))}
              <div ref={messagesEndRef} />
            </div>

            {/* Input */}
            {selectedConv.status !== 'closed' && (
              <div className="p-4 border-t border-gray-200 bg-white">
                <div className="flex gap-3 items-end">
                  <div className="flex-1">
                    <textarea
                      value={newMessage}
                      onChange={(e) => setNewMessage(e.target.value)}
                      onKeyPress={(e) => {
                        if (e.key === 'Enter' && !e.shiftKey) {
                          e.preventDefault();
                          !sending && handleSendMessage();
                        }
                      }}
                      placeholder="Type your response... (Shift+Enter for new line)"
                      rows={2}
                      className="w-full px-4 py-3 border border-gray-300 rounded-xl focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500 resize-none"
                    />
                  </div>
                  <button
                    onClick={handleSendMessage}
                    disabled={sending || !newMessage.trim()}
                    className="bg-gradient-to-r from-indigo-600 to-purple-600 text-white px-6 py-3 rounded-xl hover:from-indigo-700 hover:to-purple-700 disabled:opacity-50 disabled:cursor-not-allowed font-medium shadow-lg transition-all flex items-center gap-2"
                  >
                    <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 19l9 2-9-18-9 18 9-2zm0 0v-8" />
                    </svg>
                    Send
                  </button>
                </div>
              </div>
            )}
          </>
        ) : (
          <div className="flex items-center justify-center h-full text-gray-500 bg-gradient-to-b from-gray-50 to-white">
            <div className="text-center">
              <div className="bg-gradient-to-br from-indigo-100 to-purple-100 w-20 h-20 rounded-full flex items-center justify-center mx-auto mb-4">
                <svg className="h-10 w-10 text-indigo-600" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path
                    strokeLinecap="round"
                    strokeLinejoin="round"
                    strokeWidth={2}
                    d="M8 12h.01M12 12h.01M16 12h.01M21 12c0 4.418-4.03 8-9 8a9.863 9.863 0 01-4.255-.949L3 20l1.395-3.72C3.512 15.042 3 13.574 3 12c0-4.418 4.03-8 9-8s9 3.582 9 8z"
                  />
                </svg>
              </div>
              <p className="text-gray-700 font-medium">Select a conversation</p>
              <p className="text-gray-500 text-sm mt-1">Choose a ticket from the list to view messages</p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
