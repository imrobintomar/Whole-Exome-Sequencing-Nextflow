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
  sender_type: string;
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
    // Auto-refresh every 15 seconds
    const interval = setInterval(loadConversations, 15000);
    return () => clearInterval(interval);
  }, [statusFilter]);

  useEffect(() => {
    if (selectedConversation) {
      loadMessages(selectedConversation);
      // Auto-refresh messages every 10 seconds
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
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600"></div>
      </div>
    );
  }

  return (
    <div className="grid md:grid-cols-3 gap-6 h-[calc(100vh-300px)]">
      {/* Conversations List */}
      <div className="bg-white rounded-lg shadow overflow-hidden flex flex-col">
        <div className="p-4 border-b border-gray-200">
          <h2 className="font-semibold mb-3">Support Tickets</h2>
          <div className="flex gap-2">
            {['all', 'open', 'resolved', 'closed'].map((filter) => (
              <button
                key={filter}
                onClick={() => setStatusFilter(filter as any)}
                className={`px-3 py-1 text-xs font-medium rounded-full ${
                  statusFilter === filter
                    ? 'bg-blue-600 text-white'
                    : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
                }`}
              >
                {filter}
              </button>
            ))}
          </div>
        </div>

        <div className="flex-1 overflow-y-auto">
          {conversations.length === 0 ? (
            <div className="p-4 text-center text-gray-600 text-sm">
              No {statusFilter !== 'all' && statusFilter} conversations
            </div>
          ) : (
            conversations.map((conv) => (
              <button
                key={conv.id}
                onClick={() => setSelectedConversation(conv.id)}
                className={`w-full p-4 text-left border-b border-gray-200 hover:bg-gray-50 transition-colors ${
                  selectedConversation === conv.id ? 'bg-blue-50' : ''
                }`}
              >
                <div className="flex items-start justify-between mb-1">
                  <div className="flex-1 min-w-0">
                    <h3 className="font-medium text-gray-900 text-sm truncate">{conv.subject}</h3>
                    <p className="text-xs text-gray-600 truncate">{conv.user_email}</p>
                  </div>
                  <span
                    className={`px-2 py-0.5 text-xs font-medium rounded-full flex-shrink-0 ml-2 ${
                      conv.status === 'open'
                        ? 'bg-green-100 text-green-800'
                        : conv.status === 'resolved'
                        ? 'bg-blue-100 text-blue-800'
                        : 'bg-gray-100 text-gray-800'
                    }`}
                  >
                    {conv.status}
                  </span>
                </div>
                <p className="text-xs text-gray-500">{new Date(conv.last_message_at).toLocaleString()}</p>
              </button>
            ))
          )}
        </div>
      </div>

      {/* Messages Area */}
      <div className="md:col-span-2 bg-white rounded-lg shadow flex flex-col">
        {selectedConversation && selectedConv ? (
          <>
            <div className="p-4 border-b border-gray-200">
              <div className="flex items-center justify-between mb-2">
                <div>
                  <h2 className="font-semibold">{selectedConv.subject}</h2>
                  <p className="text-xs text-gray-600">User: {selectedConv.user_email}</p>
                </div>
                <div className="flex gap-2">
                  {selectedConv.status !== 'resolved' && (
                    <button
                      onClick={() => handleUpdateStatus(selectedConversation, 'resolved')}
                      className="text-sm px-3 py-1 bg-blue-100 text-blue-700 rounded hover:bg-blue-200"
                    >
                      Mark Resolved
                    </button>
                  )}
                  {selectedConv.status !== 'closed' && (
                    <button
                      onClick={() => handleUpdateStatus(selectedConversation, 'closed')}
                      className="text-sm px-3 py-1 bg-gray-100 text-gray-700 rounded hover:bg-gray-200"
                    >
                      Close
                    </button>
                  )}
                  {selectedConv.status === 'closed' && (
                    <button
                      onClick={() => handleUpdateStatus(selectedConversation, 'open')}
                      className="text-sm px-3 py-1 bg-green-100 text-green-700 rounded hover:bg-green-200"
                    >
                      Reopen
                    </button>
                  )}
                </div>
              </div>
              <span
                className={`inline-block px-2 py-0.5 text-xs font-medium rounded-full ${
                  selectedConv.status === 'open'
                    ? 'bg-green-100 text-green-800'
                    : selectedConv.status === 'resolved'
                    ? 'bg-blue-100 text-blue-800'
                    : 'bg-gray-100 text-gray-800'
                }`}
              >
                Status: {selectedConv.status}
              </span>
            </div>

            <div className="flex-1 overflow-y-auto p-4 space-y-4">
              {messages.map((msg) => (
                <div
                  key={msg.id}
                  className={`flex ${msg.sender_type === 'admin' ? 'justify-end' : 'justify-start'}`}
                >
                  <div
                    className={`max-w-[70%] rounded-lg p-3 ${
                      msg.sender_type === 'admin'
                        ? 'bg-blue-600 text-white'
                        : 'bg-gray-100 text-gray-900'
                    }`}
                  >
                    <p className="text-sm whitespace-pre-wrap">{msg.message}</p>
                    <p
                      className={`text-xs mt-1 ${
                        msg.sender_type === 'admin' ? 'text-blue-100' : 'text-gray-500'
                      }`}
                    >
                      {msg.sender_type === 'user' && 'Customer • '}
                      {new Date(msg.created_at).toLocaleTimeString()}
                    </p>
                  </div>
                </div>
              ))}
              <div ref={messagesEndRef} />
            </div>

            {selectedConv.status !== 'closed' && (
              <div className="p-4 border-t border-gray-200">
                <div className="flex gap-2">
                  <input
                    type="text"
                    value={newMessage}
                    onChange={(e) => setNewMessage(e.target.value)}
                    onKeyPress={(e) => e.key === 'Enter' && !sending && handleSendMessage()}
                    placeholder="Type your response..."
                    className="flex-1 px-3 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-blue-500"
                  />
                  <button
                    onClick={handleSendMessage}
                    disabled={sending || !newMessage.trim()}
                    className="bg-blue-600 text-white px-4 py-2 rounded-lg hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
                  >
                    Send
                  </button>
                </div>
              </div>
            )}
          </>
        ) : (
          <div className="flex items-center justify-center h-full text-gray-500">
            <div className="text-center">
              <svg className="h-16 w-16 mx-auto mb-4 text-gray-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M8 12h.01M12 12h.01M16 12h.01M21 12c0 4.418-4.03 8-9 8a9.863 9.863 0 01-4.255-.949L3 20l1.395-3.72C3.512 15.042 3 13.574 3 12c0-4.418 4.03-8 9-8s9 3.582 9 8z"
                />
              </svg>
              <p>Select a conversation to view messages</p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
