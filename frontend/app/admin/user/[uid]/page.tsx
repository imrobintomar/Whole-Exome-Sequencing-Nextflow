'use client';

import { use } from 'react';
import UserDetailsPage from '../../../../components/UserDetailsPage';

export default function UserDetailPage({ params }: { params: Promise<{ uid: string }> }) {
  const resolvedParams = use(params);
  return <UserDetailsPage userUid={resolvedParams.uid} />;
}
