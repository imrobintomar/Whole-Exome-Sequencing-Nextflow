'use client';

import UserDetailsPage from '../../../../components/UserDetailsPage';

export default function UserDetailPage({ params }: { params: { uid: string } }) {
  return <UserDetailsPage userUid={params.uid} />;
}
