# User Details Page - Admin Dashboard Feature

## Overview

A comprehensive user profile page for admins to view detailed information about individual users, including their full history, activity, and the ability to add notes and tags.

## Features Implemented

### 1. **User Profile Overview**
   - Complete user information (email, UID, created date)
   - Account status (active, banned, suspended)
   - Current subscription plan and status
   - Quick stats (total jobs, usage, monthly cost)
   - Visual tags with color coding

### 2. **Job History Tab**
   - Complete list of all jobs submitted by the user
   - Job details: sample name, job ID, status, current step
   - Creation and completion timestamps
   - Error messages for failed jobs
   - Color-coded status badges

### 3. **Payment History Tab**
   - Stripe webhook events related to the user
   - Payment intents, invoices, and charges
   - Amount and currency information
   - Processing status (processed/pending)
   - Chronological timeline

### 4. **Support Tickets Tab**
   - All support conversations initiated by the user
   - Subject, status, and message count
   - Last message timestamp
   - Direct link to conversation details
   - Status color coding (open, resolved, closed)

### 5. **Activity Timeline Tab**
   - Complete audit log of all actions related to the user
   - Admin actions (ban, subscription changes, usage updates)
   - User actions (job submissions, conversations)
   - Resource-level tracking
   - Metadata display for context

### 6. **Notes System**
   - Add internal notes about the user
   - Multi-line text support
   - Admin attribution (who added the note)
   - Timestamps (created and updated)
   - Historical note preservation

### 7. **Tags System**
   - Add colored tags to categorize users
   - 5 color options: blue, green, red, yellow, purple
   - Easy removal with confirmation
   - Visual display on user profile header
   - Searchable and filterable (future enhancement)

## Technical Implementation

### Backend Endpoints

#### Get User Details
```
GET /admin/users/{user_uid}/details
```
Returns comprehensive user information including:
- User profile and subscription
- Full job history
- Payment history (Stripe events)
- Support ticket history
- Activity timeline (audit logs)
- Notes and tags

#### Add Note
```
POST /admin/users/{user_uid}/notes?note_text=...
```
Adds an internal note to the user's profile.

#### Add Tag
```
POST /admin/users/{user_uid}/tags?tag_name=...&color=blue
```
Adds a colored tag to the user.

#### Remove Tag
```
DELETE /admin/users/{user_uid}/tags/{tag_id}
```
Removes a specific tag from the user.

### Database Schema

#### user_notes Table
```sql
CREATE TABLE user_notes (
    id INTEGER PRIMARY KEY,
    user_id VARCHAR(128) NOT NULL,
    admin_id VARCHAR(128) NOT NULL,
    note_text TEXT NOT NULL,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX idx_notes_user ON user_notes(user_id);
```

#### user_tags Table
```sql
CREATE TABLE user_tags (
    id INTEGER PRIMARY KEY,
    user_id VARCHAR(128) NOT NULL,
    tag_name VARCHAR(50) NOT NULL,
    color VARCHAR(20) DEFAULT 'blue',
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    created_by VARCHAR(128) NOT NULL
);
CREATE INDEX idx_tags_user ON user_tags(user_id);
CREATE INDEX idx_tags_name ON user_tags(tag_name);
```

### Frontend Components

#### UserDetailsPage Component
Location: `frontend/components/UserDetailsPage.tsx`

Features:
- Tab-based navigation
- Real-time data loading
- Interactive forms for notes and tags
- Responsive design
- Error handling and loading states

#### Navigation
Users can access the details page by clicking on any email address in the admin dashboard users table.

URL format: `/admin/user/{uid}`

## Usage Guide

### Accessing User Details

1. Navigate to the Admin Dashboard
2. Go to the "Users" tab
3. Click on any user's email address
4. You'll be redirected to their detailed profile page

### Adding Notes

1. Navigate to the "Notes" tab
2. Enter your note in the text area
3. Click "Add Note"
4. The note will appear in the list below with your admin ID and timestamp

### Adding Tags

1. On any tab, use the "Add Tag" form (also available in Profile tab)
2. Enter a tag name (e.g., "VIP", "High Priority", "Issue")
3. Select a color from the dropdown
4. Click "Add Tag"
5. Tags appear in the header for quick identification

### Removing Tags

1. Click the "×" button next to any tag in the header
2. Confirm the removal
3. The tag will be removed immediately

### Viewing Job History

1. Navigate to the "Jobs" tab
2. View all jobs with detailed information
3. Use the status colors to quickly identify job states:
   - Green: Completed
   - Blue: Running
   - Red: Failed
   - Gray: Pending

### Reviewing Payment History

1. Navigate to the "Payments" tab
2. View all Stripe events related to the user
3. Check amounts, currencies, and processing status
4. Useful for troubleshooting billing issues

### Checking Support History

1. Navigate to the "Support" tab
2. View all support tickets
3. See message counts and last activity
4. Click through to view full conversation (future enhancement)

### Reviewing Activity Timeline

1. Navigate to the "Activity" tab
2. View chronological list of all actions
3. Includes both admin actions and user actions
4. Metadata provides additional context

## Migration

To enable this feature on an existing installation:

```bash
# Run the migration script
cd /path/to/WholeExome
python migrations/add_user_notes_tags.py
```

This will create the `user_notes` and `user_tags` tables.

## API Integration

The feature uses the following API methods from `frontend/lib/api.ts`:

```typescript
// Get comprehensive user details
adminApi.getUserDetails(userUid: string): Promise<UserDetailsResponse>

// Add a note to user profile
adminApi.addUserNote(userUid: string, noteText: string): Promise<{success: boolean, note: UserNote}>

// Add a tag to user profile
adminApi.addUserTag(userUid: string, tagName: string, color: string): Promise<{success: boolean, tag: UserTag}>

// Remove a tag from user profile
adminApi.removeUserTag(userUid: string, tagId: number): Promise<{success: boolean, message: string}>
```

## Security

- All endpoints require admin authentication (`require_admin` middleware)
- Admin actions are logged in the audit log
- Notes include admin attribution
- Tags track who created them
- User IDs are validated before operations

## Future Enhancements

Potential improvements for the user details feature:

1. **Edit Notes** - Allow admins to edit existing notes
2. **Delete Notes** - Add ability to remove notes
3. **Tag Filtering** - Filter users by tags in main dashboard
4. **Quick Actions** - Add ban/suspend/upgrade buttons to profile header
5. **Export Profile** - Export complete user profile as PDF
6. **Communication** - Send direct email to user from profile
7. **Job Actions** - Retry/cancel jobs directly from user profile
8. **Usage Graphs** - Visual charts of usage over time
9. **Comparison** - Compare user metrics against averages
10. **Custom Fields** - Add custom metadata fields for users

## Files Modified/Created

### Backend
- `backend/database_extensions.py` - Added UserNote and UserTag models
- `backend/modules/admin/routes.py` - Added user details endpoints
- `backend/migrations/add_user_notes_tags.py` - Migration script

### Frontend
- `frontend/lib/api.ts` - Added API methods and interfaces
- `frontend/components/UserDetailsPage.tsx` - Main component
- `frontend/app/admin/user/[uid]/page.tsx` - Next.js page
- `frontend/components/AdminDashboard.tsx` - Added navigation links

## Testing Checklist

- [ ] User details page loads correctly
- [ ] All tabs display appropriate data
- [ ] Notes can be added successfully
- [ ] Tags can be added with different colors
- [ ] Tags can be removed
- [ ] Navigation from users table works
- [ ] Back button returns to admin dashboard
- [ ] Data refreshes after adding notes/tags
- [ ] Error handling works for invalid user IDs
- [ ] Loading states display properly
- [ ] Responsive design works on mobile
- [ ] Payment history filters correctly for user
- [ ] Activity timeline shows relevant events only

## Support

For issues or questions about this feature, please check:
1. Browser console for frontend errors
2. Backend logs for API errors
3. Database integrity (tables created, indexes present)
4. Admin authentication (user is in admin_users table)

## License

This feature is part of the WholeExome platform and follows the same license.
