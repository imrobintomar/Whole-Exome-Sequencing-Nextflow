# MinIO Cost Analysis & Deployment Strategy
## 100% Free, Open-Source Object Storage for WES SaaS Platform

**Date:** 2026-01-14
**License:** AGPL v3 (Free for self-hosted deployments)
**Total MinIO License Cost:** $0 (Forever)

---

## 1. MinIO Licensing Clarification

### ✅ What's FREE (AGPL v3 - Open Source)

MinIO is **completely free** for self-hosted deployments with **no limitations**:

| Feature | Free (AGPL v3) | Commercial License |
|---------|----------------|-------------------|
| **Storage Capacity** | ✅ **Unlimited** (200 TB+) | ✅ Unlimited |
| **Number of Nodes** | ✅ **Unlimited** | ✅ Unlimited |
| **S3-Compatible API** | ✅ **Full Access** | ✅ Full Access |
| **Encryption (SSE-S3, KMS)** | ✅ **Included** | ✅ Included |
| **IAM, Bucket Policies** | ✅ **Included** | ✅ Included |
| **Lifecycle Management** | ✅ **Included** | ✅ Included |
| **Replication** | ✅ **Included** | ✅ Included |
| **Versioning** | ✅ **Included** | ✅ Included |
| **Object Locking (WORM)** | ✅ **Included** | ✅ Included |
| **Prometheus Metrics** | ✅ **Included** | ✅ Included |
| **Console UI** | ✅ **Included** | ✅ Included |
| **SDK Support** | ✅ **All Languages** | ✅ All Languages |
| **Nextflow Integration** | ✅ **S3-compatible** | ✅ S3-compatible |
| **Community Support** | ✅ Slack, GitHub, Docs | ❌ |
| **Enterprise Support** | ❌ | ✅ 24/7 Support |
| **Indemnification** | ❌ | ✅ Legal Protection |

### 📜 AGPL v3 Requirements

**You can use MinIO for FREE if you:**
- Self-host on your own infrastructure
- Don't modify MinIO source code (use as-is)
- If you modify MinIO code, you must open-source your changes

**Perfect for:**
- ✅ Genomics SaaS platforms (like yours)
- ✅ Healthcare/HIPAA applications
- ✅ Research institutions
- ✅ Startups
- ✅ Enterprise internal deployments

**When you DON'T need commercial license:**
- You're using MinIO as object storage backend
- You're not selling MinIO as a product
- You're not embedding MinIO in proprietary software
- You're using standard S3 APIs

---

## 2. Total Cost of Ownership (TCO)

### Your Current Setup: 200 TB Storage Available

**Hardware Costs (One-Time):**

| Item | Description | Cost | Notes |
|------|-------------|------|-------|
| **Storage** | 200 TB existing storage | **$0** | Already owned |
| **Server** | Existing server hardware | **$0** | Already owned |
| **MinIO Software** | AGPL v3 License | **$0** | Forever free |
| **Installation** | Self-service setup | **$0** | 1-2 days |
| **Configuration** | Bucket setup, policies | **$0** | 2-3 hours |
| **Total Upfront Cost** | | **$0** | |

### Operating Costs (Monthly)

| Item | Cost/Month | Notes |
|------|------------|-------|
| **MinIO License** | **$0** | AGPL v3 - Free |
| **Electricity** | $50-150 | Depends on server power draw |
| **Internet** | $0 | Included in existing connection |
| **Maintenance** | $0 | Self-managed |
| **Monitoring** | $0 | Prometheus + Grafana (free) |
| **Backups** | $0 | On-premises replication |
| **Total Operating Cost** | **$50-150** | Just electricity |

### Cost Comparison: MinIO vs Cloud Storage (Annual)

**Assumptions:**
- 100 WES jobs/month = 1,200 jobs/year
- 50 GB per job (FASTQ + BAM + VCF + results)
- 60 TB total storage (growing)
- 120 TB annual egress (users downloading results)

| Provider | Storage Cost | Egress Cost | Annual Total |
|----------|--------------|-------------|--------------|
| **AWS S3 Standard** | $1,380/year | $11,040/year | **$12,420/year** |
| **AWS S3 Glacier** | $240/year | $11,040/year | **$11,280/year** |
| **Google Cloud Storage** | $1,260/year | $9,600/year | **$10,860/year** |
| **Azure Blob** | $1,200/year | $10,200/year | **$11,400/year** |
| **MinIO (Self-Hosted)** | $0 | $0 | **$600-1,800/year** ⚡ |

**Annual Savings with MinIO:** **$9,000 - $11,820**

**5-Year Savings:** **$45,000 - $59,100**

**ROI:** ∞ (infinite) - No licensing costs, just electricity

---

## 3. Deployment Options for Your 200 TB Storage

### Option A: Single-Node Deployment (Simplest)

**Best for:** Development, small-scale production (<50 users)

```
┌────────────────────────────────────────────┐
│         Single MinIO Server                 │
│                                            │
│  Storage: 200 TB (your existing drives)    │
│  RAM: 32-64 GB (recommended)              │
│  CPU: 8-16 cores                          │
│  Network: 10 Gbps (optional)              │
│                                            │
│  Redundancy: RAID 6 (hardware level)      │
│  Backup: Nightly snapshots to NAS         │
└────────────────────────────────────────────┘
```

**Pros:**
- ✅ Simple setup (30 minutes)
- ✅ Low complexity
- ✅ No distributed coordination needed

**Cons:**
- ❌ Single point of failure
- ❌ No horizontal scaling

**Setup:**
```bash
# Install MinIO (single binary)
wget https://dl.min.io/server/minio/release/linux-amd64/minio
chmod +x minio
sudo mv minio /usr/local/bin/

# Start MinIO pointing to your 200 TB storage
export MINIO_ROOT_USER=admin
export MINIO_ROOT_PASSWORD=YourSecurePassword123!

minio server /mnt/200tb-storage --console-address ":9001"
```

**That's it! MinIO is running on your 200 TB storage.**

---

### Option B: Distributed Multi-Node (High Availability)

**Best for:** Production with 100+ concurrent users

If you have **multiple servers or drives**, split your 200 TB:

```
┌─────────────────────────────────────────────────────────────┐
│                  Distributed MinIO Cluster                   │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Node 1: 50 TB    Node 2: 50 TB    Node 3: 50 TB    Node 4: 50 TB │
│    [Drive 1]        [Drive 1]        [Drive 1]        [Drive 1]   │
│    [Drive 2]        [Drive 2]        [Drive 2]        [Drive 2]   │
│                                                              │
│  Erasure Coding: EC:4 (N/2 redundancy)                      │
│  Can lose 2 nodes and still operate                         │
└─────────────────────────────────────────────────────────────┘
```

**Setup:**
```bash
# On each node, run MinIO pointing to all nodes
export MINIO_ROOT_USER=admin
export MINIO_ROOT_PASSWORD=YourSecurePassword123!

minio server \
  http://node{1...4}/mnt/disk{1...2}/minio \
  --console-address ":9001"
```

**Pros:**
- ✅ High availability (node failure tolerance)
- ✅ Horizontal scaling
- ✅ Better performance (parallel I/O)

**Cons:**
- ⚠️ More complex setup
- ⚠️ Requires multiple servers

---

## 4. Recommended Deployment for Your Use Case

### 🎯 **Recommended: Single-Node with RAID 6 + Nightly Backups**

Given your 200 TB storage, start with:

```
┌───────────────────────────────────────────────────────────┐
│                     Production Setup                       │
├───────────────────────────────────────────────────────────┤
│                                                            │
│  MinIO Server (Single Node)                               │
│  ├── Storage: 200 TB RAID 6 (existing)                    │
│  ├── RAM: 64 GB (for caching metadata)                    │
│  ├── CPU: 16 cores (for parallel uploads)                 │
│  ├── Network: 10 Gbps (optional, improves throughput)     │
│  └── OS: Ubuntu 22.04 LTS                                 │
│                                                            │
│  Backup Strategy:                                         │
│  ├── Daily: MinIO → NAS (rsync/rclone)                   │
│  ├── Weekly: Full snapshot                                │
│  └── Monthly: Offsite backup                              │
│                                                            │
│  Monitoring:                                               │
│  ├── Prometheus: Metrics collection                       │
│  ├── Grafana: Dashboards                                  │
│  └── Alertmanager: Disk usage, errors                     │
│                                                            │
│  Total Cost: $0 (software) + $100/month (electricity)     │
└───────────────────────────────────────────────────────────┘
```

**Capacity Planning:**
- 200 TB raw storage
- ~180 TB usable (after RAID 6 overhead)
- ~3,600 WES jobs (at 50 GB per job)
- With lifecycle policies: **>10,000 jobs** (archival/deletion)

---

## 5. Installation Guide (Production-Ready)

### Step 1: Install MinIO Server (5 minutes)

```bash
#!/bin/bash
# install-minio.sh - Production MinIO installation

# Download MinIO binary
wget https://dl.min.io/server/minio/release/linux-amd64/minio
chmod +x minio
sudo mv minio /usr/local/bin/

# Download MinIO client (mc)
wget https://dl.min.io/client/mc/release/linux-amd64/mc
chmod +x mc
sudo mv mc /usr/local/bin/

# Create minio user
sudo useradd -r -s /sbin/nologin minio

# Create data directory (point to your 200 TB mount)
DATA_DIR="/mnt/200tb-storage/minio"
sudo mkdir -p ${DATA_DIR}
sudo chown -R minio:minio ${DATA_DIR}

echo "✓ MinIO binaries installed"
```

### Step 2: Configure SystemD Service (Production)

```bash
# Create systemd service
sudo tee /etc/systemd/system/minio.service << 'EOF'
[Unit]
Description=MinIO Object Storage
Documentation=https://docs.min.io
After=network-online.target
Wants=network-online.target

[Service]
Type=notify
WorkingDirectory=/usr/local/
User=minio
Group=minio

# MinIO credentials (CHANGE THESE!)
Environment="MINIO_ROOT_USER=admin"
Environment="MINIO_ROOT_PASSWORD=ChangeThisToAStrongPassword123!"

# Storage location (your 200 TB)
Environment="MINIO_VOLUMES=/mnt/200tb-storage/minio"

# Server options
Environment="MINIO_OPTS=--console-address :9001"

# Performance tuning
Environment="MINIO_API_REQUESTS_MAX=10000"
Environment="MINIO_API_REQUESTS_DEADLINE=10s"

# Enable encryption at rest (optional)
# Environment="MINIO_KMS_SECRET_KEY=my-minio-key:REPLACE-WITH-32-CHAR-KEY"

# Start MinIO
ExecStartPre=/bin/bash -c "if [ -z \"${MINIO_VOLUMES}\" ]; then echo \"MINIO_VOLUMES not set\"; exit 1; fi"
ExecStart=/usr/local/bin/minio server $MINIO_OPTS $MINIO_VOLUMES

# Restart policy
Restart=always
RestartSec=5
LimitNOFILE=65536
TasksMax=infinity
TimeoutStopSec=infinity
SendSIGKILL=no

[Install]
WantedBy=multi-user.target
EOF

# Reload and start
sudo systemctl daemon-reload
sudo systemctl enable minio
sudo systemctl start minio

# Check status
sudo systemctl status minio
```

### Step 3: Verify Installation

```bash
# Check MinIO is running
curl http://localhost:9000/minio/health/live

# Expected output: HTTP 200 OK

# Access MinIO Console
echo "MinIO Console: http://your-server-ip:9001"
echo "Username: admin"
echo "Password: ChangeThisToAStrongPassword123!"
```

### Step 4: Configure MinIO Client

```bash
# Configure mc client
mc alias set myminio http://localhost:9000 admin ChangeThisToAStrongPassword123!

# Test connection
mc admin info myminio

# Expected output:
# ●  localhost:9000
#    Uptime: X hours
#    Version: RELEASE.XXXX-XX-XXTXX-XX-XXZ
#    Storage: 200 TB Free
```

### Step 5: Create Buckets

```bash
#!/bin/bash
# setup-buckets.sh - Initialize all WES buckets

ALIAS="myminio"

echo "Creating buckets..."

# Create buckets
mc mb ${ALIAS}/wes-raw-data
mc mb ${ALIAS}/wes-intermediate
mc mb ${ALIAS}/wes-results
mc mb ${ALIAS}/wes-archives
mc mb ${ALIAS}/wes-reference
mc mb ${ALIAS}/wes-logs

echo "✓ Buckets created"

# Enable versioning (HIPAA compliance)
echo "Enabling versioning..."
mc version enable ${ALIAS}/wes-raw-data
mc version enable ${ALIAS}/wes-results
mc version enable ${ALIAS}/wes-archives

# Enable encryption
echo "Enabling encryption..."
mc encrypt set sse-s3 ${ALIAS}/wes-raw-data
mc encrypt set sse-s3 ${ALIAS}/wes-intermediate
mc encrypt set sse-s3 ${ALIAS}/wes-results
mc encrypt set sse-s3 ${ALIAS}/wes-archives

# Set lifecycle policies
echo "Configuring lifecycle policies..."

# Raw data: Delete after 7 days
cat > /tmp/raw-lifecycle.json << 'POLICY'
{
  "Rules": [
    {
      "ID": "DeleteAfter7Days",
      "Status": "Enabled",
      "Expiration": {
        "Days": 7
      }
    }
  ]
}
POLICY
mc ilm import ${ALIAS}/wes-raw-data < /tmp/raw-lifecycle.json

# Intermediate: Delete after 30 days
cat > /tmp/intermediate-lifecycle.json << 'POLICY'
{
  "Rules": [
    {
      "ID": "DeleteAfter30Days",
      "Status": "Enabled",
      "Expiration": {
        "Days": 30
      }
    }
  ]
}
POLICY
mc ilm import ${ALIAS}/wes-intermediate < /tmp/intermediate-lifecycle.json

# Results: Keep 90 days, then archive
cat > /tmp/results-lifecycle.json << 'POLICY'
{
  "Rules": [
    {
      "ID": "DeleteAfter90Days",
      "Status": "Enabled",
      "Expiration": {
        "Days": 90
      }
    }
  ]
}
POLICY
mc ilm import ${ALIAS}/wes-results < /tmp/results-lifecycle.json

echo "✓ Lifecycle policies configured"
echo "✓ MinIO setup complete!"

# Display bucket info
mc ls ${ALIAS}
```

---

## 6. Upload Reference Genome to MinIO

Your hg38 reference and ANNOVAR databases should be stored in MinIO:

```bash
#!/bin/bash
# upload-reference-data.sh - Upload genomic reference files

ALIAS="myminio"
BUCKET="wes-reference"

echo "Uploading reference genome..."

# Upload hg38 reference
mc cp /media/drprabudh/m1/hg38/hg38.fa ${ALIAS}/${BUCKET}/hg38/
mc cp /media/drprabudh/m1/hg38/hg38.fa.fai ${ALIAS}/${BUCKET}/hg38/
mc cp /media/drprabudh/m1/hg38/hg38.dict ${ALIAS}/${BUCKET}/hg38/

# Upload BWA indices
mc cp /media/drprabudh/m1/hg38/hg38.fa.amb ${ALIAS}/${BUCKET}/hg38/
mc cp /media/drprabudh/m1/hg38/hg38.fa.ann ${ALIAS}/${BUCKET}/hg38/
mc cp /media/drprabudh/m1/hg38/hg38.fa.bwt ${ALIAS}/${BUCKET}/hg38/
mc cp /media/drprabudh/m1/hg38/hg38.fa.pac ${ALIAS}/${BUCKET}/hg38/
mc cp /media/drprabudh/m1/hg38/hg38.fa.sa ${ALIAS}/${BUCKET}/hg38/

echo "Uploading ANNOVAR databases..."

# Upload ANNOVAR databases (recursive)
mc cp --recursive /media/drprabudh/m1/annovar/hg38_humandb/ \
  ${ALIAS}/${BUCKET}/annovar/hg38_humandb/

echo "Uploading known sites VCFs..."

# Upload known sites
mc cp /media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz \
  ${ALIAS}/${BUCKET}/known_sites/
mc cp /media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
  ${ALIAS}/${BUCKET}/known_sites/

mc cp /media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  ${ALIAS}/${BUCKET}/known_sites/
mc cp /media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
  ${ALIAS}/${BUCKET}/known_sites/

echo "✓ Reference data uploaded to MinIO"

# Verify
mc ls ${ALIAS}/${BUCKET}/hg38/
mc ls ${ALIAS}/${BUCKET}/annovar/hg38_humandb/ | head
mc ls ${ALIAS}/${BUCKET}/known_sites/
```

**Storage Used:**
- hg38 reference + indices: ~3.5 GB
- ANNOVAR databases: ~50 GB
- Known sites VCFs: ~2 GB
- **Total:** ~55 GB (out of your 200 TB)

---

## 7. Performance Benchmarking

### Expected Performance on Your Hardware

**Single-Node MinIO (200 TB RAID 6):**

| Operation | Speed | Notes |
|-----------|-------|-------|
| Upload (large file) | 100-500 MB/s | Limited by disk write speed |
| Download (large file) | 200-800 MB/s | Limited by disk read speed |
| Small file upload | 1,000-5,000 ops/s | Depends on IOPS |
| Parallel uploads (10 clients) | 800 MB/s - 2 GB/s | Network becomes bottleneck |
| Metadata operations | 10,000-50,000 ops/s | RAM cached |

**10 GB FASTQ Upload Test:**
```bash
# Create 10 GB test file
dd if=/dev/urandom of=/tmp/test-10gb.fastq.gz bs=1M count=10240

# Upload with mc
time mc cp /tmp/test-10gb.fastq.gz myminio/wes-raw-data/test/

# Expected: 20-100 seconds (100-500 MB/s)
```

### Optimization Tips

**1. Increase File Descriptor Limits:**
```bash
# /etc/security/limits.conf
minio soft nofile 65536
minio hard nofile 65536

# /etc/systemd/system/minio.service
LimitNOFILE=65536
```

**2. Tune Kernel Parameters:**
```bash
# /etc/sysctl.conf
net.core.rmem_max = 134217728
net.core.wmem_max = 134217728
net.ipv4.tcp_rmem = 4096 87380 67108864
net.ipv4.tcp_wmem = 4096 65536 67108864
net.core.netdev_max_backlog = 5000

# Apply
sudo sysctl -p
```

**3. Use SSD for Metadata (Optional):**
If you have SSDs, use them for MinIO metadata:
```bash
# Move .minio.sys to SSD
ln -s /mnt/ssd/.minio.sys /mnt/200tb-storage/minio/.minio.sys
```

---

## 8. Monitoring Setup (Free Tools)

### Prometheus + Grafana Dashboard

**Install Prometheus:**
```bash
# Download Prometheus
wget https://github.com/prometheus/prometheus/releases/download/v2.45.0/prometheus-2.45.0.linux-amd64.tar.gz
tar xvf prometheus-2.45.0.linux-amd64.tar.gz
cd prometheus-2.45.0.linux-amd64

# Configure MinIO as target
cat > prometheus.yml << 'EOF'
global:
  scrape_interval: 15s

scrape_configs:
  - job_name: 'minio'
    bearer_token: YOUR_MINIO_PROMETHEUS_TOKEN
    metrics_path: /minio/v2/metrics/cluster
    scheme: http
    static_configs:
      - targets: ['localhost:9000']
EOF

# Start Prometheus
./prometheus --config.file=prometheus.yml
```

**Generate MinIO Prometheus Token:**
```bash
mc admin prometheus generate myminio
# Copy the bearer_token and add to prometheus.yml
```

**Install Grafana:**
```bash
sudo apt-get install -y grafana
sudo systemctl enable grafana-server
sudo systemctl start grafana-server

# Access: http://localhost:3000
# Default: admin/admin
```

**Import MinIO Dashboard:**
1. Open Grafana → Dashboards → Import
2. Enter Dashboard ID: `13502` (Official MinIO Dashboard)
3. Select Prometheus data source
4. Click Import

**Metrics to Monitor:**
- Storage usage per bucket
- Upload/download throughput
- API request rate
- Error rate
- Disk I/O utilization

---

## 9. Backup Strategy (Free)

### Option A: Bucket Replication (MinIO-to-MinIO)

If you have a second MinIO instance (DR site):

```bash
# Add remote site
mc admin bucket remote add myminio/wes-results \
  https://dr-site.example.com/wes-results-replica \
  --service replication \
  --access-key DR_ACCESS_KEY \
  --secret-key DR_SECRET_KEY

# Enable replication
mc replicate add myminio/wes-results \
  --remote-bucket wes-results-replica \
  --replicate "delete,delete-marker"

# Verify
mc replicate status myminio/wes-results
```

### Option B: Rclone to External Storage

```bash
# Install rclone (free)
sudo apt install rclone

# Configure rclone for MinIO
rclone config

# Daily backup script
cat > /opt/scripts/minio-backup.sh << 'EOF'
#!/bin/bash
DATE=$(date +%Y%m%d)

# Sync critical buckets to NAS/external storage
rclone sync myminio:wes-results /mnt/nas/minio-backup/wes-results/ \
  --progress \
  --checksum

rclone sync myminio:wes-archives /mnt/nas/minio-backup/wes-archives/ \
  --progress \
  --checksum

echo "✓ Backup completed: ${DATE}"
EOF

chmod +x /opt/scripts/minio-backup.sh

# Add to cron (daily at 2 AM)
(crontab -l 2>/dev/null; echo "0 2 * * * /opt/scripts/minio-backup.sh >> /var/log/minio-backup.log 2>&1") | crontab -
```

---

## 10. Security Hardening (Free)

### Enable TLS/HTTPS

```bash
# Generate Let's Encrypt certificate (free)
sudo apt install certbot
sudo certbot certonly --standalone -d minio.yourdomain.com

# Copy certs to MinIO directory
sudo mkdir -p /etc/minio/certs
sudo cp /etc/letsencrypt/live/minio.yourdomain.com/fullchain.pem /etc/minio/certs/public.crt
sudo cp /etc/letsencrypt/live/minio.yourdomain.com/privkey.pem /etc/minio/certs/private.key
sudo chown -R minio:minio /etc/minio/certs

# Update systemd service
sudo systemctl edit minio
# Add: Environment="MINIO_OPTS=--certs-dir /etc/minio/certs --console-address :9001"

# Restart MinIO
sudo systemctl restart minio

# Verify HTTPS
curl https://minio.yourdomain.com:9000/minio/health/live
```

### Firewall Rules

```bash
# Allow only necessary ports
sudo ufw allow 22/tcp    # SSH
sudo ufw allow 9000/tcp  # MinIO API
sudo ufw allow 9001/tcp  # MinIO Console
sudo ufw enable

# Restrict API access to application servers only
sudo ufw allow from YOUR_APP_SERVER_IP to any port 9000
```

---

## 11. Final Cost Summary

### Total Cost Breakdown

| Component | Cost | Frequency | Notes |
|-----------|------|-----------|-------|
| **Software Costs** | | | |
| MinIO AGPL v3 License | **$0** | Forever | Open source |
| Prometheus | $0 | Forever | Open source |
| Grafana | $0 | Forever | Open source |
| Rclone | $0 | Forever | Open source |
| **Subtotal Software** | **$0** | | |
| | | | |
| **Hardware Costs** | | | |
| 200 TB Storage | $0 | N/A | Already owned |
| Server | $0 | N/A | Already owned |
| **Subtotal Hardware** | **$0** | | |
| | | | |
| **Operating Costs** | | | |
| Electricity (500W server) | $50-150 | Monthly | ~$100/month avg |
| Internet/Bandwidth | $0 | Monthly | Existing connection |
| Maintenance | $0 | Monthly | Self-managed |
| **Subtotal Operating** | **$100** | **Monthly** | |
| | | | |
| **TOTAL FIRST YEAR** | **$1,200** | | Just electricity |
| **TOTAL ANNUAL (ongoing)** | **$1,200** | | Just electricity |

### Comparison to AWS S3 (5-Year TCO)

| Solution | Year 1 | Year 2-5 | 5-Year Total |
|----------|--------|----------|--------------|
| **MinIO (Self-Hosted)** | $1,200 | $4,800 | **$6,000** |
| **AWS S3** | $12,420 | $49,680 | **$62,100** |
| **Savings with MinIO** | | | **$56,100** |

**ROI:** 10x return on investment over 5 years

---

## 12. Quick Start Checklist

- [ ] Download MinIO binary (1 minute)
- [ ] Create systemd service (5 minutes)
- [ ] Start MinIO server (1 minute)
- [ ] Access console at http://server:9001 (1 minute)
- [ ] Create buckets (2 minutes)
- [ ] Configure lifecycle policies (3 minutes)
- [ ] Upload reference genome (30 minutes)
- [ ] Test upload/download (5 minutes)
- [ ] Configure Prometheus monitoring (15 minutes)
- [ ] Set up daily backups (10 minutes)
- [ ] Enable HTTPS (15 minutes)

**Total Setup Time:** ~2 hours

---

## 13. Support Resources (All Free)

- 📖 **Official Docs:** https://docs.min.io
- 💬 **Slack Community:** https://slack.min.io
- 🐛 **GitHub Issues:** https://github.com/minio/minio/issues
- 📹 **YouTube Tutorials:** https://www.youtube.com/c/MinioInc
- 📚 **Cookbook:** https://github.com/minio/cookbook

---

## Conclusion

**MinIO is 100% free for self-hosted deployments** under AGPL v3 license.

With your existing 200 TB storage, you can:
- ✅ Deploy enterprise-grade object storage at **$0 licensing cost**
- ✅ Save **$56,000+ over 5 years** vs AWS S3
- ✅ Store **10,000+ WES jobs** with lifecycle management
- ✅ Achieve **clinical-grade HIPAA compliance** with encryption
- ✅ Scale to petabytes if needed (no license restrictions)

**Next Step:** Run the installation script and you'll have MinIO running in 2 hours!

---

**License:** This document is free to use and modify.
**MinIO License:** AGPL v3 - https://github.com/minio/minio/blob/master/LICENSE
