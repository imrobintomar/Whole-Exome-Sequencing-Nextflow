# Backend Directory

## ⚠️ IMPORTANT: How to Run the Pipeline

### From this directory, use:
```bash
./run.sh -resume
```

### DON'T use:
```bash
nextflow run ../main.nf  # ❌ Missing JVM settings, will cause StackOverflowError
```

## Why?

The `run.sh` script sets the proper Java memory settings to prevent StackOverflowError:
```bash
export NXF_OPTS="-Xms2g -Xmx8g -Xss4m"
```

## See Also
- [QUICK_FIX.md](../QUICK_FIX.md) - Detailed fix guide
- [TROUBLESHOOTING.md](../TROUBLESHOOTING.md) - Full troubleshooting guide
