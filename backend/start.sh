#!/bin/bash
# WES Backend Startup Script

export PYTHONPATH="/home/drprabudh/.local/lib/python3.10/site-packages:$PYTHONPATH"
python3 -m uvicorn main:app --reload
