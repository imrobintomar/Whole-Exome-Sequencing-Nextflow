#!/bin/bash

echo "=========================================="
echo "  WES Pipeline - Dependency Checker"
echo "=========================================="
echo ""

# Function to check if a command exists
check_command() {
    local cmd=$1
    local name=$2

    if command -v $cmd &> /dev/null; then
        local version=$($ Human: thew pipeline use conda create environment