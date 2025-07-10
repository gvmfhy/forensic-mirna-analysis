#!/bin/bash

# Fix R library paths on macOS ARM64

echo "Setting up R environment variables for macOS ARM64..."

# Export paths for current session
export LDFLAGS="-L/opt/homebrew/opt/krb5/lib"
export CPPFLAGS="-I/opt/homebrew/opt/krb5/include"
export PKG_CONFIG_PATH="/opt/homebrew/opt/krb5/lib/pkgconfig"

# Create R environment file
cat > ~/.Renviron << EOF
LDFLAGS=-L/opt/homebrew/opt/krb5/lib
CPPFLAGS=-I/opt/homebrew/opt/krb5/include
PKG_CONFIG_PATH=/opt/homebrew/opt/krb5/lib/pkgconfig
EOF

echo "Environment variables set. Now trying to install curl..."

# Try installing curl package with correct paths
R --slave << EOF
Sys.setenv(LDFLAGS = "-L/opt/homebrew/opt/krb5/lib")
Sys.setenv(CPPFLAGS = "-I/opt/homebrew/opt/krb5/include")
install.packages("curl", repos = "https://cloud.r-project.org", type = "source")
EOF