#!/bin/bash

# Setup proper R environment for macOS ARM64 with Homebrew

echo "Setting up R environment for oligo installation..."

# 1. Export the dynamic library fallback path for runtime
export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/opt/krb5/lib:$DYLD_FALLBACK_LIBRARY_PATH"

# 2. Create clean R 4.5 library directory
mkdir -p ~/Library/R/4.5/library

# 3. Clean up any partial installations
echo "Cleaning up partial installations..."
rm -rf /opt/homebrew/lib/R/4.5/site-library/00LOCK-*
rm -rf /opt/homebrew/lib/R/4.5/site-library/curl
rm -rf /opt/homebrew/lib/R/4.5/site-library/oligo
rm -rf /opt/homebrew/lib/R/4.5/site-library/oligoClasses

# 4. Install curl R package first
echo "Installing curl R package..."
R --vanilla -q -e 'install.packages("curl", type="source", repos="https://cloud.r-project.org", Ncpus = parallel::detectCores()-1)'

# 5. Check if curl installed successfully
echo "Checking curl installation..."
R --vanilla -q -e 'if ("curl" %in% installed.packages()[,"Package"]) cat("SUCCESS: curl is installed\n") else cat("ERROR: curl installation failed\n")'

# 6. If curl succeeded, install BiocManager and oligo
if R --vanilla -q -e 'q(status = as.numeric(!"curl" %in% installed.packages()[,"Package"]))' > /dev/null 2>&1; then
    echo "Installing BiocManager and oligo..."
    R --vanilla -q -e '
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager");
    BiocManager::install(version = "3.21");
    BiocManager::install(c("oligo","limma"), 
                         Ncpus = parallel::detectCores()-1, 
                         ask = FALSE)
    '
    
    # Check installation
    R --vanilla -q -e '
    pkgs <- c("oligo", "limma")
    installed <- pkgs %in% installed.packages()[,"Package"]
    cat("\nPackage installation status:\n")
    for (i in seq_along(pkgs)) {
        cat(sprintf("  %s: %s\n", pkgs[i], 
                    ifelse(installed[i], "INSTALLED", "FAILED")))
    }
    '
else
    echo "ERROR: curl installation failed. Cannot proceed with oligo installation."
    echo "Please check the error messages above."
fi

# 7. Add to .Rprofile for future sessions
echo "Updating .Rprofile..."
if ! grep -q "getRversion.*4.5" ~/.Rprofile 2>/dev/null; then
    echo 'if (getRversion() >= "4.5.0") .libPaths("~/Library/R/4.5/library")' >> ~/.Rprofile
    echo "Added R 4.5 library path to .Rprofile"
fi

echo "Setup complete!"