#!/bin/bash

# Start R WebAssembly Server for RJMCMC Widget
# This script starts the R WebAssembly server that can run the rjmc package

echo "Starting R WebAssembly Server for RJMCMC Widget..."

# Check if R is installed
if ! command -v R &> /dev/null; then
    echo "Error: R is not installed. Please install R first."
    exit 1
fi

# Check if webr package is installed
if ! Rscript -e "library(webr)" 2>/dev/null; then
    echo "Installing webr package..."
    Rscript -e "install.packages('webr', repos='https://webr.r-wasm.org/latest')"
fi

# Check if rjmc package is available
if [ ! -d "/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/software/mcmc/rjmc" ]; then
    echo "Error: rjmc package directory not found."
    echo "Please ensure the rjmc package is available at the expected location."
    exit 1
fi

# Start the R WebAssembly server
echo "Starting R WebAssembly server..."
echo "The widget will be available at http://localhost:9090"
echo "Press Ctrl+C to stop the server"

# Start Python HTTP server in background
python3 -m http.server 9090 &
HTTP_PID=$!

# Start R WebAssembly server
Rscript r_wasm_server.R &
R_PID=$!

# Function to cleanup on exit
cleanup() {
    echo "Stopping servers..."
    kill $HTTP_PID 2>/dev/null
    kill $R_PID 2>/dev/null
    exit 0
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM

# Wait for processes
wait
