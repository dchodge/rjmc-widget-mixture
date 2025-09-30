#!/bin/bash

# RJMCMC WebAssembly Build Script
# This script compiles the C++ code to WebAssembly using Emscripten

set -e

echo "🚀 Building RJMCMC WebAssembly Widget..."

# Check if Emscripten is installed
if ! command -v emcc &> /dev/null; then
    echo "❌ Emscripten not found. Please install Emscripten first:"
    echo "   https://emscripten.org/docs/getting_started/downloads.html"
    exit 1
fi

# Create build directory
mkdir -p build
cd build

# Configure with CMake
echo "📦 Configuring with CMake..."
emcmake cmake .. -DCMAKE_BUILD_TYPE=Release

# Build the project
echo "🔨 Building WebAssembly module..."
emmake make

# Copy the generated files to the root directory
echo "📋 Copying generated files..."
cp rjmc_wasm.js ../rjmc_wasm.js
cp rjmc_wasm.wasm ../rjmc_wasm.wasm

# Create a simple HTTP server script for testing
echo "🌐 Creating test server..."
cat > ../serve.py << 'EOF'
#!/usr/bin/env python3
import http.server
import socketserver
import os

PORT = 8000

class MyHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        super().end_headers()

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    with socketserver.TCPServer(("", PORT), MyHTTPRequestHandler) as httpd:
        print(f"🌐 Server running at http://localhost:{PORT}")
        print("📱 Open the widget in your browser!")
        httpd.serve_forever()
EOF

chmod +x ../serve.py

echo "✅ Build complete!"
echo ""
echo "🎯 To run the widget:"
echo "   cd /Users/davidhodgson/Dropbox/Mac\ \(3\)/Documents/research/software/mcmc/rjmc_widget"
echo "   python3 serve.py"
echo "   Then open http://localhost:8000 in your browser"
echo ""
echo "📁 Generated files:"
echo "   - rjmc_wasm.js (WebAssembly JavaScript glue code)"
echo "   - rjmc_wasm.wasm (WebAssembly binary)"
echo "   - serve.py (Test server script)"
