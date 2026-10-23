#!/usr/bin/env python3
"""
Local dev server for m2g HTML viewer.
Serves static files + handles /api/run for mmseqs commands.
Usage: python3 server.py [port]  (default: 8080)
Run from the directory where mmseqs was executed (e.g., examples/mmseqs2go/).
"""
import http.server
import subprocess
import json
import sys
import os

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8080

class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        # Normalize double-leading-slash (e.g. //storage2/...) to absolute path
        path = self.path.split('?')[0]
        abs_path = '/' + path.lstrip('/')
        if os.path.isfile(abs_path):
            try:
                with open(abs_path, 'rb') as f:
                    data = f.read()
                self.send_response(200)
                self.send_header('Content-Type', 'application/octet-stream')
                self.send_header('Content-Length', str(len(data)))
                self.end_headers()
                self.wfile.write(data)
            except Exception as e:
                self.send_response(500)
                self.end_headers()
        else:
            super().do_GET()

    def do_POST(self):
        if self.path == '/api/run':
            length = int(self.headers.get('Content-Length', 0))
            body = json.loads(self.rfile.read(length))
            cmd = body.get('cmd', '')
            cwd = body.get('cwd') or None
            try:
                result = subprocess.run(
                    cmd, shell=True, capture_output=True, text=True, timeout=30,
                    cwd=cwd
                )
                output = result.stdout
                if result.stderr:
                    output += '\n[stderr]\n' + result.stderr
            except subprocess.TimeoutExpired:
                output = '[Error] Command timed out (30s)'
            except Exception as e:
                output = f'[Error] {e}'
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            self.wfile.write(json.dumps({'output': output}).encode())
        else:
            self.send_response(405)
            self.end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.end_headers()

    def log_message(self, fmt, *args):
        print(f'  {self.address_string()} - {fmt % args}')

if __name__ == '__main__':
    os.chdir(os.getcwd())
    print(f'Serving at http://localhost:{PORT}')
    print(f'Run from: {os.getcwd()}')
    http.server.test(HandlerClass=Handler, port=PORT, bind='127.0.0.1')
