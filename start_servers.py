import subprocess
import sys
import time
import os

python_exe = sys.executable
app_dir = os.path.dirname(os.path.abspath(__file__))

print(f"Using Python: {python_exe}")
print(f"Working Dir: {app_dir}")

# 1. Start FastAPI Backend
print("Starting FastAPI Backend on http://127.0.0.1:8000 ...")
backend_proc = subprocess.Popen(
    [
        python_exe, "-m", "uvicorn", "backend.main:app",
        "--host", "127.0.0.1",
        "--port", "8000"
    ],
    cwd=app_dir
)

# Wait a brief moment for backend to initialize
time.sleep(2)

# 2. Start Streamlit Frontend
print("Starting Streamlit Frontend on http://localhost:8501 ...")
frontend_proc = subprocess.Popen(
    [
        python_exe, "-m", "streamlit", "run", "frontend/app.py",
        "--server.port", "8501",
        "--server.headless", "true",
        "--browser.serverAddress", "localhost",
        "--server.enableCORS", "false",
        "--server.enableXsrfProtection", "false"
    ],
    cwd=app_dir
)

print("\n" + "="*60)
print("SUCCESS: Both services are now running live!")
print("Frontend URL: http://localhost:8501")
print("Backend API URL: http://127.0.0.1:8000")
print("Backend Swagger Docs: http://127.0.0.1:8000/docs")
print("="*60 + "\n")

try:
    while True:
        time.sleep(1)
        if backend_proc.poll() is not None:
            print("Backend exited unexpectedly!")
            break
        if frontend_proc.poll() is not None:
            print("Frontend exited unexpectedly!")
            break
except KeyboardInterrupt:
    print("Shutting down servers...")
    backend_proc.terminate()
    frontend_proc.terminate()
