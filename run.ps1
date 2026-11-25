# Universal Windows PowerShell Launcher
Write-Host "=========================================================================" -ForegroundColor Cyan
Write-Host "Starting Protein Mutation Analysis & Tertiary Structure Comparison Suite" -ForegroundColor Green
Write-Host "=========================================================================" -ForegroundColor Cyan

# Detect Python command
$pythonCmd = "python"
try {
    & python --version 2>$null
} catch {
    try {
        & py --version 2>$null
        $pythonCmd = "py"
    } catch {
        Write-Host "[ERROR] Python was not found in PATH!" -ForegroundColor Red
        Write-Host "Please install Python 3.10+ from python.org and check 'Add Python to PATH'." -ForegroundColor Yellow
        Pause
        Exit
    }
}

Write-Host "Using Python command: $pythonCmd" -ForegroundColor Gray

# Optional package install
Write-Host "`n[1/3] Ensuring required packages are installed..." -ForegroundColor Yellow
& $pythonCmd -m pip install -r requirements.txt

# Start FastAPI backend
Write-Host "[2/3] Starting FastAPI Backend on http://127.0.0.1:8000 ..." -ForegroundColor Yellow
Start-Process -FilePath $pythonCmd -ArgumentList "-m uvicorn backend.main:app --host 127.0.0.1 --port 8000" -NoNewWindow

Start-Sleep -Seconds 2

# Start Streamlit frontend
Write-Host "[3/3] Starting Streamlit Frontend on http://localhost:8501 ..." -ForegroundColor Yellow
Write-Host "`nOpening application in your default web browser...`n" -ForegroundColor Green
& $pythonCmd -m streamlit run frontend/app.py --server.port 8501
