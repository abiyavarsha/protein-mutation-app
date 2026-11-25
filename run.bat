@echo off
setlocal enabledelayedexpansion

echo =========================================================================
echo Starting Protein Mutation Analysis & Tertiary Structure Comparison Suite
echo =========================================================================

:: Detect Python
set PYTHON_CMD=python
python --version >nul 2>&1
if %errorlevel% neq 0 (
    py --version >nul 2>&1
    if %errorlevel% equ 0 (
        set PYTHON_CMD=py
    ) else (
        echo [ERROR] Python was not found in PATH!
        echo Please install Python 3.10+ and make sure to check "Add Python to PATH".
        pause
        exit /b 1
    )
)

echo Detected Python: %PYTHON_CMD%

:: Ensure dependencies are installed
echo [1/3] Checking requirements...
%PYTHON_CMD% -m pip install -r requirements.txt >nul 2>&1

:: Launch servers
echo [2/3] Starting FastAPI Backend on http://127.0.0.1:8000 ...
start /B "" %PYTHON_CMD% -m uvicorn backend.main:app --host 127.0.0.1 --port 8000

:: Small delay
timeout /t 2 /nobreak >nul

echo [3/3] Starting Streamlit Web App on http://localhost:8501 ...
echo.
echo Application opening in your default web browser...
echo.
%PYTHON_CMD% -m streamlit run frontend/app.py --server.port 8501

pause
