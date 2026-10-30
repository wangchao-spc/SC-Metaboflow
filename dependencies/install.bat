@echo off
setlocal
chcp 65001 >nul
title SC-MetaboFlow Setup

echo ================================================
echo   SC-MetaboFlow (English) - One-Click Setup
echo ================================================
echo.

REM --- 1. Check Python ---
where python >nul 2>nul
if errorlevel 1 (
    echo [ERROR] Python is not installed or not in PATH.
    echo.
    echo   1. Download Python from https://www.python.org/downloads/
    echo   2. Run the installer and CHECK "Add python.exe to PATH"
    echo   3. Run this setup again
    echo.
    pause
    exit /b 1
)

python -c "import sys; sys.exit(0 if sys.version_info >= (3,12) else 1)"
if errorlevel 1 (
    echo [ERROR] Python 3.12 or newer is required. Found:
    python --version
    echo.
    pause
    exit /b 1
)
echo [OK] Python found:
python --version

REM --- 2. Create virtual environment ---
".venv\Scripts\python.exe" -c "import sys" >nul 2>nul
if errorlevel 1 (
    echo.
    echo Recreating virtual environment .venv ...
    rmdir /s /q .venv 2>nul
    python -m venv .venv
    if errorlevel 1 (
        echo [ERROR] Failed to create .venv
        pause
        exit /b 1
    )
) else (
    echo [OK] Virtual environment .venv already exists
)

REM --- 3. Install dependencies ---
echo.
echo Installing dependencies (first time takes a few minutes)...
".venv\Scripts\python.exe" -m pip install --upgrade pip >nul 2>nul
".venv\Scripts\python.exe" -m pip install -r requirements_full.txt
if errorlevel 1 (
    echo [ERROR] Dependency installation failed. Check the internet connection.
    pause
    exit /b 1
)

REM --- 4. Create desktop shortcut ---
echo.
echo Creating desktop shortcut...
powershell -NoProfile -Command "$ws=New-Object -ComObject WScript.Shell; $lnk=$ws.CreateShortcut([Environment]::GetFolderPath('Desktop')+'\SC-MetaboFlow.lnk'); $lnk.TargetPath='%~dp0.venv\Scripts\pythonw.exe'; $lnk.Arguments='launcher.pyw'; $lnk.WorkingDirectory='%~dp0'; $lnk.Save()"

echo.
echo ================================================
echo   Setup complete. Launching SC-MetaboFlow...
echo ================================================
start "" "%~dp0.venv\Scripts\pythonw.exe" "%~dp0launcher.pyw"
endlocal
