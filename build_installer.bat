@echo off
setlocal

:: ============================================================
:: CoolSolve Windows Build + Installer Script
::
:: Requirements:
::   - Visual Studio 2022 with C++ Desktop Development workload
::   - Node.js (https://nodejs.org) for the React GUI frontend
::   - Python (e.g. Anaconda) for CoolProp build scripts
::   - NSIS (https://nsis.sourceforge.net) for the installer
::
:: Update the paths in Section 1 if your installation differs.
:: ============================================================

:: --- Section 1: Configuration ---
set "ANACONDA_DIR=E:\Anaconda"
set "CONDA_ENV=fast_numpy_env"
set "CMAKE_EXE=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
set "MAKENSIS_EXE=C:\Program Files (x86)\NSIS\makensis.exe"
if not exist "%MAKENSIS_EXE%" set "MAKENSIS_EXE=C:\Program Files\NSIS\makensis.exe"

:: Root of the project (directory containing this bat file)
set "PROJECT_DIR=%~dp0"
set "PROJECT_DIR=%PROJECT_DIR:~0,-1%"

:: --- Section 2: Add Anaconda + scripts to PATH ---
echo [1/5] Setting up environment...
set "PATH=%ANACONDA_DIR%;%ANACONDA_DIR%\Scripts;%ANACONDA_DIR%\Library\bin;%PATH%"
call "%ANACONDA_DIR%\Scripts\activate.bat" %CONDA_ENV%

:: --- Section 3: Build the React GUI frontend ---
echo [2/5] Building React GUI frontend (npm run build)...
where node >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo WARNING: Node.js not found in PATH. Skipping GUI frontend build.
    echo          The server will launch in dev mode ^(GUI served from filesystem^).
    echo          To embed the GUI, install Node.js ^(https://nodejs.org^) and re-run.
    goto :skip_npm
)
cd "%PROJECT_DIR%\gui"
npm ci
if %ERRORLEVEL% neq 0 (
    echo Error running npm ci.
    pause
    exit /b %ERRORLEVEL%
)
npm run build
if %ERRORLEVEL% neq 0 (
    echo Error building the GUI frontend.
    pause
    exit /b %ERRORLEVEL%
)
:skip_npm

:: --- Section 4: Configure CMake ---
echo [3/5] Configuring CMake...
if not exist "%PROJECT_DIR%\build" mkdir "%PROJECT_DIR%\build"
cd "%PROJECT_DIR%\build"
"%CMAKE_EXE%" "%PROJECT_DIR%" -G "Visual Studio 17 2022" -A x64 -DCOOLSOLVE_BUILD_GUI=ON -DCOOLSOLVE_BUILD_TESTS=OFF
if %ERRORLEVEL% neq 0 (
    echo Error during CMake configuration.
    pause
    exit /b %ERRORLEVEL%
)

:: --- Section 5: Compile ---
echo [4/5] Compiling CoolSolve with GUI (Release)...
"%CMAKE_EXE%" --build . --config Release --target coolsolve
if %ERRORLEVEL% neq 0 (
    echo Error during compilation.
    pause
    exit /b %ERRORLEVEL%
)

:: --- Section 6: Generate installer ---
echo [5/5] Generating Windows Installer...
cd "%PROJECT_DIR%"
if exist "%MAKENSIS_EXE%" (
    "%MAKENSIS_EXE%" "%PROJECT_DIR%\coolsolve.nsi"
    if %ERRORLEVEL% neq 0 (
        echo Error during installer generation.
        pause
        exit /b %ERRORLEVEL%
    )
    echo.
    echo SUCCESS! CoolSolve_Installer.exe has been created in:
    echo   %PROJECT_DIR%
) else (
    echo.
    echo WARNING: makensis.exe not found at %MAKENSIS_EXE%
    echo Install NSIS from https://nsis.sourceforge.net and rerun,
    echo or run manually:  makensis "%PROJECT_DIR%\coolsolve.nsi"
)

pause
endlocal
