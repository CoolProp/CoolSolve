@echo off
setlocal

:: ============================================================
:: CoolSolve Windows Build + Installer Script
::
:: Requirements:
::   - Visual Studio 2022 with C++ Desktop Development workload
::   - Node.js (https://nodejs.org) for the React GUI frontend
::   - Python 3 (e.g. Anaconda) — needed by CoolProp's CMake scripts
::   - NSIS (https://nsis.sourceforge.net) for the installer
::
:: Update the paths in Section 1 if your installation differs.
:: ============================================================

:: --- Section 1: Configuration ---
:: Path to your Anaconda installation. Used to activate the conda env
:: below AND as one of the Python lookup locations. Harmless if absent.
set "ANACONDA_DIR=E:\Anaconda"
set "CONDA_ENV=fast_numpy_env"

:: CMake shipped with Visual Studio 2022 (Build Tools or Community).
set "CMAKE_EXE=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"

:: NSIS installer compiler.
set "MAKENSIS_EXE=C:\Program Files (x86)\NSIS\makensis.exe"
if not exist "%MAKENSIS_EXE%" set "MAKENSIS_EXE=C:\Program Files\NSIS\makensis.exe"

:: Override PYTHON_EXE to force a specific Python interpreter, e.g.:
::   set "PYTHON_EXE=C:\Anaconda3\python.exe"
:: Otherwise the script auto-discovers Python (see :find_python below).

:: Root of the project (directory containing this bat file)
set "PROJECT_DIR=%~dp0"
set "PROJECT_DIR=%PROJECT_DIR:~0,-1%"

:: --- Section 2: Environment setup (Anaconda PATH + Python discovery) ---
echo [1/5] Setting up environment...

:: Prepend Anaconda to PATH (no-op if ANACONDA_DIR does not exist).
set "PATH=%ANACONDA_DIR%;%ANACONDA_DIR%\Scripts;%ANACONDA_DIR%\Library\bin;%PATH%"

:: Activate the conda environment if Anaconda is installed (sets
:: CONDA_PREFIX etc.). Harmless if activate.bat is missing.
if exist "%ANACONDA_DIR%\Scripts\activate.bat" (
    call "%ANACONDA_DIR%\Scripts\activate.bat" %CONDA_ENV%
) else (
    echo WARNING: Anaconda not found at "%ANACONDA_DIR%" — relying on Python auto-discovery.
)

:: CoolProp's CMake calls find_package(Python) at configure time and
:: aborts if no interpreter is found. Resolve Python once here and pass
:: it explicitly via -DPython_EXECUTABLE so the build does not silently
:: depend on the right conda env being active.
call :find_python
if errorlevel 1 (
    echo.
    echo ERROR: Python interpreter not found. CoolProp's CMake scripts
    echo        require Python 3 at configure time.
    echo.
    echo        I tried, in order:
    echo          1. The PYTHON_EXE environment variable ^(full path override^)
    echo          2. "python" on PATH
    echo          3. "%ANACONDA_DIR%\envs\%CONDA_ENV%\python.exe"
    echo          4. "%ANACONDA_DIR%\python.exe"
    echo          5. Common Anaconda / Miniconda / Python.org install locations
    echo.
    echo        Fix: set the PYTHON_EXE environment variable to the full
    echo        path of your python.exe, then re-run this script. Example:
    echo.
    echo          set "PYTHON_EXE=C:\Anaconda3\python.exe"
    echo          build_installer.bat
    echo.
    echo        Or install Python 3 from https://www.python.org/downloads/
    echo        or Anaconda from https://www.anaconda.com/download/
    echo.
    pause
    exit /b 1
)
echo        Using Python: "%PYTHON_EXE%"
"%PYTHON_EXE%" --version

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
call npm ci
if %ERRORLEVEL% neq 0 (
    echo Error running npm ci.
    pause
    exit /b %ERRORLEVEL%
)
call npm run build
if %ERRORLEVEL% neq 0 (
    echo Error building the GUI frontend.
    pause
    exit /b %ERRORLEVEL%
)
:skip_npm

:: --- Section 4: Configure CMake ---
:: -DPython_EXECUTABLE makes FindPython use the interpreter resolved above,
:: independent of which conda env (if any) is currently active.
echo [3/5] Configuring CMake...
if not exist "%PROJECT_DIR%\build" mkdir "%PROJECT_DIR%\build"
cd "%PROJECT_DIR%\build"
"%CMAKE_EXE%" "%PROJECT_DIR%" -G "Visual Studio 17 2022" -A x64 -DCOOLSOLVE_BUILD_GUI=ON -DCOOLSOLVE_BUILD_TESTS=OFF -DPython_EXECUTABLE="%PYTHON_EXE%"
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
    echo WARNING: makensis.exe not found at "%MAKENSIS_EXE%"
    echo Install NSIS from https://nsis.sourceforge.net and rerun,
    echo or run manually:  makensis "%PROJECT_DIR%\coolsolve.nsi"
)

pause
endlocal
exit /b 0


:: ============================================================
:: Subroutine: locate the Python interpreter.
:: Sets the PYTHON_EXE variable and returns 0 on success, 1 on failure.
:: Search order:
::   1. The PYTHON_EXE environment variable (user override, full path)
::   2. "python" on PATH (via "where python" — benefits from the PATH
::      prepend done in Section 2)
::   3. %ANACONDA_DIR%\envs\%CONDA_ENV%\python.exe  (the conda env)
::   4. %ANACONDA_DIR%\python.exe                    (Anaconda base)
::   5. Common Anaconda / Miniconda / Python.org install paths
:: ============================================================
:find_python

:: 1. User-provided override (must be a full path that exists).
if defined PYTHON_EXE if exist "%PYTHON_EXE%" exit /b 0

:: 2. "python" found on PATH (first match wins).
for /f "delims=" %%i in ('where python 2^>nul') do (
    set "PYTHON_EXE=%%i"
    exit /b 0
)

:: 3. The conda env named by CONDA_ENV under ANACONDA_DIR.
if exist "%ANACONDA_DIR%\envs\%CONDA_ENV%\python.exe" (
    set "PYTHON_EXE=%ANACONDA_DIR%\envs\%CONDA_ENV%\python.exe"
    exit /b 0
)

:: 4. Anaconda base environment.
if exist "%ANACONDA_DIR%\python.exe" (
    set "PYTHON_EXE=%ANACONDA_DIR%\python.exe"
    exit /b 0
)

:: 5. Common install locations for Anaconda, Miniconda, and python.org.
for %%P in (
    "C:\Anaconda\python.exe"
    "C:\Anaconda3\python.exe"
    "C:\ProgramData\Anaconda3\python.exe"
    "C:\ProgramData\miniconda3\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python313\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
    "%USERPROFILE%\Anaconda3\python.exe"
    "%USERPROFILE%\miniconda3\python.exe"
    "C:\Python313\python.exe"
    "C:\Python312\python.exe"
    "C:\Python311\python.exe"
    "C:\Python310\python.exe"
) do (
    if exist %%~P (
        set "PYTHON_EXE=%%~P"
        exit /b 0
    )
)

:: Not found anywhere.
exit /b 1
