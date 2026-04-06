!include "MUI2.nsh"

Name "CoolSolve"
OutFile "CoolSolve_Installer.exe"
InstallDir "$PROGRAMFILES64\CoolSolve"
InstallDirRegKey HKCU "Software\CoolSolve" ""

RequestExecutionLevel admin

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "LICENSE"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

!insertmacro MUI_LANGUAGE "English"

Section "CoolSolve" Section_Main
    SetOutPath "$INSTDIR"
    File "build\Release\coolsolve.exe"
    ; If we had external DLLs, they would go here

    ; Copy example models
    SetOutPath "$INSTDIR\examples"
    File "examples\*.eescode"
    File "examples\*.initials"
    SetOutPath "$INSTDIR"
    
    WriteUninstaller "$INSTDIR\Uninstall.exe"
    
    ; Add to Start Menu
    CreateShortcut "$SMPROGRAMS\CoolSolve.lnk" "$INSTDIR\coolsolve.exe" "" "$INSTDIR\coolsolve.exe" 0
    
    ; File Association for .eescode files
    WriteRegStr HKCR ".eescode" "" "CoolSolve.EESCode"
    WriteRegStr HKCR "CoolSolve.EESCode" "" "CoolSolve Equation Set"
    WriteRegStr HKCR "CoolSolve.EESCode\DefaultIcon" "" "$INSTDIR\coolsolve.exe,0"
    WriteRegStr HKCR "CoolSolve.EESCode\shell\open\command" "" '"$INSTDIR\coolsolve.exe" --gui "%1"'
    
    ; Add to Add/Remove Programs
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoolSolve" "DisplayName" "CoolSolve"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoolSolve" "UninstallString" '"$INSTDIR\Uninstall.exe"'
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoolSolve" "DisplayIcon" "$INSTDIR\coolsolve.exe"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoolSolve" "Publisher" "CoolSolve"
SectionEnd

Section "Uninstall"
    Delete "$INSTDIR\coolsolve.exe"
    Delete "$INSTDIR\Uninstall.exe"
    RMDir /r "$INSTDIR\examples"
    RMDir "$INSTDIR"
    
    Delete "$SMPROGRAMS\CoolSolve.lnk"
    
    DeleteRegKey HKCR ".eescode"
    DeleteRegKey HKCR "CoolSolve.EESCode"
    DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoolSolve"
SectionEnd
