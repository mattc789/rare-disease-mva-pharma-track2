@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -STA -ExecutionPolicy Bypass -File "%~dp0JURY_START_V31.ps1" -Action menu
if errorlevel 1 (
  echo.
  echo The jury launcher returned an error. Review the message above.
  pause
)
endlocal
