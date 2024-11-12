@echo off
setlocal

REM Pfad des Verzeichnisses, in dem die Batch-Datei liegt
set "WIN_PATH=%~dp0"

REM Entfernt den abschließenden Backslash, falls vorhanden
if "%WIN_PATH:~-1%"=="\" set "WIN_PATH=%WIN_PATH:~0,-1%"

REM Manuelle Konvertierung des Windows-Pfads in einen WSL-Pfad
set "WSL_PATH=/mnt/%WIN_PATH:~0,1%%WIN_PATH:~2%"
set "WSL_PATH=%WSL_PATH:\=/%"

REM Konvertiere Laufwerksbuchstaben zu Kleinbuchstaben (optional)
set "WSL_DRIVE_LETTER=%WSL_PATH:~5,1%"
for %%A in ("%WSL_DRIVE_LETTER%") do set "WSL_DRIVE_LETTER=%%~lA"
set "WSL_PATH=/mnt/%WSL_DRIVE_LETTER%%WSL_PATH:~6%"

REM Ausgabe des WSL-Pfads zur Überprüfung
echo WSL_PATH ist: "%WSL_PATH%"

REM Öffnet VSCode im WSL-Modus im Verzeichnis der Batch-Datei
wsl.exe code "%WSL_PATH%"

if errorlevel 1 (
    echo Fehler beim Öffnen von Visual Studio Code.
) else (
    echo Visual Studio Code im WSL-Modus erfolgreich gestartet.
)

endlocal
