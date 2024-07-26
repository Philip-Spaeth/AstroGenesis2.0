@echo off
setlocal

REM Check if ffmpeg is installed
ffmpeg -version >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo ffmpeg is not installed or not in the system path. Please install ffmpeg to use this script.
    exit /b 1
)

REM Set the framerate and output file
set "framerate=30"
set "output_file=output_video.mp4"

REM Ensure the BMP files are numerically ordered
set "input_pattern=Picture_%%d.bmp"

REM Convert all BMP files to a video with high quality and maintain aspect ratio
ffmpeg -framerate %framerate% -i %input_pattern% -vf "scale=-2:1080:flags=lanczos" -c:v libx264 -preset slow -crf 18 -pix_fmt yuv420p %output_file%

REM Delete each BMP file after use
REM for %%f in (Picture_*.bmp) do del %%f

echo Conversion complete.
pause
