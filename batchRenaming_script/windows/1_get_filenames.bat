:: This script takes all filenames with *FCS pattern and store them in a text file
ECHO OFF
setlocal enabledelayedexpansion

set "pattern= "
set "replace=_"

for %%I in (*.fcs) do (
    set "file=%%~I"
    ren "%%I" "!file:%pattern%=%replace%!"
)

dir /b *.FCS > list_filenames.txt
ECHO 'File generated: list_filenames.txt'
PAUSE