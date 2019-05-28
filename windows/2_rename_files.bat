:: This script takes list_filenames.txt with new names added an rename all files.
ECHO OFF
for /F "tokens=*" %%A in (list_filenames.txt) do ren %%A
::for /F "tokens=1,2 delims=," %%g in (list_filenames.txt) do ren %%g %%h
ECHO 'Files renamed!'
PAUSE