cat list_filenames.txt | tr -d '\r' > tmp && mv tmp list_filenames.txt #if text file has been modified in Windows, line endings are different and problems
xargs -a list_filenames.txt -n 2 mv

echo 'Files renamed!'
