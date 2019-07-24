for i in *fcs; do mv "$i" "${i// /_}"; done

ls *fcs > list_filenames.txt
echo "File generated: list_filenames.txt"
