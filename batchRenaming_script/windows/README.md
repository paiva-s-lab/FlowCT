# Batch file renaming (for Windows) 

**IMPORTANT!** These two files (`1_get_filenames.bat` + `2_rename_files.bat`) **MUST** be within the folder with all FCS files to be processed.

1. Run `1_get_filenames.bat`. This file performs two operations:
	- Replaces all spaces in filenames with underscores.
	- Creates a text file called _list_filenames.txt_ with all FCS files whithin this folder.

2. Open the text file in Excel and, in the next column to that one with filenames, write the new name for each FCS file. Do not create a header in Excel. Save it as txt (without change its name)!!

			filename1	newname1
			filename2	newname2
					...
					...

3. Run `2_rename_files.bat` and automatically all files will be renamed to new names.
