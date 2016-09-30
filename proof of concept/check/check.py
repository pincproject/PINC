#!/usr/bin/env python

import sys
import re
import tempfile
import os

"""
This script checks for certain common violations of coding practices, and
suggests correcting them.

The patterns to match for using regular expressions (regex) is given in the
variable 'pattern', one pattern to match against per line. A suggested
replacement is given in 'replacement' and a human readable message of what's
wrong is given in 'message'. The lines of these three variables correspond.
Since regexes are typically quite incomprehensible, it's best to read first in
'message' what it does.

For this script to work, each line in 'pattern' must have exactly one capturing
group (this is a limitation of the script)
"""

pattern = re.compile(	"( \t)"
						"|(^    )"
						"|([\t ]$)"
						"|(?://(?P<a>[^ /@<]))"
						"|(?:/\*(?P<b>[^ \*\n]))"
						"|(//[ ]{2,})"
						"|(/\*[ ]{2,})")

replacement = [	"\t",
				"\t",
				"",
				"// \g<a>",
				"/* \g<b>",
				"// ",
				"/* "]

message = [	"space before hard tab",
			"soft tab (spaces) instead of hard tab",
			"white space at end of line",
			"// not followed by space, /, @ or <",
			"/* not followed by space, * or newline",
			"// followed by multiple spaces",
			"/* followed by multiple spaces"]

def repl(match):
	index = match.lastindex-1
	return match.expand(replacement[index])

for fileName in sys.argv[1:]:
	file = open(fileName,'r')
	mistakes = 0

	for i, line in enumerate(file):
		for match in re.finditer(pattern, line):

			# Matched pattern group, unescaping escaped characters.
			# matched = match.group().encode('unicode-escape').decode()
			num = match.lastindex-1

			print("%s:%d:%d: warning: %s (%s)"%
				(fileName,i+1,match.start(),message[num],line[:-1].strip()))

			mistakes += 1

	file.close()

	print("%d mistakes in %s"%(mistakes,fileName))

	if mistakes==0:
		fix = 'NO'
	else:
		fix = raw_input("Auto-fix %s? (Y)es, (N)o or "
						"(S)tep through each change (default)? "%fileName)

	if not fix.upper() in ['N','NO']:
		file = open(fileName,'r')
		tempFile = open(fileName+".autofix",'w')

		if fix.upper() in ['Y','YES']:
			for i, line in enumerate(file):
				newLine = pattern.sub(repl,line)
				tempFile.write(newLine)
		else:
			for i, line in enumerate(file):

				newLine = pattern.sub(repl,line)

				if line != newLine:

					match = re.search(pattern,line)
					num = match.lastindex

					rep = raw_input("Replace (%s):\n\t%swith:\n\t%s"
									"(Y)es (default) or (N)o? "
									%(message[num],line,newLine))

					if rep.upper() in ['N','NO']:
						tempFile.write(line)
					else:
						tempFile.write(newLine)
				else:
					tempFile.write(line)

		tempFileName = tempFile.name
		tempFile.close()
		file.close()

		print("All mistakes fixed")
		overwrite = raw_input(	"(O)verwrite %s by auto-fixed file (default), "
								"store to (N)ew file or (K)eep backup of old "
								"file? "%fileName)

		if overwrite.upper() in ['N','NEW']:
			pass
		elif overwrite.upper() in ['K','KEEP']:
			os.rename(filename,filename+".pre-autofix")
			os.rename(tempFileName,fileName)
		else:
			os.remove(fileName)
			os.rename(tempFileName,fileName)
