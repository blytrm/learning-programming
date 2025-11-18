# print weirdly named files
	cat ./- 
	cat -- ./"--spaces in this filename--"

# find human readable file
	file ./-*
	cat ./"_"

# find file: human-readable, 1033 bytes in size, not executable
	find . -type f -size 1033c ! -executable

# find file: owned by user bandit7 + owned by group bandit6 + 33 bytes in size
find / -type f -user bandit7 -group bandit6 -size 33c 2>/dev/null 
	# 2>/dev/null redirects error messages > cleaner output
	cat /var/lib/dpkg/info/bandit7.password
		morbNTDkSW6jIlUc0ymOdMaLnOlFVAaj

# lvl 7 : The password for the next level is stored in the file data.txt next to the word millionth
du -b data.txt # file is huge > cant print
cat data.txt | grep "millionth" # prints the line
# dfwvzFQi4mU0wfNbFOe9RoWskMLg7eEc

# lvl 8: The password for the next level is stored in the file data.txt and is the only line of text that occurs only once
sort data.txt | uniq -c # can see the answer but also prints other occurances
sort data.txt | uniq -u # single occurances
	# 4CKMh1JI91bUIZZPXDqGanal4x

# lvl 9: The password for the next level is stored in the file data.txt in one of the few human-readable strings, preceded by several ‘=’ characters.
strings data.txt | grep "==*" # strings = finds human readable strings in files
	# FGUW5ilLVJrxX9kMYMmlN4MgbpfMiqey
