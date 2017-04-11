import string, re

# Written by Robert Belshaw (School of Biomedical & Healthcare Sciences, University of Plymouth) & Aris Katzourakis (Department of Zoology, University of Oxford)
# For more information and to cite see Belshaw, R & Katzourakis, A (2005) BlastAlign: a program that uses blast to align problematic nucleotide sequences. Bioinformatics 21:122-123.
# Please send any comments to robert.belshaw@plymouth.ac.uk or aris.katzourakis@zoo.ox.ac.uk

file = open('blast_out', 'r')
buffer = file.readlines()

def Calculate_hits():
	Number_of_landmarks = len(Permanent_dictionary[KeyList[0]]) # use legth of first entry
	counter = 1
	while counter < Number_of_landmarks: # Less than because list starts from zero
		number_of_hits = 0
		for item in KeyList:
			list = Permanent_dictionary[item]
			landmark = list[counter]
			if landmark != '*':
				number_of_hits = number_of_hits + 1
		List_of_hits.append(number_of_hits)
		counter = counter +1
	return List_of_hits
	
def doInsertRoutine(list, value):
	no_ast = 0
	old_diff = 0	
	switch = 0	
	for item in list:
		if item == '*':		
			no_ast = no_ast+1
		else:
			new_diff = (item - value)*(item - value)
			if item < value:
				no_ast = 0
			else:
				i = list.index(item)
				if new_diff > old_diff:
					i = i-no_ast
					list.insert(i, value)
				else:
					list.insert(i, value)
				switch = 1
				break
		old_diff	= new_diff		
	if switch == 0:
		no_ast = 0 
		for item in list:
			if item == '*':		
				no_ast = no_ast+1
			else:
				no_ast = 0
		i = len(list) - no_ast # Finds position before any trailing asterisks
		list.insert(i, value)			
	return list, i

def go_through_Library(Library_dictionary, tempKey, LandmarkInsertPos):
	tempKeyList = []
	for item in KeyList:		
		tempKeyList.append(item)
	tempKeyList.remove(tempKey)
	for item in tempKeyList:
		tempList = []
		for subitem in Permanent_dictionary[item]:
			tempList.append(subitem)
		if Library_dictionary.has_key(item):
			tempList.insert(LandmarkInsertPos, Library_dictionary[item])
			Permanent_dictionary[item] = tempList
		else:
			tempList.insert(LandmarkInsertPos, '*')
			Permanent_dictionary[item] = tempList

def process_previous_block(tempKey, tempValue, Library_dictionary):
	landmark = 0
	tempList = []
	for item in (Permanent_dictionary[tempKey]):
		tempList.append(item)
	for item in (Permanent_dictionary[tempKey]):
		if item != '*':
			if (tempValue >= item-30) and (tempValue <= item+30):
				landmark = 1
			else:
				pass
	if landmark == 0:
		theAnswer = doInsertRoutine(tempList, tempValue)
		tempList  = theAnswer[0]
		LandmarkInsertPos = theAnswer[1]	
		Permanent_dictionary[tempKey] = tempList
		go_through_Library(Library_dictionary, tempKey, LandmarkInsertPos)

def makeOutFile():
	theOutFile = open('blast_out_python', 'w')
	theOutFile.write('\t\t') # Just to line up entries for ease of viewing
	for item in List_of_hits:
		theOutFile.write('%s\t' %item)
	theOutFile.write('\n')
	for item in KeyList:
		theOutFile.write('%s\t' %item)
		for listItem in Permanent_dictionary[item]:
			theOutFile.write('%s\t' %listItem)
		theOutFile.write('\n')

Query_dictionary = {}
Library_dictionary = {}
Permanent_dictionary = {}
KeyList = []
list = [0]
List_of_hits = [] # To note whether entries are unique or not

for line in buffer:
	if line[0] == '*':
		entry = ""
		entry = line[1:-1]
		Permanent_dictionary[entry] = list
		KeyList.append(entry)
n=0
previousKey = "null" # Need in case have identical sequences & then need to avoid unassigned variable

for line in buffer:
	tempList = []
	if line[0:5] == 'Query':
		if n >= 1:
			process_previous_block(QueryKey, QueryValue, Library_dictionary)		
		Library_dictionary = {}
		line = string.split(line)
		QueryKey = (line[0])[5:]
		QueryValue = string.atoi(line[1])
		Query_dictionary[QueryKey] = QueryValue
		n=n+1
	elif line[0:7] == 'Library':
		line = string.split(line)
		LibraryKey = (line[0])[7:]
		LibraryValue = string.atoi(line[1])
		if LibraryKey != QueryKey:
			if previousKey == LibraryKey:	
				previousDist = (previousValue-QueryValue)*(previousValue-QueryValue)
				currentDist = (LibraryValue-QueryValue)*(LibraryValue-QueryValue)
				if currentDist < previousDist:
					Library_dictionary[LibraryKey] = LibraryValue
			else:
				Library_dictionary[LibraryKey] = LibraryValue			
		previousKey = (line[0])[7:]
		previousValue = string.atoi(line[1])
	
Calculate_hits() 
makeOutFile()
		

