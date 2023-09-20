Protein_Seq= input("Put in the Protein Sequence:")

filename = 'OneLetterCode.dat'
fi = open (filename, 'w')
for index,letter in enumerate(Protein_Seq,1):
	print(letter, file = fi)
fi.close()
