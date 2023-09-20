 # Convert One Letter Amino Acid Sequence to Three Letter Sequence

seq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP'} 

oneletter_seq = open("OneLetterCode.dat", "r").read()
seq_length = len(oneletter_seq)
print ("One Letter Sequence:", oneletter_seq)
print ("Length of the protein:", seq_length)

threeletter_seq = open("ThreeLetterCode.dat", "w")
for i in range (seq_length):
    seq_key = oneletter_seq[i]
    if (seq_key != '\n'):
        print (i+1, seq_key, seq[seq_key])
        print (seq[seq_key], file = threeletter_seq)

threeletter_seq.close()
print ("Conversion Complete!")
