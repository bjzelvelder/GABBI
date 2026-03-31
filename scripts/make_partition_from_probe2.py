#! /usr/bin/python3

import sys

def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

probe=str(sys.argv[2])
min_size=int(sys.argv[3])
sortie=sys.argv[4]
fout=open(sortie,'w')

list_good_part=[]
list_bad_part=[]
for name, seq in read_fasta(open(sys.argv[1])):
	if name == ">"+probe:
		c=0
		probe_start=""
		probe_end=""
		after_probe_start=""
		after_probe_end=""
		before_probe_start=""
		before_probe_end=""
		done=""
		for i in seq.upper():
			c+=1
			if c == 1 and (i == "N" or i == "-"):
				before_probe_start=c
			elif c == 1 and i != "N" and i != "-":
				probe_start=c
			elif i != "N" and i != "-" and before_probe_start != "" and after_probe_start == "" and probe_start == "":
				before_probe_end=c-1
				probe_start=c
				
			elif i != "N" and i != "-" and after_probe_start == "" and probe_start != "":
				if c == len(seq):
					probe_end=c
					after_probe_start == ""
					after_probe_end == ""
				continue
				
			elif (i == "N" or i == "-") and after_probe_start == "" and probe_start != "":
				after_probe_start=c
				probe_end=c-1
				done="?"
			elif i != "N" and i != "-" and after_probe_start != "" and done == "?":
				after_probe_start=""
				after_probe_end == ""
				probe_end=c
			elif (i == "N" or i == "-") and (before_probe_start != "" or after_probe_start != "") and c == len(seq):
				if before_probe_start != "" and before_probe_end == "":
					before_probe_end=c
				elif after_probe_start != "":
					after_probe_end=c
			elif i != "N" and i != "-" and after_probe_start != "" and c == len(seq):
				probe_end=c
				after_probe_start == ""
				after_probe_end == ""
		

		if before_probe_start != "" and before_probe_end != "":
			if before_probe_end - before_probe_start + 1 >= min_size:
				part1=str(before_probe_start)+"-"+str(before_probe_end)
				fout.write("part1 = "+part1+"\n")
		if probe_start != "" and probe_end != "":
			part2=str(probe_start)+"-"+str(probe_end)
			fout.write("part2 = "+part2+"\n")
			"""fout.write("DNA, part2_pos1 = "+str(probe_start)+" - "+str(probe_end)+"\\3"+"\n")
			fout.write("DNA, part2_pos2 = "+str(probe_start+1)+" - "+str(probe_end)+"\\3"+"\n")
			fout.write("DNA, part2_pos3 = "+str(probe_start+2)+" - "+str(probe_end)+"\\3"+"\n")"""
		if after_probe_start != "" and after_probe_end != "":
			if after_probe_end - after_probe_start +1 >= min_size:
				part3=str(after_probe_start)+"-"+str(after_probe_end)
				fout.write("part3 = "+part3+"\n")

				
				
"""


c=0
for i in list_bad_part:
	c+=1
	fout.write("DNA, part"+str(c)+" = "+i+"\n")



for i in list_good_part:
	c+=1
	start=int(i.split("-")[0])
	end=int(i.split("-")[1])
	fout.write("DNA, part"+str(c)+" = " + str(start)+ " - " + str(end)+ "\\3\nDNA, part"+str(c+1)+" = " + str(start+1)+ " - " + str(end)+ "\\3\nDNA, part"+str(c+2)+" = " + str(start+2)+ " - " + str(end)+ "\\3\n")
	c=c+2	
"""
					
