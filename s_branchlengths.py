#!/bin/env python3

#this script parses phylip outfiles 
#usage: python3 s_branchlengths.py <lineageID>

import sys

lineage=sys.argv[1]
filename=lineage+'_outfile' #PHYLIP outfile name 

file=open(filename, 'r')

##########read through file and collect data##################

entries=[] ##list of lists 
finals=[] ## end branches 
endnode='off'

on='F'
for line in file:
	line_list=line.strip().split(' ')
	if 'requires' in line_list:
		num=float(line_list[-1]) ## record total number SNPs 

	if 'between' in line_list: #activate next stage 
		on='T'

	if on=='T' and line_list[0]!='From' and line_list[0]!='':

		templist=[]

		n=len(line_list)

		endnode='off'

		for i in range(n):
			if line_list[i]!='':
				templist.append(line_list[i]) ##append node numbers and branch lengths 
			if line_list[i][0:3]=='MAB':  ##if line records end of branch 
				endnode='on'

		if endnode=='on':
			finals.append(templist)  ##record of branch ends (aka isolate IDs)
		
		entries.append(templist)	


	elif line_list[0]=='From':
		on='F'

		break


file.close()



#####reformat templist##########
n=len(entries)

newids={} #reformatted templist 

for i in range(2,n):
	one=entries[i][0]
	two=entries[i][1]
	val=entries[i][2]

	newid=one+'_'+two ##1_2, 3_4, etc 

	newids[newid]=val

#print(newids)

########go through data and add branch lengths #####################

branchlengths={} ## dictionary of isolate ID: branchlength

m=len(finals)

for i in range(m): ## start by skipping first 2 entries (index 0 and 1 )
	temp=finals[i]
	start=temp[0]
	iso=temp[1]
	dist=temp[2]

	branchlengths[iso]=float(dist) ##initialize everything 

	if start!='1':
		counter=1 #if start=2, counter=1. 

		while counter!=int(start):
			tester=str(counter)+'_'+start

			if tester in newids:

				branchlengths[iso]+=float(newids[tester])
			counter+=1


#print(branchlengths)

#print('ref:',lineage,',', num, sep='')

SNPs=[]

for item in branchlengths:
	temp=branchlengths[item]*num 
	#print(item, temp, sep=',')  ##gives root-to-tip branch lengths (to node 1) for each isolate, in SNPs

	SNPs.append(temp)

#print(SNPs)

##remove outgroup data from SNPs 
SNPs.remove(max(SNPs))

##dMRCA can be calculated by averaging all branch lengths 
dmrca=sum(SNPs)/len(SNPs)

print(dmrca)


