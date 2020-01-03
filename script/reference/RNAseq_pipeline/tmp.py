
Out={}
with open('D:/my documents/refGene.txt') as refFile:
  for line in refFile:
    col=line.strip().split('\t')
    #GeneLength= txEnd- txStart+1
    GeneLen=abs(int(col[5])-int(col[4]))+1
    #get every start or end location of each exon
    exonStarts=col[9].strip().split(',')
    exonEnds=col[10].strip().split(',')
    #calculate the length of full transcript(the sum of exon length)
    exonLen=0
    for i in range(0,len(exonStarts)-1):
      exon=abs(int(exonEnds[i])-int(exonStarts[i]))+1
      exonLen += exon
    #calculate the sum of intron length
    intronLen=0
    for i in range(0,len(exonStarts)-2):
      intron= abs(int(exonStarts[i+1])-int(exonEnds[i]))-1
      intronLen += intron
   #output as dictionary format
    Out[col[1]] =[GeneLen,exonLen,intronLen]
#Output the results   
OutFile=open('RefGene_Stat.txt','w')
#add header to the file
OutFile.write('RefSeqID'+'\t'+'GeneLen'+'\t'+'TxLen'+'\t'+'IntronLen'+'\n')
for key in Out.keys():
  OutFile.write(key+'\t'+str(Out[key][0])+'\t'+str(Out[key][1])+'\t'+str(Out[key][2])+'\n')
  
OutFile.close()