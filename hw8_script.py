'''
reads in a file of read counts, finds genes of interest, and performs hierarchical clustering on these data;
accepts filename arguments from the commandline for the counts file (option "-i") and a GFF file (option "-g")
'''
# Name:Yuan Shen  
# UCID:YS348
###########################################################################

# using optparse to control commandline input
import optparse

# default value is the file name
parser=optparse.OptionParser()
parser.add_option('-i',action='store',type='string',default='NextGenRaw.txt')
parser.add_option('-g',action='store',type='string',default='athaliana_lab10.gff')

options,args=parser.parse_args()

######################################################################
print 'data loading...'

import pysam

# load gff file
f = open(options.g,'r')  

d=[i.strip().split('\t') for i in f.readlines()]

# 'genes' is a list with only gene records 
genes=[d[i] for i in range(len(d)-1) if d[i][2]=='gene']  



import sys
from math import log
import numpy as np
from hcluster import *
from scipy import stats
import array
from scipy import stats
from math import log
from numpy import mean


# load .txt counts file
filename = options.i  
f = open( filename , 'r' )
lines = f.readlines()
header = lines[0]
lines = lines[1:]
f.close()

genenames=['']*len( lines )
genedata=['']*len( lines )
for i in range( len( lines ) ):
        templine = lines[i].strip().split( '\t' )
        genenames[i] = templine[0]
        genedata[i] = [float( j ) for j in templine[1:]]
        
#####################################################################

print 'filtering data... (depends on your internet condition, this process may take as long as 40 minutes %>_<% )' 

gene_interest_name=[]
gene_interest_data=[]
for i in range( len( genenames ) ):
    for k in range (len(genes)):
		#  extract gene names from gff file
        gene_id=genes[k][8].strip().split(';')[0][3:]   
        if gene_id==genenames[i]:
			# N is gene length
            N=int(genes[k][4])-int(genes[k][3])    
            RPKM=[]
            RPKM_C=[]
            RPKM_T=[]
            Reads=[]
            Reads_C=[]
            Reads_T=[]

            for j in range(4):
				# Rc1 is reads of each gene
                Rc1=genedata[i][j]     
				# T is the total number of reads for each column
                T=sum(u[j] for u in genedata )  
				# RPKM is a list of 4 RPKM values for each gene
                RPKM.append((Rc1/(N/1e3))/(T/1e6))  
                Reads.append(Rc1)
			# 2 RPKM values for the control of each gene
            RPKM_C=RPKM[0:2]  
			# 2 RPKM values for the treatment of each gene
            RPKM_T=RPKM[2:4]  
			# 2 counts for the control of each gene
            Reads_C=Reads[0:2]  
			# 2 counts for the treatment of each gene
            Reads_T=Reads[2:4]  
			
            
			# to avoid the 'math domain error'
            if (mean(Reads_T)==0)or mean(Reads_C)==0:  
				# iterate gene of interest
                if  (mean(Reads_C)!=0) & ((stats.ttest_rel(RPKM_C,RPKM_T))[1]<0.05):  
                        gene_interest_name.append(genenames[i])
                        gene_interest_data.append(genedata[i])
            elif (((stats.ttest_rel(RPKM_C,RPKM_T))[1]<0.05) & ((log(mean(Reads_C)/mean(Reads_T))>1) | (log(mean(Reads_C)/mean(Reads_T))<-1))):
                gene_interest_name.append(genenames[i])
                gene_interest_data.append(genedata[i])

######################################################################
print 'clustering...'

genedatadist = pdist( gene_interest_data , metric = 'correlation' )

genedatalink = linkage( genedatadist , 'average' , 'correlation' )

# filter the linkage matrix
for i in range( len( genedatalink ) ):   
    for j in range( len( genedatalink[i] ) ):
        if abs( genedatalink[i][j] ) < 1e-7:
            genedatalink[i][j] = 0

# the value of t determine the number of clusters, the default value is 0.1 
t=0.1  

while len(set(fcluster( genedatalink , t , criterion = 'distance' )))<20 or len(set(fcluster( genedatalink , t , criterion = 'distance' )))>30:
        cluster_number=len(set(fcluster( genedatalink , t , criterion = 'distance' )))
        if cluster_number<20:
				# a smaller t indicates a greater cluster_numberv
                t=t-0.001  
        else:
                t=t+0.001
        print cluster_number

mycluster=fcluster( genedatalink , t , criterion = 'distance' )

text = 'genename'+' '+'cluster'+'\n'
for i in range( len( genedatalink ) ):
    text += gene_interest_name[i] +' '+ str( mycluster[i] )+'\n'

print text
        
        
#############################################################################



