'''
Consider last week’s as signment with the counts data in NextGenRaw.txt, specifically the filtering step. Create a list of the p-values (from t-tests) and log ratios for every gene. Create a volcano plot (scatter plot) of the log ratios vs the p-values. Add lines onto the plot that indicate the .05 p-value cutoff and log ratios greater than 1 or less than -1 (simply adding lines is sufficient). Make sure to add a title and label the axes!
'''

# Name:Yuan Shen  
# UCID:YS348
###########################################################################

# using optparse to control commandline input
import optparse

# default value is file name
parser=optparse.OptionParser()
parser.add_option('-i',action='store',type='string',default='NextGenRaw.txt')
parser.add_option('-g',action='store',type='string',default='athaliana_lab10.gff')

options,args=parser.parse_args()

##########################################################################
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
        
#####################################################################################33

print 'filtering data... (depends on your internet condition, this process may take as long as 40 minutes %>_<% )' 

gene_p_value=[]
gene_log_ratio=[]
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
				# T is the total number of reads of each column
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
                    pass
                
                    
            else:
                gene_p_value.append(stats.ttest_rel(RPKM_C,RPKM_T)[1])
                gene_log_ratio.append((log(mean(Reads_C)/mean(Reads_T))))

import matplotlib.pyplot as plt

# range of table
plt.axis([-4,4,0,3])   

# create the plots
plt.plot(gene_log_ratio,gene_p_value,'ro')   

# add cutoff line with y=0.05
plt.plot([-4,4],[0.05,0.05],'k-')    

# add cutoff line with x=-1
plt.plot([-1,-1],[0,3],'k-')      

# add cutoff line with x=1
plt.plot([1,1],[0,3],'k-')       

plt.xlabel('log2_ratio')

plt.ylabel('P_value')

plt.savefig('hw9_fig1'+'.pdf')

plt.title('counts data in NextGenRaw')

plt.show()

