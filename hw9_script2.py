
# using optparse to control commandline input
import optparse

# default value equals to each corresponding file name
parser=optparse.OptionParser()
parser.add_option('-i',action='store',type='string',default='NextGenClusters.txt')
options,args=parser.parse_args()

# load the raw data
f = open( options.i , 'r' )
lines = f.readlines()
header = lines[0]
lines = lines[1:]
f.close()
cluster_no=[]
for i in lines:
    cluster_no.append(int(i[10:12]))
cluster_no=max(cluster_no)

# to select unique element from a list
clusters=[]   
for i in range(cluster_no):
    clusters.append([])
for i in range( len( lines ) ):
    templine = lines[i].strip().split( ' ' )
    clusters[int(templine[1])-1].append(templine)   

###########################

import numpy as np
import pylab as plt



for i in range(cluster_no):
    plt.figure()
    plt.hold = True
    boxes=[]
    for k in np.arange(4):
        col=[]
        for j in clusters[i]:
            col.append(float(j[k+2]))
        boxes.append(col)
    plt.boxplot(boxes,vert=0)
    plt.yticks([1,2,3,4],['Control1SE','Control2SE','Nitrate1SE','Nitrate2SE'])
    plt.xlabel('Counts')
    plt.ylabel('Experiment')
    plt.title('Boxplots of Cluster'+str(i+1))
    plt.savefig('hw9_cluster'+str(i+1)+'.pdf')

                       



######################





