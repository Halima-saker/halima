from random import randint
import random
from collections import defaultdict
from collections import OrderedDict
from json import loads, dumps
import matplotlib.pyplot as plt
import itertools
from itertools import repeat
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import spline
import numpy
from matplotlib import pyplot
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

#nombre of intervals
N=10

#nombre of datatrack
n=6

#min and max length of interval
min_=300
max_=1000

#generate intervals length
intervals=[random.randrange(min_,max_) for _ in range(0, N)]
#####intervals=[59,252,142,790,213,1474,1208,171,118,338]


#generate datatracks ###datatracks  (random decimal => x="decimal_nbrs", x.split("\t"), [float(e) for e in x]
datatrack=[]
for i in range(n):
	datatrack.append([random.random() for _ in range(0, N)])
#####datatrack=[]
#####datatrack.append([0.8,3.3,6.6,8,7.1,5.9,6.9,4.9,7.5,3.5])
    #####datatrack.append([2.4,4.7,7.8,6.9,5.7,6.8,8.1,8.8,7.6,3.6])



#generate random intervals sequence
alphabet = list('ATGC')
for i in range(1,len(intervals)):
    f = open('C:/Users/Halima/Desktop/phd_script/seq/seq'+str(i), 'w')
    dna = [random.choice(alphabet) for i in range(intervals[i])]
    dna = ''.join(dna)
    f.write(">"+str(i)+"\n"+dna+"\n")
    f.close()


#define start and end for each interval
inter=[]
p=0
for i in range(len(intervals)):
	inter.append([p+1,p+intervals[i]])
	p=p+intervals[i]


print inter   
    

#plot intervals
radius = list(itertools.chain(*inter))
area={}
b=0
#for k in range(len(datatrack)):
    #print k
for i in inter:
    plt.plot(i, [datatrack[0][b], datatrack[0][b]])
    b=b+1
    
#plt.show()
    


    
for k in range(len(datatrack)):
    print (k+1)
    area[k+1]=[x for item in datatrack[k] for x in repeat(item, 2)]

for i in range(len(datatrack)):
    plt.plot(radius, area[i+1], label="data track"+str(i+1))
    plt.xlabel('position')
    plt.ylabel('values')
    plt.title('data track'+str(i+1))
   # plt.show()

test=[]
for i in range(len(intervals)):
	test.append([inter[i],[[row[i] for row in datatrack]]]) 
print test

#noise data
noise={}
noise1={}
noise2={}
dna_pos=range(1,sum(intervals)+1)
radius_noise = dna_pos

for i in range(len(datatrack)):
    noise[i+1]=[(x+ np.random.normal(0,0.1)) for item in range(len(datatrack[i])) for x in repeat(datatrack[i][item], inter[item][1]-inter[item][0]+1)]
    noise1[i+1]=[(x+ np.random.normal(0,0.5)) for item in range(len(datatrack[i])) for x in repeat(datatrack[i][item], inter[item][1]-inter[item][0]+1)]
    noise2[i+1]=[(x+ np.random.normal(0,0.6)) for item in range(len(datatrack[i])) for x in repeat(datatrack[i][item], inter[item][1]-inter[item][0]+1)]
    
  #  for v in intervals:        
    #    noise.append()
    
for i in range(len(datatrack)):
    plt.plot(radius_noise, noise2[i+1], label="data track"+str(i+1)+" with noise")
    plt.xlabel('position')
    plt.ylabel('values')
    plt.axhline(y=0.5, xmin=1000, xmax=100000)
    plt.title('data track'+str(i+1))
    #plt.show()


#smoothing data for datatracks

x = np.array(dna_pos)

for i in range(1,len(datatrack)+1):
    y = np.array(noise[i])
    y1=np.array(noise1[i])
    y2=np.array(noise2[i])
    yhat = savitzky_golay(y, 51 ,2) # window size 73, polynomial order 2
    yhat1=savitzky_golay(y1, 51 ,2)
    yhat2=savitzky_golay(y2, 51 ,2)
    yhat.tofile("C:/Users/Halima/Desktop/phd_script/data1/data_"+str(i)+".txt", sep="\t")
    yhat1.tofile("C:/Users/Halima/Desktop/phd_script/data2/data_"+str(i)+".txt", sep="\t")
    yhat2.tofile("C:/Users/Halima/Desktop/phd_script/data3/data_"+str(i)+".txt", sep="\t")
   
    plt.figure(figsize=(9,6))
    plt.plot(x,y,alpha=.4,label='noisy data '+str(i))
    plt.plot(x,yhat, color='red',label='filtered')
    plt.legend()
    #plt.show()

    plt.figure(figsize=(9,6))
    plt.plot(x,y,alpha=.4,label='noisy data '+str(i))
    plt.plot(x,yhat1, color='red',label='filtered')
    plt.legend()
    #plt.show()

    plt.figure(figsize=(9,6))
    plt.plot(x,y,alpha=.4,label='noisy data '+str(i))
    plt.plot(x,yhat2, color='red',label='filtered')
    plt.legend()
    #plt.show()
	
#exit()

#hh = open('C:/Python27/phd/sequences/values', 'w')

#l = [15, 18, 2, 36, 12, 78, 5, 6, 9]
#print reduce(lambda x, y: x + y, l) / len(l)
lista=[]
seg=5
r=0.002

for i in range(0,len(yhat)):
    
    lista.sort(reverse=False)
   # print lista
    if i<len(yhat)-1 and abs(yhat[i+1]-yhat[i]) > r:
        
        if len(lista)>=seg:
            
            if abs((yhat[i+1]-yhat[i]))>lista[0][0]:
              #  print "yessss"+str(len(lista))
                lista[0]=[abs((yhat[i+1]-yhat[i])),i]
                lista.sort(reverse=False)
        else:
            print "noooo"
            lista.append([abs((yhat[i+1]-yhat[i])),i])
            
               

