#from __future__ import division
import sys
import numpy as np					
import matplotlib as mpl
import matplotlib.pyplot as plt

def convert2float(raw_data):
	i = 0; 
	j = 0;
	rows = raw_data.count('\n')
	columns = (raw_data.count('\t')/rows)
	data = np.zeros(shape=(rows,columns))
	for line in (raw_data.split('\n')[:-1]):
		for word in (line.split('\t')[:-1]):
			data[i][j] = float(word)
			j += 1
		i += 1
		j = 0;
	return data
	

life_raw = open('lifetime.dat','r')
life_data = convert2float(life_raw.read())
print 'Life data loaded!'

#change_raw = open('change.dat','r')
#change_data = convert2float(change_raw.read())
#print 'Change data loaded!'

plt.figure(1)
ax = plt.axes([0,0,1,1], frameon=False)
ax.set_axis_off()
plt.imshow(life_data,interpolation='nearest',aspect='auto',cmap=mpl.cm.get_cmap('RdBu_r'))

plt.figure(2)
ax = plt.axes([0,0,1,1], frameon=False)
ax.set_axis_off()
plt.imshow(life_data,interpolation='nearest',aspect='auto',cmap=mpl.cm.get_cmap('RdBu'))

#plt.figure(2)
#plt.plot(change_data[:,0],change_data[:,1],'.')

plt.show()



