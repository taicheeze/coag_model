#%%
import matplotlib.pyplot as plt
import numpy
import subprocess
n=58
subprocess.call('sed -i "s/r_constant\[{0}\] =.*/r_constant[{0}] = 2.1000e-03;/g" ../VolumeReactions.h'.format(n), shell = True)
subprocess.call("make -C ../ coagulationo3", shell = True)
subprocess.call("make -C ../ runcoagulation", shell = True)
times = []
meanconc = []
for i in range (1,301):
    f = open('../output/concentrations_{0}.m'.format(i))

    f.readline()

    times.append(float(f.readline()))

    f.readline()
    f.readline()
    vals = []
    vals2 = []

    for line in f.readlines():
        try:
            val = float(line)
            vals.append(val)
        except:
            break

    f.close()

    vals = numpy.array(vals)
    vals2 = numpy.reshape(vals, (25, 50))
    vals3 = vals2[0:15, 9:20]
    meanval = numpy.mean(vals3)
    meanconc.append(meanval)
maxval=numpy.amax(meanconc)
float(maxval)
endval=meanconc[-1]
float(endval)
f = open('./outfile.dat','a')
if maxval > 2e-3 and endval < maxval/3:
    f.write("1\n")
else:
    f.write("0\n")
f.close()
#plt.imshow(vals2)
