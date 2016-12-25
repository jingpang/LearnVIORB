import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

filepath = '/home/jp/opensourcecode/OpenSourceORBVIO/tmp/';
filename = filepath+'KeyFrameNavStateTrajectory.txt';

NS = np.loadtxt(filename);

time = NS[:,0];

fig = plt.figure(1);
ax = fig.add_subplot(111, projection='3d');
ax.plot(NS[:,1],NS[:,2],NS[:,3]);

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Position in initial frame');
plt.savefig(filepath+'Pw.eps', format="eps");


bg = NS[:,11:14];
ba = NS[:,14:17];
dbg = NS[:,17:20];
dba = NS[:,20:23];

fig2 = plt.figure(2);
p21, =plt.plot(time-time[0],(bg+dbg)[:,0]);
p22, =plt.plot(time-time[0],(bg+dbg)[:,1]);
p23, =plt.plot(time-time[0],(bg+dbg)[:,2]);
plt.legend([p21,p22,p23],["bgx","bgy","bgz"]);
plt.title('gyr bias');
plt.savefig(filepath+'biasgyr.eps', format="eps");

fig3 = plt.figure(3);
p31, =plt.plot(time-time[0],(ba+dba)[:,0]);
p32, =plt.plot(time-time[0],(ba+dba)[:,1]);
p33, =plt.plot(time-time[0],(ba+dba)[:,2]);
plt.legend([p31,p32,p33],["bax","bay","baz"]);
plt.title('acc bias');
plt.savefig(filepath+'biasacc.eps', format="eps");

v = NS[:,4:7];
fig4 = plt.figure(4);
p41, =plt.plot(time-time[0],v[:,0]);
p42, =plt.plot(time-time[0],v[:,1]);
p43, =plt.plot(time-time[0],v[:,2]);
plt.legend([p41,p42,p43],["x","y","z"]);
plt.title('Velocity in initial frame');
plt.savefig(filepath+'Vw.eps', format="eps");

plt.show();

#biasa = np.loadtxt(filepath+'biasa.txt');
#plt.figure(1);
#p11, =plt.plot(biasa[:,0]-biasa[0,0],biasa[:,1]);
#p12, =plt.plot(biasa[:,0]-biasa[0,0],biasa[:,2]);
#p13, =plt.plot(biasa[:,0]-biasa[0,0],biasa[:,3]);
#plt.title('bias-acc');
#plt.legend([p11,p12,p13],["x","y","z"]);
#plt.savefig(filepath+"biasa.eps", format="eps")
##plt.legend(p12,'y');
##plt.legend(p13,'z');

#scale = np.loadtxt(filepath+'scale.txt');
#plt.figure(2);
#[p21,p22] = plt.plot(scale[:,0]-scale[0,0],scale[:,1:3]);
#plt.title('scale');
#plt.legend([p21,p22],['aftopt','befopt']);
#plt.savefig(filepath+'/scale.eps', format="eps")

#condnum = np.loadtxt(filepath+'condnum.txt');
#plt.figure(3);
#plt.plot(condnum[:,0]-condnum[0,0],condnum[:,1]/condnum[:,6]);
#plt.title('condnum');
#plt.savefig(filepath+'condnum.eps', format="eps")


#biasg = np.loadtxt(filepath+'biasg.txt');
#plt.figure(4);
#p41, =plt.plot(biasg[:,0]-biasg[0,0],biasg[:,1]);
#p42, =plt.plot(biasg[:,0]-biasg[0,0],biasg[:,2]);
#p43, =plt.plot(biasg[:,0]-biasg[0,0],biasg[:,3]);
#plt.title('bias-gyr');
#plt.legend([p41,p42,p43],["x","y","z"]);
#plt.savefig(filepath+"biasg.eps", format="eps")

#plt.show();
