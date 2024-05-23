#*************************************************************************** 
#
# created by: Nathan Alexander (March 2024)
#
# The purpose of this program is to extract and then animate numerical data
#
#***************************************************************************
     
import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.display import HTML
from io import StringIO
import pandas as pd

#-- extract constants from fort.q000
f = np.genfromtxt(r'/home/npalex/Lid-driven-cavity/_output/fort.q0000', usecols=np.arange(0,1))

#-- define constants
mx = int(f[2])                # total number of cells in x
my = int(f[3])                # total number of cells in y

xlow = f[4]                   # lower limit of spatial domain
ylow = f[5]
dx = f[6]                     # grid spacing in x
dy = f[7]                     # grid spacing in y

print("dx = ", dx)

#-- extract the number of output times from claw.data
g = np.genfromtxt(r'/home/npalex/Lid-driven-cavity/_output/claw.data', usecols=np.arange(0,1))
steps = int(g[9])                                 # number of time steps
meqn = 4                                          # number of governing equations

print("steps = ", steps)

#-- print parameters
with open(r'/home/npalex/Lid-driven-cavity/_output/setprob.data', 'r') as f:
    print(f.read())

#-- initialize arrays 
q = np.zeros((meqn, steps+1, mx, my))                  

#-- exract data from files for each time step
for k in range(0, meqn):
    for i in range(0, steps+1):
       
        if i<10:
            f = pd.read_csv(r'/home/npalex/Lid-driven-cavity/_output/fort.q000' + str(i), skiprows=8, dtype=float, names = None, header = None, sep = '\s+', usecols=np.arange(0,meqn)) #sep = '    '
            
        elif i<100:
            f = pd.read_csv(r'/home/npalex/Lid-driven-cavity/_output/fort.q00' + str(i), skiprows=8, dtype=float, names = None, header = None, sep = '\s+', usecols=np.arange(0,meqn))

        elif i<1000:
            f = pd.read_csv(r'/home/npalex/Lid-driven-cavity/_output/fort.q0' + str(i), skiprows=8, dtype=float, names = None, header = None, sep = '\s+', usecols=np.arange(0,meqn))

        else:
            f = pd.read_csv(r'/home/npalex/Lid-driven-cavity/_output/fort.q' + str(i), skiprows=8, dtype=float, names = None, header = None, sep = '\s+', usecols=np.arange(0,meqn))

        f = f.values
        f = f.T

        for j in range(0, mx):
            q[k, i, j, :] = f[k, j*mx:(j+1)*mx]
            
#------------------------------------------------
#-- plot results as an animation using matplotlib
#------------------------------------------------
q_step = 2                           # quiver plot spacing

#-- define array of cell centers
x = np.arange(xlow + dx/2, 1, dx)
y = np.arange(ylow + dy/2, 1, dy)
xgrid, ygrid = np.meshgrid(x, y)

fig = plt.figure(figsize=(6,6))
ax = plt.axes(xlim=(-0.02, 1.02), ylim=(-0.02, 1.1))            # creates axes at specifed limits, (gca() not required)

#------------------------------------------------
#-- plot vector field and pressure distribution
#------------------------------------------------
#-- define function, which is an argument for the method animation.FuncAnimation() and is called for each frame
def fplot(frame_number):
    
    ax.clear()
    
    #-- plot pressure distribution
    ax.pcolormesh(xgrid, ygrid, q[2,frame_number,:,:], cmap = 'coolwarm',vmin =-.1, vmax = .1) #--cmap = 'GnBu'
    
    #-- plot velocity distribution
    ax.quiver(x[::q_step], y[::q_step], q[0,frame_number,::q_step,::q_step], q[1,frame_number,::q_step,::q_step], scale_units = 'xy', scale = 1)

    #-- set axex limits and tick marks
    ax.set_xticks(np.arange(0, 1.2, .2))
    ax.set_yticks(np.arange(0, 1.2, .2))
    ax.set_xlim(left = 0, right = 1.)
    ax.set_ylim(bottom = 0, top = 1.)
    
    ax.set_ylabel('$y$', fontsize = 16)
    ax.set_xlabel('$x$', fontsize = 16)
    
    return()

#-- generate animation
anim = animation.FuncAnimation(fig = fig, func = fplot, frames=int(steps), interval=80, repeat=False)
plt.close()                        #-- removes residual plot at final time
#HTML(anim.to_jshtml())             #-- print animation in jupyter notebook

#-- save animation as an html file
with open("Results_vector_field_DCU_CFL_0_9.html", "w") as f:
    print(anim.to_html5_video(embed_limit=None), file=f)
    
#------------------------------------------------
#-- plot streamlines
#------------------------------------------------
fig = plt.figure(figsize = (7.5, 6), layout = "tight")
ax = plt.axes(xlim=(-0.02, 1.02), ylim=(-0.02, 1.1))            # creates axes at specifed limits, (gca() not required)
q_step = 1                           # quiver plot spacing

#-- generate colorbar
ax.axis('square')
p = ax.pcolormesh(xgrid, ygrid, q[2,0,:,:], cmap = 'coolwarm',vmin =-.06, vmax = .06)
cbar = fig.colorbar(p, fraction = 0.045, pad = 0.05, ticks = np.arange(-0.06, .08, 0.02))
cbar.set_label('$p/(\u03C1 U^2)$', size = 20)

#-- define function, which is an argument for the method animation.FuncAnimation() and is called for each frame
def fplot(frame_number):
    
    ax.clear()
    
    #-- plot pressure distribution
    ax.pcolormesh(xgrid, ygrid, q[2,frame_number,:,:], cmap = 'coolwarm',vmin =-.06, vmax = .06)
    
    #-- plot velocity distribution
    ax.streamplot(xgrid, ygrid, q[0,frame_number,::q_step,::q_step], q[1,frame_number,::q_step,::q_step], color = "black", linewidth = .3, density=.7, broken_streamlines=False)  

    #-- set axex limits and tick marks
    ax.set_xticks(np.arange(0, 1.2, .2))
    ax.set_yticks(np.arange(0, 1.2, .2))
    ax.set_xlim(left = 0, right = 1.)
    ax.set_ylim(bottom = 0, top = 1.)
    
    ax.set_ylabel('$y$', fontsize = 16)
    ax.set_xlabel('$x$', fontsize = 16)
    
    return()

#-- generate animation
anim = animation.FuncAnimation(fig = fig, func = fplot, frames=int(steps), interval=80, repeat=False)
plt.close()                        #-- removes residual plot at final time
#HTML(anim.to_jshtml())             #-- print animation in jupyter notebook

#-- save animation as an html file
with open("Results_streamlines_DCU_CFL_0_9.html", "w") as f:
    print(anim.to_html5_video(embed_limit=None), file=f)