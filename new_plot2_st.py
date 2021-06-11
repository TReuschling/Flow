import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import math
import pandas as pd

#st. set_page_config(layout="wide")
col1, col2 = st.beta_columns([1,2])
#with col1:
xy_Size = st.sidebar.slider('Size of Model', 1, 100, 50, 1)

xmin = -xy_Size
xmax = xy_Size
ymin = -xy_Size
ymax = xy_Size

#with col1:
x_loca = st.sidebar.slider('x-location of source/sink', xmin, xmax, 0, 1)
y_loca = st.sidebar.slider('y-location of source/sink', ymin, ymax, 0, 1)

#source
#x0 = 0
#y0 = .905
x0 = x_loca
y0 = y_loca
#Q = 50
#with col1:
Q = st.sidebar.slider('Q', -50, 50, 50, 1)

#Qx0 = 0.1
#Qy0 = 0
#with col1:
Qx0 = st.sidebar.slider('Qx-flow', -1.0, 1.0, 0.1, 0.0001)
Qy0 = st.sidebar.slider('Qy-flow', -1.0, 1.0, 0.0, 0.0001)



# mesh generation
xvec = np.linspace(xmin,xmax,100)
yvec = np.linspace(ymin,ymax,100)
[x,y] = np.meshgrid(xvec,yvec)

# processing
r = np.sqrt((x-x0)**2+(y-y0)**2)
phi = -Qx0*x - Qy0*y + (Q/(2*np.pi))*np.log(r); # potential
#post-processing

#surf (x,y,phi)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(x, y, phi, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=10)
#plt.show()
#with col2:
st.pyplot(fig)
