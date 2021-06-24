import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

#---------------STREAMLIT--------------------------------------------------
st.set_page_config(layout="wide")
st.title("Potential and Flow Visualization")
col1, col2 = st.beta_columns([1,1])


#--------------SIDEBAR-----------------------------------------------------
st.sidebar.title("Parameters")
X1_para = st.sidebar.slider("Well x-cordinate", 0., 200., 99., 1.)
Y1_para = st.sidebar.slider("Well y-cordinate", 0., 200., 50., 1.)

Q_para = st.sidebar.slider("Pumping / recharge rates (Slider * 1.e-4))", -10., 10., 1., 0.1)
K_para = st.sidebar.slider("Hydraulic conductivity (Slider * 5.e-5))", 0., 10., 1., 0.1)
Por_para = st.sidebar.slider("Porosity", 0., 1., 0.25, 0.01)
Qx_para = st.sidebar.slider("Baseflow in x-direction (Slider * 1.e-10))", -100., 100., 1., 1)

#------------------VARIABLES------------------------------------------------
H = 6.             # thickness [L]
h0 = 5.5           # reference piezometric head [L] 
K = K_para * 5.e-5          # hydraulic conductivity [L/T] 
por = Por_para         # porosity []   old 0.25
Qx0 = Qx_para * 1.e-10       # baseflow in x-direction [L^2/T] was 1.e-6 before
Qy0 = 0            # baseflow in y-direction [L^2/T]
# Wells
x0 = [X1_para, 145]      # x-coordinates well position [L] [99, 145]
y0 = [Y1_para, 78]       # y-coordinates well position [L] [50, 78
Q = Q_para * np.array([1.e-4, 1.e-4])  # pumping / recharge rates [L^3/T]
R = [0.3, 0.2]      # well radius [L]
# Mesh
xmin = 0           # minimum x-position of mesh [L]
xmax = 200         # maximum x-position of mesh [L]
ymin = 0           # minimum y-position of mesh [L]
ymax = 200         # maximum y-position of mesh [L]
# Reference point position in mesh
iref = 1
jref = 1
# Graphical output options
gsurfh = 1         # piezometric head surface plot
gcontf = 10        # no. filled contour lines (=0: none)
gquiv = 0         # arrow field plot
gflowp_fit = 1     # flowpaths forward in time
gflowp_bit = 0     # no. flowpaths backward in time (=0: none)
gflowp_dot = 0     # flowpaths with dots indicating speed
gstream = 0       # streamfunction plot            10
#----------------------------------------execution-------------------------------
xvec = np.linspace(xmin,xmax,100)
yvec = np.linspace(ymin,ymax,100)
[x,y] = np.meshgrid(xvec,yvec)                     # mesh
phi = -Qx0*x - Qy0*y;                              # baseflow potential
for j in range(0, Q.ndim):                               # old: for j = 1:size(Q,2)
    r = np.sqrt((x-x0[j])*(x-x0[j])+(y-y0[j])*(y-y0[j]))   # distances to well
    phi = phi + (Q[j]/(2*np.pi))*np.log(r)              # potential
if h0 > H:
    phi0 = -phi(iref,jref) + K*H*h0 - 0.5*K*H*H 
else:
    phi0 = -phi[iref,jref] + 0.5*K*h0*h0;          # reference potential                                                 
hc = 0.5*H+(1/K/H)*(phi+phi0)                     # head confined
hu = np.sqrt((2/K)*(phi+phi0))                      # head unconfined
phicrit = phi0 + 0.5*K*H*H                        # transition confined / unconfined
confined = (phi>=phicrit)                         # confined / unconfined indicator
h = confined*hc+ ~confined*hu                    # head
psi = -Qx0*y + Qy0*x
for j in range(0, Q.ndim):
    psi = psi + (Q[j]/(np.pi+np.pi))*np.arctan2((y-y0[j]),(x-x0[j]));  # streamfunction 
#---------------------------------------display messages-------------------
#if all(all(confined)):
#    print('aquifer confined')
#else:
#    if all(all(~confined)): 
#        print('aquifer unconfined') 
#    else:
#        print('aquifer partially confined and unconfined')     
#if any(any(h<0)): 
#    print('aquifer falls partially dry') 
#    h = max(0, h)
#--------------------------------------------------------------------------
[u,v] = np.gradient(-phi)                           # discharge vector  
Hh = confined*H + ~confined*h                   # aquifer depth  
u = u/Hh/(xvec[2]-xvec[1])/por; v = v/Hh/(yvec[2]-yvec[1])/por
#--------------------------------------graphical output--------------------
if gsurfh: 
    #plt.figure()
    
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(x, y, h,
    cmap=cm.coolwarm,
    linewidth=0,
    antialiased=True)                         # surface 


    #fig.colorbar(surf, shrink=0.5, aspect=10)
if gcontf or gquiv or gflowp_fit or gflowp_bit or gflowp_dot or gstream:
    #plt.figure()
    fig2 = plt.figure()
    #fig, ax = plt.subplots(2)
    #st.pyplot(fig)
if gcontf:                                          # filled contours  
    #colormap(winter); 
    plt.contourf(x,y,h,gcontf)                  #old contourf(x,y,h,gcontf,'w')
    #colorbar
if gquiv:
    plt.quiver(x,y,u,v)                          # arrow field // quiver(x,y,u,v,'y') 
if gflowp_fit:                                      # flowpaths 
    xstart = []
    ystart = []
    for i in range(100):
        if v[1,i] > 0:
            xstart = [xstart, xvec[i]]
            ystart = [ystart, yvec[1]]
        if v[99,i] < 0:
            xstart = [xstart, xvec[i]]
            ystart = [ystart, yvec[99]]
        if u[i,1] > 0:
            xstart = [xstart, xvec[1]]
            ystart = [ystart, yvec[i]]
        if u[i,99] < 0:
            xstart = [xstart, xvec[99]]
            ystart = [ystart, yvec[i]]
    #fig, ax.streamplot(x,y,u,v)
    h = plt.streamplot(x,y,u,v)#,xstart,ystart)
    plt.streamplot(x,y,u,v)#,xstart,ystart)
    #set(h,'Color') 
if gflowp_bit:          
    for j in range(0, Q.ndim):
        if Q[j]>0:           # only for pumping wells
            xstart = x0[j] + R[j]*np.cos(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit) 
            ystart = y0[j] + R[j]*np.sin(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit)
            seed_points = np.array([xstart,ystart])
            #fig, ax.streamplot()
            h = plt.streamplot(x,y,-u,-v,start_points=seed_points.T)
            #set (h,'Color')
#if gflowp_dot:
#    [verts ~] = streamslice(x,y,u,v,gflowp_dot)
#    sc = 10/np.mean(np.mean(np.sqrt(u*u+v*v)),)
#    iverts = interpstreamspeed(x,y,u,v,verts,sc)  
#    h = streamline(iverts)
#    set (h,'Marker','.','Color','y','MarkerSize',18)
if gstream:
    plt.contour(x,y,psi,gstream)#,'k','LineWidth',1)
#plt.show()
with col1:
    st.header("3D-Plot")
    st.pyplot(fig)
with col2:
    st.header("Surfaceplot")
    st.markdown('')
    st.markdown('')
    st.markdown('')
    st.pyplot(fig2)