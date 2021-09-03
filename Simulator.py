import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
from matplotlib.ticker import StrMethodFormatter
import numpy as np
#import math
import pandas as pd
#from bokeh.plotting import figure
import base64
import pathlib



st. set_page_config(layout="wide")
##my_expander = st.sidebar.beta_expander('Navigation')
##with my_expander:
nav = st.sidebar.radio("Navigation", ["Welcome", "Wells Documentation", "Wells", " River Documentation", "River", "No Flow Documentation", "No Flow", "Special Documentation", "Special"])
if nav == "Welcome":
    st.title("Welcome to the prototype")

if nav == "Wells":
    st.title("Wells in  uniform flowfield")
    col1, col2, col3, col4, col5= st.beta_columns([20,1,20,1,20])
    #--------------SLIDER-----------------------------------------------------
    with col1:
        st.subheader("Well 1")
        X1_para = st.slider("x-cordinate [m]", 1., 199., 99., 1.)
        Y1_para = st.slider("y-cordinate [m]", 1., 199., 50., 1.)
        Q_para1 = st.number_input("Pumping / recharge rate Q1 in [m\u00B3/s] (Input * 1.e-4))", -30., 30., 1., 0.1)

    with col3:
        st.subheader("Well 2")
        X2_para = st.slider("x-cordinate [m]", 1., 199., 170., 1.)
        Y2_para = st.slider("y-cordinate [m]", 1., 199., 125., 1.)
        Q_para2 = st.number_input("Pumping / recharge rate Q2 in [m\u00B3/s] (Input * 1.e-4))", -30., 30., 0., 0.1)

    with col5:
        st.subheader("Parameters")
        #H_para = st.slider("Thickness of Aquifer [m])", 5., 10., 8., 0.1)
        #h0_para = st.slider("Reference piezometric head [m])", 5., 10., 8., 0.1)
        K_para = st.slider("Hydraulic conductivity [m/s] (Slider * 5.e-5))", 0.1, 1000., 1., 1.)
        Por_para = st.slider("Porosity", 0., 1., 0.25, 0.01)
        Qx_para = st.slider("Baseflow in x-direction [m\u00B2/s] (Slider * 1.e-10))", -10000., 10000., 0., 0.1)
        #Qy_para = st.slider("Baseflow in y-direction [m\u00B2/s] (Slider * 1.e-10))", -10., 10., 0., 0.1)
    #------------------VARIABLES------------------------------------------------
    H = 9.                                     # thickness [L]
    h0 = 9.5                                    # reference piezometric head [L] 
    K = K_para * 5.e-5                          # hydraulic conductivity [L/T] 
    por = Por_para                              # porosity []   old 0.25
    Qx0 = Qx_para * 1.e-10                      # baseflow in x-direction [L^2/T] was 1.e-6 before
    #Qy0 = Qy_para * 1.e-10                                    # baseflow in y-direction [L^2/T]
    Qy0 = 0
    # Wells
    xwell = np.array([X1_para, X2_para])        # x-coordinates well position [L] [99, 145]
    ywell = np.array([Y1_para, Y2_para])        # y-coordinates well position [L] [50, 78
    Qwell = np.array([Q_para1 * 1.e-4, Q_para2 * 1.e-4])   # pumping / recharge rates [L^3/T]
    R = [0.3, 0.2]                              # well radius [L]
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
    gquiv = 0          # arrow field plot
    gflowp_fit = 0     # flowpaths forward in time
    gflowp_bit = 0     # no. flowpaths backward in time (=0: none)
    gflowp_dot = 1     # flowpaths with dots indicating speed
    gstream = 25        # streamfunction plot            10
    #----------------------------------------execution-------------------------------
    xvec = np.linspace(xmin, xmax, 100)
    yvec = np.linspace(ymin, ymax, 100)
    [x, y] = np.meshgrid(xvec, yvec)                        # mesh
    phi = -Qx0 * x - Qy0 * y                                # baseflow potential
    psi = -Qx0 * y + Qy0 * x
    for i in range(0, xwell.size):                          # old version was: for i = 1:size(xwell,2)
        #r = np.sqrt((x - xwell[i]) * (x - xwell[i]) + (y - ywell[i]) * (y - ywell[i]))
        r = np.sqrt((x - xwell[i])**2 + (y - ywell[i])**2)
        phi = phi + (Qwell[i] / (2 * np.pi)) * np.log(r)    # potential
        psi = psi + (Qwell[i]/ (2 * np.pi)) * np.arctan2((y - ywell[i]), (x - xwell[i]))
    if h0 > H:
        phi0 = -phi[iref, jref] + K * H * h0 - 0.5 * K * H * H 
    else:
        phi0 = -phi[iref, jref] + 0.5 * K * h0 * h0          # reference potential                                                 
    hc = 0.5 * H + (1 / K / H) * (phi + phi0)                     # head confined
    hu = np.sqrt((2 / K) * (phi + phi0))                      # head unconfined
    phicrit = phi0 + 0.5 * K * H * H                        # transition confined / unconfined
    confined = (phi >= phicrit)                         # confined / unconfined indicator
    h = confined * hc+ ~confined * hu                    # head

    [u,v] = np.gradient(-phi)                           # discharge vector  
    Hh = confined * H + ~confined * h                   # aquifer depth  
    u = u / Hh / (xvec[2] - xvec[1]) / por
    v = v / Hh / (yvec[2] - yvec[1]) / por
    #--------------------------------------graphical output--------------------
    if gsurfh: 
        #plt.figure()
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(x, y, h,
                            cmap=cm.coolwarm,
                            linewidth=0,
                            antialiased=True)                         # surface 
        plt.gca().zaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        #ax.set_title('3-Plot')
        # set z (h)-axis invisible for better visibilty if colorbar is active
        #if fig.colorbar:
        #    ax.set_zticks([])
        #else:
        ax.set_zlabel('drawdown [m]')
        fig.colorbar(surf, shrink=.8, ax=[ax], location = "left") # ax=[ax], location='left' for left side

        #fig.colorbar(surf, shrink=0.5, aspect=10)
    if gcontf or gquiv or gflowp_fit or gflowp_bit or gflowp_dot or gstream:
        #fig2 = plt.figure()
        fig2, ax = plt.subplots()
        contour = plt.contour(x, y, h,
            gcontf,
            cmap = cm.Blues)
            #colors=['#808080', '#A0A0A0', '#C0C0C0'], extend='both')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        # fig2, ax = plt.subplots()
        contour2 = plt.contour(x ,y, psi,
            gstream,
            #cmap = cm.binary
            colors=['#808080', '#808080', '#808080'], extend='both')
        contour.cmap.set_over('red')
        contour.cmap.set_under('blue')
        contour.changed()
        #plt.clabel(contour, inline=1, fontsize=10)
        labels = ['Streamline', 'Potentialline']
        contour2.collections[8].set_label(labels[0])
        contour.collections[7].set_label(labels[1])

        plt.legend(loc='upper left')
    # if gstream:
    #     plt.contour(x ,y, psi,
    #             gstream,
    #             cmmap = cm.Greys)
        #set(psi,'LineColor','none')
    #fig.colorbar(contour, shrink=.8, ax=[ax], location = "left")    
    # if gcontf:                                          # filled contours  
    #     #colormap(winter); 
    #     plt.contour(x, y, h,
    #                 gcontf)         
    #fig.colorbar(contour, shrink=.8, ax=[ax], location = "left")                       
        #colorbar

    if gquiv:
        plt.quiver(x,y,v,u)                          # arrow field // quiver(x,y,u,v,'y') 
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
        h = plt.streamplot(x,y,u,v,)#,xstart,ystart)
        plt.streamplot(x,y,u,v,color='b')#,xstart,ystart)
        #set(h,'Color') 
    if gflowp_bit:          
        for j in range(0, Qwell.size):
            if Qwell[j]>0:           # only for pumping wells
                xstart = xwell[j] + R[j]*np.cos(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit) 
                ystart = ywell[j] + R[j]*np.sin(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit)
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
    #plt.show()
    with col3:
        st.header("Drawdown-Plot")
        st.pyplot(fig)
    with col1:
        st.header("Flowfield-Plot")
        st.markdown('')
        st.pyplot(fig2)

    dfh = pd.DataFrame(data = h)
    dfh_rounded = dfh.round(decimals = 3)

#------------------DOWNLOAD CSV FILE----------------------------------------------------------------------#
    csv = dfh_rounded.to_csv(sep="\t", index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}">CSV File for head</a>'
    with col5:
        st.markdown('')
        st.markdown('')
        st.markdown('')
        st.markdown('')
        st.markdown('')
        st.markdown('')
        '''
        **Results:**  
        (right-click and save as name.csv)  
        '''
        st.markdown(href, unsafe_allow_html=True)
    
    st.markdown('')
    st.markdown('')
    '''
    **Sourcefile:** https://github.com/TReuschling/Flow (open-source CC BY 4.0)  
    **Disclaimer:** Authors of code are not responsible for obtained results.
    '''
#-----------------------------RIVER------------------------------------------#
#-----------------------------RIVER------------------------------------------#
#-----------------------------RIVER------------------------------------------#
if nav == "River":
    st.title("Potential and Flow Visualization")
    col1, col2 = st.beta_columns([1,1])


    #--------------INPUT PARAMETERS-----------------------------------------------------
    with col1:
        X1_para = st.number_input("Well1 x-cordinate [m]", 1., 199., 100., 1.)
        Y1_para = st.slider("Well1 y-cordinate [m]", 1., 199., 100., 1.)
        Q_para1 = st.slider("Pumping / recharge rate1 in [m\u00B3/s] (Slider * 1.e-4))", -100., 100., 8., 0.1)
    with col2:
        X2_para = st.slider("Well2 x-cordinate [m]", 1., 199., 100., 1.)
        Y2_para = st.slider("Well2 y-cordinate [m]", 1., 199., 150., 1.)
        Q_para2 = st.slider("Pumping / recharge rate2 in [m\u00B3/s] (Slider * 1.e-4))", -10., 10., 0., 0.1)

    Qx_para = st.slider("Baseflow in x-direction [m\u00B2/s] (Slider * 1.e-10))", 0, 100000, 32000, 10)
    K_para = st.slider("Hydraulic conductivity [m/s] (Slider * 5.e-5))", 0.1, 10., 1., 0.1)
    Por_para = st.slider("Porosity", 0., 1., 0.25, 0.01)


    #------------------VARIABLES------------------------------------------------
    H = 10.                                      # thickness [L]
    h0 = 9.5                                    # reference piezometric head [L] 
    K = K_para * 5.e-5                          # hydraulic conductivity [L/T] 
    por = Por_para                              # porosity []   old 0.25
    Qx0 = Qx_para * 1.e-10                      # baseflow in x-direction [L^2/T] was 1.e-6 before
    #Qx0 = 0
    Qy0 = 0                                     # baseflow in y-direction [L^2/T]
    # Mesh
    xmin = 0           # minimum x-position of mesh [L]
    xmax = 200         # maximum x-position of mesh [L]
    ymin = 0           # minimum y-position of mesh [L]
    ymax = 200         # maximum y-position of mesh [L]
    # Wells
    xwell = np.array([X1_para, 2 * xmax - X1_para, X2_para, 2 * xmax + X2_para])        # x-coordinates well position [L] [99, 145]
    ywell = np.array([Y1_para, Y1_para, Y2_para, Y2_para])        # y-coordinates well position [L] [50, 78
    Qwell = np.array([Q_para1 * 1.e-4, -Q_para1 * 1.e-4, Q_para2 * 1.e-4, -Q_para2 * 1.e-4])   # pumping / recharge rates [L^3/T]
    R = [0.2, 0.2, 0.2, 0.2]                              # well radius [L]
    # Second mesh for grap
    xmax2 = ymax * 2         # maximum x-position of mesh [L]
    ymax2 = xmax * 2        # maximum y-position of mesh [L]
    # Reference point position in mesh
    iref = 1
    jref = 1
    # Graphical output options
    gsurfh = 1         # piezometric head surface plot
    gcontf = 10       # no. filled contour lines (=0: none)
    gquiv = 1          # arrow field plot
    gflowp_fit = 0     # flowpaths forward in time
    gflowp_bit = 0     # no. flowpaths backward in time (=0: none)
    gflowp_dot = 1     # flowpaths with dots indicating speed
    gstream = 25        # streamfunction plot            10
    #----------------------------------------execution-------------------------------
    xvec = np.linspace(xmin, xmax, 100)
    yvec = np.linspace(ymin, ymax, 100)
    [x, y] = np.meshgrid(xvec, yvec)                        # mesh

    #xvec2 = np.linspace(xmin, xmax2, 100)
    #yvec2 = np.linspace(ymin, ymax2, 100)
    #[x2, y2] = np.meshgrid(xvec2, yvec2)   
    #x2[x2 > 100] = np.nan                                     # mesh

    phi = -Qx0 * x - Qy0 * y                                # baseflow potential
    psi = -Qx0 * y + Qy0 * x

    for i in range(0, xwell.size):                       
        r = np.sqrt((x - xwell[i])**2 + (y - ywell[i])**2)
        phi = phi + (Qwell[i] / (2 * np.pi)) * np.log(r)    
        psi = psi + (Qwell[i]/ (2 * np.pi)) * np.arctan2((y - ywell[i]), (x - xwell[i]))

    if h0 > H:
        phi0 = -phi(iref, jref) + K * H * h0 - 0.5 * K * H * H 
    else:
        phi0 = -phi[iref, jref] + 0.5 * K * h0 * h0          # reference potential                                                 
    hc = 0.5 * H + (1 / K / H) * (phi + phi0)                     # head confined
    hu = np.sqrt((2 / K) * (phi + phi0))                      # head unconfined
    phicrit = phi0 + 0.5 * K * H * H                        # transition confined / unconfined
    confined = (phi >= phicrit)                         # confined / unconfined indicator
    h = confined * hc+ ~confined * hu                    # head
    #--------------------------------------------------------------------------
    [u,v] = np.gradient(-phi)                           # discharge vector  
    Hh = confined * H + ~confined * h                   # aquifer depth  
    u = u / Hh / (xvec[2] - xvec[1]) / por
    v = v / Hh / (yvec[2] - yvec[1]) / por
    #--------------------------------------graphical output--------------------
    if gsurfh: 
        #plt.figure()
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(x, y, h,
                            cmap=cm.coolwarm,
                            linewidth=0.1,
                            antialiased=True)
        #fig.colorbar(surf, ax=ax, shrink=.8)
        #ax.margins(x= 0, y=0)
        #plt.xlim(0,200)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        # set z (h)-axis invisible for better visibilty if colorbar is active
        #if fig.colorbar:
        #    ax.set_zticks([])
        #else:
        ax.set_zlabel('drawdown [m]')
        fig.colorbar(surf, shrink=.8, ax=[ax], location = "left") # ax=[ax], location='left' for left side


        #fig.colorbar(surf, shrink=0.5, aspect=10)
    if gcontf or gquiv or gflowp_fit or gflowp_bit or gflowp_dot or gstream:
        #plt.figure()
        fig2 = plt.figure()
        #fig, ax = plt.subplots(2)
        #st.pyplot(fig)
    if gcontf:                                          # filled contours  
        fig3, ax = plt.subplots()
        contour = plt.contour(x, y, h,
                    gcontf)
        ax.set_title('Contour Plot')
        #ax.margins(x= -0.25, y=0)
        #plt.xlim(0,200)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        # Create a Rectangle patch
        rect = Rectangle((198, 0),2,200,linewidth=1,edgecolor='b',facecolor='b', zorder=2) # zorder makes is intransparent
        # Add the patch to the Axes
        ax.add_patch(rect)
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
        h = plt.streamplot(x,y,u,v,)#,xstart,ystart)
        plt.streamplot(x,y,u,v,color='b')#,xstart,ystart)
        #set(h,'Color') 
    if gflowp_bit:          
        for j in range(0, Qwell.size):
            if Qwell[j]>0:           # only for pumping wells
                xstart = xwell[j] + R[j]*np.cos(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit) 
                ystart = ywell[j] + R[j]*np.sin(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit)
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

    with col1:
        st.header("3D-Plot")
        st.pyplot(fig)
    with col2:
        st.header("Surfaceplot")
        st.pyplot(fig3)

    # df = pd.DataFrame(data = phi)
    # df.round(decimals = 3)
    # df.to_csv(r'/Users/tassiloreuschling/Uni/Bachelorarbeit/Data_aus_python/phi.csv', index = False)

    # df1 = pd.DataFrame(data = psi)
    # df1.round(decimals = 3)
    # df1.to_csv(r'/Users/tassiloreuschling/Uni/Bachelorarbeit/Data_aus_python/psi.csv', index = False)

    dfh = pd.DataFrame(data = h)
    dfh_rounded = dfh.round(decimals = 3)

#------------------DOWNLOAD CSV FILE----------------------------------------------------------------------#
    csv = dfh_rounded.to_csv(sep="\t", index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}">Download CSV File</a> (right-click and save as &lt;some_name&gt;.csv)'
    st.markdown(href, unsafe_allow_html=True)


    # # HACK This only works when we've installed streamlit with pipenv, so the
    # # permissions during install are the same as the running process
    # STREAMLIT_STATIC_PATH = pathlib.Path(st.__path__[0]) / 'static'
    # # We create a downloads directory within the streamlit static asset directory
    # # and we write output files to it
    # DOWNLOADS_PATH = (STREAMLIT_STATIC_PATH / "downloads")
    # if not DOWNLOADS_PATH.is_dir():
    #     DOWNLOADS_PATH.mkdir()

    # def main():
    #     st.markdown("Download from [downloads/mydata.csv](downloads/mydata.csv)")
    #     mydataframe = pd.DataFrame.from_dict({'col_1': [3, 2, 1, 0], 'col_2': ['a', 'b', 'c', 'd']})
    #     mydataframe.to_csv(str(DOWNLOADS_PATH / "mydata.csv"), index=False)

    # if __name__ == "__main__":
    #     main()

    ratio = Qwell[0] / ( Qx0 * np.pi * ((2 * xmax - X1_para)- X1_para)/2)

    with col2:
        st.write(ratio)
    print(ratio)


#-----------------------------NO FLOW----------------------------------------#
#-----------------------------NO FLOW----------------------------------------#
#-----------------------------NO FLOW----------------------------------------#

if nav == "No Flow":
    col1, col2 = st.beta_columns([1,1])


    #--------------SIDEBAR-----------------------------------------------------
    with col1:
        X1_para = st.slider("Well1 x-cordinate [m]", 1., 200., 100., 1.)
        Y1_para = st.slider("Well1 y-cordinate [m]", 1., 200., 100., 1.)
        Q_para1 = st.slider("Pumping / recharge rate1 in [m\u00B3/s] (Slider * 1.e-4))", -100., 100., 8., 0.1)
    with col2:
        X2_para = st.slider("Well2 x-cordinate [m]", 1., 200., 100., 1.)
        Y2_para = st.slider("Well2 y-cordinate [m]", 1., 200., 150., 1.)
        Q_para2 = st.slider("Pumping / recharge rate2 in [m\u00B3/s] (Slider * 1.e-4))", -10., 10., 0., 0.1)

    K_para = st.slider("Hydraulic conductivity [m/s] (Slider * 5.e-5))", 0.1, 10., 1., 0.1)
    Por_para = st.slider("Porosity", 0., 1., 0.25, 0.01)


    #------------------VARIABLES------------------------------------------------
    H = 10.                                      # thickness [L]
    h0 = 9.5                                    # reference piezometric head [L] 
    K = K_para * 5.e-5                          # hydraulic conductivity [L/T] 
    por = Por_para                              # porosity []   old 0.25
    Qx0 = 0                    # baseflow in x-direction [L^2/T] was 1.e-6 before
    Qy0 = 0                                     # baseflow in y-direction [L^2/T]
    # Mesh
    xmin = 0           # minimum x-position of mesh [L]
    xmax = 200         # maximum x-position of mesh [L]
    ymin = 0           # minimum y-position of mesh [L]
    ymax = 200         # maximum y-position of mesh [L]
    # Wells
    xwell = np.array([X1_para, 2 * xmax - X1_para, X2_para, 2 * xmax + X2_para])        # x-coordinates well position [L] [99, 145]
    ywell = np.array([Y1_para, Y1_para, Y2_para, Y2_para])        # y-coordinates well position [L] [50, 78
    Qwell = np.array([Q_para1 * 1.e-4, Q_para1 * 1.e-4, Q_para2 * 1.e-4, Q_para2 * 1.e-4])   # pumping / recharge rates [L^3/T]
    R = [0.2, 0.2, 0.2, 0.2]                              # well radius [L]
    # Second mesh for grap
    xmax2 = ymax * 2         # maximum x-position of mesh [L]
    ymax2 = xmax * 2        # maximum y-position of mesh [L]
    # Reference point position in mesh
    iref = 1
    jref = 1
    # Graphical output options
    gsurfh = 1         # piezometric head surface plot
    gcontf = 10       # no. filled contour lines (=0: none)
    gquiv = 1          # arrow field plot
    gflowp_fit = 0     # flowpaths forward in time
    gflowp_bit = 0     # no. flowpaths backward in time (=0: none)
    gflowp_dot = 1     # flowpaths with dots indicating speed
    gstream = 25        # streamfunction plot            10
    #----------------------------------------execution-------------------------------
    xvec = np.linspace(xmin, xmax, 100)
    yvec = np.linspace(ymin, ymax, 100)
    [x, y] = np.meshgrid(xvec, yvec)                        # mesh

    #xvec2 = np.linspace(xmin, xmax2, 100)
    #yvec2 = np.linspace(ymin, ymax2, 100)
    #[x2, y2] = np.meshgrid(xvec2, yvec2)   
    #x2[x2 > 100] = np.nan                                     # mesh

    phi = -Qx0 * x - Qy0 * y                                # baseflow potential
    psi = -Qx0 * y + Qy0 * x

    for i in range(0, xwell.size):                       
        r = np.sqrt((x - xwell[i])**2 + (y - ywell[i])**2)
        phi = phi + (Qwell[i] / (2 * np.pi)) * np.log(r)    
        psi = psi + (Qwell[i]/ (2 * np.pi)) * np.arctan2((y - ywell[i]), (x - xwell[i]))

    if h0 > H:
        phi0 = -phi(iref, jref) + K * H * h0 - 0.5 * K * H * H 
    else:
        phi0 = -phi[iref, jref] + 0.5 * K * h0 * h0          # reference potential                                                 
    hc = 0.5 * H + (1 / K / H) * (phi + phi0)                     # head confined
    hu = np.sqrt((2 / K) * (phi + phi0))                      # head unconfined
    phicrit = phi0 + 0.5 * K * H * H                        # transition confined / unconfined
    confined = (phi >= phicrit)                         # confined / unconfined indicator
    h = confined * hc+ ~confined * hu                    # head
    #--------------------------------------------------------------------------
    [u,v] = np.gradient(-phi)                           # discharge vector  
    Hh = confined * H + ~confined * h                   # aquifer depth  
    u = u / Hh / (xvec[2] - xvec[1]) / por
    v = v / Hh / (yvec[2] - yvec[1]) / por
    #--------------------------------------graphical output--------------------
    if gsurfh: 
        #plt.figure()
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(x, y, h,
                            cmap=cm.coolwarm,
                            linewidth=0.1,
                            antialiased=True)
        #fig.colorbar(surf, ax=ax, shrink=.8)
        #ax.margins(x= 0, y=0)
        #plt.xlim(0,200)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        # set z (h)-axis invisible for better visibilty if colorbar is active
        #if fig.colorbar:
        #    ax.set_zticks([])
        #else:
        ax.set_zlabel('drawdown [m]')
        fig.colorbar(surf, shrink=.8, ax=[ax], location = "left") # ax=[ax], location='left' for left side


        #fig.colorbar(surf, shrink=0.5, aspect=10)
    if gcontf or gquiv or gflowp_fit or gflowp_bit or gflowp_dot or gstream:
        #plt.figure()
        fig2 = plt.figure()
        #fig, ax = plt.subplots(2)
        #st.pyplot(fig)
    if gcontf:                                          # filled contours  
        fig3, ax = plt.subplots()
        contour = plt.contour(x, y, h,
                    gcontf)
        ax.set_title('Contour Plot')
        #ax.margins(x= -0.25, y=0)
        #plt.xlim(0,200)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        # Create a Rectangle patch
        rect = Rectangle((198, 0),2,200,linewidth=1,edgecolor='k',facecolor='k', zorder=2) # zorder makes is intransparent
        # Add the patch to the Axes
        ax.add_patch(rect)
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
        h = plt.streamplot(x,y,u,v,)#,xstart,ystart)
        plt.streamplot(x,y,u,v,color='b')#,xstart,ystart)
        #set(h,'Color') 
    if gflowp_bit:          
        for j in range(0, Qwell.size):
            if Qwell[j]>0:           # only for pumping wells
                xstart = xwell[j] + R[j]*np.cos(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit) 
                ystart = ywell[j] + R[j]*np.sin(2*np.pi*np.array([1,1,gflowp_bit])/gflowp_bit)
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

    with col1:
        st.header("3D-Plot")
        st.pyplot(fig)
    with col2:
        st.header("Surfaceplot")
        st.pyplot(fig3)

    # df = pd.DataFrame(data = phi)
    # df.round(decimals = 3)
    # df.to_csv(r'/Users/tassiloreuschling/Uni/Bachelorarbeit/Data_aus_python/phi.csv', index = False)

    # df1 = pd.DataFrame(data = psi)
    # df1.round(decimals = 3)
    # df1.to_csv(r'/Users/tassiloreuschling/Uni/Bachelorarbeit/Data_aus_python/psi.csv', index = False)

    # dfh = pd.DataFrame(data = h)
    # dfh_rounded = dfh.round(decimals = 3)
    # dfh_rounded.to_csv(r'/Users/tassiloreuschling/Uni/Bachelorarbeit/Data_aus_python/h.csv', index = False)

    # ratio = Qwell[0] / ( Qx0 * np.pi * ((2 * xmax - X1_para)- X1_para)/2)

    # with col2:
    #     st.write(ratio)
