import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import drag_functions as dfunc

def plot_data(Lambda,ax):
    
    if int(Lambda)==5:
        df5 = pd.read_csv("C:\\Users\\lucca\\Desktop\\Friction_Biofilms\\Analysis\\data\\Inner_fraction_Lambda5", sep='\t')
        time5 = np.arange(0,len(df5["avg"])*0.1 - 0.05,0.1)
        ax.errorbar(time5,df5["avg"],yerr=df5["stdev"],c=dfunc.colour_Lambda(5.0),label=r"sim $\Lambda = 5$",linewidth=2)
    if int(Lambda)==10:
        df10 = pd.read_csv("C:\\Users\\lucca\\Desktop\\Friction_Biofilms\\Analysis\\data\\Inner_fraction_Lambda10", sep='\t')
        time10 = np.arange(0,len(df10["avg"])*0.1 - 0.05,0.1)
        ax.errorbar(time10,df10["avg"],yerr=df10["stdev"],c=dfunc.colour_Lambda(10.0),label=r"sim $\Lambda = 10$",linewidth=2)
        

def R_g(t,R0=10,tau=1.8):
    return R0*(2**(t/tau))

def s_green_evolution(time,R0=10,tau=1.8):
    sgreen = 2*np.pi*R_g(time,R0=R0,tau=tau)
    return sgreen

def dsb_dt_sum(t,sblue,mu=5,av_l=5,R0=10,phi0=np.pi/2,tau=1.8):
    n = sblue/(2*av_l)
    dsdt = 1
    for i in range(1,int(n)):
        if i*av_l/(R_g(t,R0=R0,tau=tau)) < phi0:
            dsdt += np.cos(i*av_l/(R_g(t,R0=R0,tau=tau)))
    return dsdt*2*mu

def dsb_dt_int(t,sblue,mu=5,av_l=5,R0=10,tau=1.8):
    Rg = R_g(t,R0=R0,tau=tau)
    imax = sblue/(2*av_l)
    return 2*mu*(Rg/av_l)*np.sin(imax*av_l/Rg)

def s_blue_evolution(time,s0,dt,mu=5,av_l=5,R0=10,tau=1.8):
    n_timesteps = len(time)
    sblue = np.zeros(n_timesteps)
    sblue[0] = s0
    for t in range(1,n_timesteps):
        dsdt = dsb_dt_int(time[t],sblue[t-1],mu=mu,av_l=av_l,R0=R0,tau=tau)
        sblue[t] = sblue[t-1] +dsdt*dt
    return sblue

def evolve_time(time,frac_s0=0.05,dt=0.0001,mu=5,av_l=5,R0=10,tau=1.8):
    sgreen = s_green_evolution(time,R0=R0,tau=tau)
    s0 = sgreen[0]*frac_s0
    sblue = s_blue_evolution(time,s0,dt,mu=mu,av_l=av_l,R0=R0,tau=tau)
    return sgreen,sblue



def define_axis(axes):
    axes[0].set_ylabel("$s_{blue}/s_{green}$",fontsize = 20)
    axes[0].set_xlabel(r'Time after collision, $t - t_{collision}$ (h)',fontsize = 20)
    axes[0].legend()


    axes[1].set_ylabel("$s_{blue}$, $s_{green}$",fontsize = 20)
    axes[1].set_xlabel(r'Time after collision, $t - t_{collision}$ (h)',fontsize = 20)


def plot_average_inner_fraction(t0,tend=6,Lambda=10,ax=None,line=True):
    dt = 0.0001
    time = np.arange(0,tend,dt)

    mu=5
    av_l=5
    tau=1.8
    s0 = 0.1

    fraction = np.zeros((len(time),len(t0)))
    for i in range(0,len(t0)):
        r0 = np.sqrt(2)*2**(t0[i]/tau)
        sgreen,sblue = evolve_time(time,frac_s0=s0,dt=dt,mu=mu,av_l=av_l,R0=r0,tau=tau)
        fraction[:,i] = sblue/sgreen
    
    av_fraction = np.mean(fraction,axis=1)
    err_fraction = np.std(fraction,axis=1)/np.sqrt(len(t0))

    if ax ==None:
        plt.plot(time,av_fraction,"-",c="k",label="analytical model")
        plt.fill_between(time, av_fraction-err_fraction, av_fraction+err_fraction, facecolor="k",alpha=0.2)
    else:
        if line:
            ax.plot(time,av_fraction,"-",c="k",label="analytical model")
            ax.fill_between(time, av_fraction-err_fraction, av_fraction+err_fraction, facecolor="k",alpha=0.2)
        plot_data(Lambda,ax)