from math import dist

import matplotlib.pyplot as plt
import numpy as np
import drag_functions as dfunc
import utilities as ut
import pandas as pd
import os
import re
import glob


#plots quantity vs time
def plot_count(data_dir,channels=False, width=60):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir, channels, width)
    Lambdas = sorted(Lambdas)

    if Lambdas[1] == 1.0:
        plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(Lambdas[0]))
        return [Lambdas[0]]
        

    if len(Lambdas) != 1.0:
        plt.plot(time_steps, Lambda2_counts, 'o', color=dfunc.colour_Lambda(Lambdas[1]))
        return [Lambdas[1]]
    plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(Lambdas[0]))
    return Lambdas

def plot_count_tot(data_dir,channels=False, width=60):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir, channels, width)
    Lambdas = sorted(Lambdas)
    for Lambda in Lambdas:
        if Lambda != 1:
            plt.plot(time_steps, np.asarray(Lambda1_counts)+np.asarray(Lambda2_counts), color=dfunc.colour_Lambda(Lambda))
            return [Lambda]


def plot_rg_linear(data_dir):
    time_steps, Rg_tot, Lambdas, Rg1, Rg2 = dfunc.RgLambda_time(data_dir)

    if len(Lambdas) != 1:
        params,err = dfunc.RgLambda_time_est(Rg1)
        print(f"Behaviour for Lambda={Lambdas[0]}: log2(R_g) = ({round(params[0],3)}+-{round(err[0],3)}) + t/({round(params[1],3)}+-{round(err[1],3)})")
        plt.plot(time_steps, Rg1, '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[0]))
        plt.plot(time_steps, dfunc.Rg_Lambda_growth(np.array(time_steps), *params), '-', color=dfunc.colour_Lambda(Lambdas[0]))

        params,err = dfunc.RgLambda_time_est(Rg2)
        print(f"Behaviour for Lambda={Lambdas[1]}: log2(R_g) = ({round(params[0],3)}+-{round(err[0],3)}) + t/({round(params[1],3)}+-{round(err[1],3)})")
        plt.plot(time_steps, Rg2, '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[1]))
        plt.plot(time_steps, dfunc.Rg_Lambda_growth(np.array(time_steps), *params), '-', color=dfunc.colour_Lambda(Lambdas[1]))    

        params_tot,err_tot = dfunc.RgLambda_time_est(Rg_tot,total=True)
        print(f"Behaviour for total colony: log2(R_g) =({round(params_tot[0],3)}+-{round(err_tot[0],3)}) + t/({round(params_tot[1],3)}+-{round(err_tot[1],3)})")
        plt.plot(time_steps, Rg_tot, '.', markersize=6, color="k")#, label=" Counts $\Lambda = 1$")
        plt.plot(time_steps, dfunc.Rg_Lambda_growth(np.array(time_steps), *params_tot), '-', color="k")
 
    else:
        params_tot,err_tot = dfunc.RgLambda_time_est(Rg_tot)
        print(f"Behaviour for $\Lambda={Lambdas[0]}$: $log2(R_g) =({round(params_tot[0],3)}+-{round(err_tot[0],3)}) + t/({round(params_tot[1],3)}+-{round(err_tot[1],3)})")
        plt.plot(time_steps, np.log2(Rg_tot), '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[0]))#, label=" Counts $\Lambda = 1$")
        plt.plot(time_steps, dfunc.Rg_Lambda_growth(np.array(time_steps), *params_tot), '-', color=dfunc.colour_Lambda(Lambdas[0]))
    return Lambdas

def plot_GR_linear(data_dir):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)
    
    params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(Lambda1_counts, time_step=0.1)
    print(f"Behaviour for $\Lambda={Lambdas[0]}$: $log_2 N(t) =({round(params_Lambda1[1],3)}+-{round(err_Lambda1[1],3)}) + t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
    plt.plot(time_steps, np.log2(Lambda1_counts), '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[0]))#, label=" Counts $\Lambda = 1$")
    plt.plot(time_steps, dfunc.doubling_linear_growth(np.array(time_steps), *params_Lambda1), '-', color=dfunc.colour_Lambda(1.0))#, label="Fit: $log_2(N)="+str(round(params_Lambda1[0], 2)) +" + t/"+str(round(params_Lambda1[1], 2)) +"$")


    if len(Lambdas) != 1:
        params_Lambda2, err_Lambda2 = dfunc.estimate_growth_rate(Lambda2_counts, time_step=0.1)
        print(f"Behaviour for $\Lambda={Lambdas[1]}$: $log_2 N(t) =({round(params_Lambda2[1],3)}+-{round(err_Lambda2[1],3)}) + t/({round(params_Lambda2[0],3)}+-{round(err_Lambda2[0],3)})")
        plt.plot(time_steps, np.log2(Lambda2_counts), '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[1]))#, label="$Counts \Lambda =" + str(Lambdas[Lambdas != 1.0][0])+"$")
        plt.plot(time_steps, dfunc.doubling_linear_growth(np.array(time_steps), *params_Lambda2), '-', color=dfunc.colour_Lambda(Lambdas[1]))#, label="Fit: $log_2(N)="+str(round(params_Lambdanot1[0], 2)) +" + t/"+str(round(params_Lambdanot1[1], 2)) +"$")

    plt.ylabel("Count $log_2 [N(t)]$")
    return Lambdas

def plot_GR_exp(data_dir):
    time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
    
    if len(Lambda1_counts)>0:
        params_Lambda1, err_Lambda1 = dfunc.estimate_exp_growth_rate(Lambda1_counts, time_step=0.1)
        print(f"Behaviour for $\Lambda = 1$: $N(t) = 2^t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
        if len(Lambda1_counts) == len(time_steps):
            plt.plot(time_steps, Lambda1_counts, '.', markersize=6, color=dfunc.colour_Lambda(1.0))#, label=" Counts $\Lambda = 1$")
            plt.plot(time_steps, dfunc.doubling_exp_growth(np.array(time_steps), params_Lambda1[0]), '-', color=dfunc.colour_Lambda(1.0))#, label="Fit: $log_2(N)="+str(round(params_Lambda1[0], 2)) +" + t/"+str(round(params_Lambda1[1], 2)) +"$")
            
    if len(Lambdanot1_counts)>0:
        params_Lambdanot1, err_Lambdanot1 = dfunc.estimate_exp_growth_rate(Lambdanot1_counts, time_step=0.1)
        print(f"Behaviour for $\Lambda={Lambdas[Lambdas != 1.0][0]}$: $N(t) = 2^t/({round(params_Lambdanot1[0],3)}+-{round(err_Lambdanot1[0],3)})")
        if len(Lambdanot1_counts) == len(time_steps):
            plt.plot(time_steps, Lambdanot1_counts, '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]))#, label="$Counts \Lambda =" + str(Lambdas[Lambdas != 1.0][0])+"$")
            plt.plot(time_steps, dfunc.doubling_exp_growth(np.array(time_steps), params_Lambdanot1[0]), '-', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]))#, label="Fit: $log_2(N)="+str(round(params_Lambdanot1[0], 2)) +" + t/"+str(round(params_Lambdanot1[1], 2)) +"$")
    plt.ylabel("Count $N(t)$")
    return Lambdas

def plot_shape_asphericity_time(data_dir):

    """
    2D Shape asphericity in this case is defined:
    For eig1**2, eig2**2 being the 2 eigenvalues of the gyration tensor,
    diff =  (eig1-eig2)**2 /(eig1+eig2)**2 
    
    definition described in https://pubs.aip.org/aip/jcp/article/160/1/014906/2932465/Onset-of-glassiness-in-two-dimensional-ring
    """
    time_steps, Lambdas, asph1, asph2 = dfunc.calc_asphericity(data_dir)

    if len(Lambdas) == 1:
        plt.plot(time_steps, asph1, 'o', color=dfunc.colour_Lambda(Lambdas[0]))
        return Lambdas

    else: #If only one Lambda, skip the rest
        plt.plot(time_steps, asph1, 'o', color=dfunc.colour_Lambda(Lambdas[0]))
        plt.plot(time_steps, asph2, 'o', color=dfunc.colour_Lambda(Lambdas[1]))
        return Lambdas

def plot_delta_asph_time(data_dir):
    time_steps, Lambdas, asph1, asph2 = dfunc.calc_asphericity(data_dir)
    dasph = np.asarray(asph1) -np.asarray(asph2)

    if len(Lambdas) == 1:
        plt.plot(time_steps, asph1, color=dfunc.colour_Lambda(Lambdas[0]))
        print("OBS: ONLY ONE $\LAMBDA$, CANNOT CALCULATE THE DIFFERENCE")
        return Lambdas

    else: # color is the higher Lambda.
        plt.plot(time_steps, dasph, color=dfunc.colour_Lambda(Lambdas[1]))
        return Lambdas

def plot_stress_time(data_dir):
    files = get_file_paths(data_dir)

    time_steps = []

    Lambda1_stress_perp = []
    Lambda1_stress_par = []
    Lambda1_stress_shear1 = []
    Lambda1_stress_shear2 = []
    err_Lamba1 =np.zeros((len(files),4))

    Lambdanot1_stress_shear1 = []
    Lambdanot1_stress_shear2 = []
    Lambdanot1_stress_par = []
    Lambdanot1_stress_perp = []
    err_Lambanot1 =np.zeros((len(files),4))

    i = 0
    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        df = pd.read_csv(file_path, sep="\t")
        time_steps.append(time_step * 0.1)
        dfunc.find_stress(df)
        Lambdas = dfunc.find_Lambdas(df)
        for Lambda in Lambdas:
            if Lambda == 1.0:
                cells = dfunc.find_Lambda_cells(df)
                par,perp,tau1,tau2 = dfunc.find_stress(cells)

                Lambda1_stress_perp.append(abs(perp.mean()))
                err_Lamba1[i,0] = abs(perp.std())/np.sqrt(len(perp))

                Lambda1_stress_par.append(abs(par.mean()))
                err_Lamba1[i,1] = abs(par.std())/np.sqrt(len(par))

                Lambda1_stress_shear1.append(abs(tau1.mean()))
                err_Lamba1[i,2] = abs(tau1.std())/np.sqrt(len(tau1))

                Lambda1_stress_shear2.append(abs(tau2.mean()))
                err_Lamba1[i,3] = abs(tau2.std())/np.sqrt(len(tau2))

            else:
                cells = dfunc.find_Lambda_cells(df,Lambda=Lambda)
                par,perp,tau1,tau2 = dfunc.find_stress(cells)
                Lambdanot1_stress_perp.append(abs(perp.mean()))
                err_Lambanot1[i,0] = abs(perp.std())/np.sqrt(len(perp))

                Lambdanot1_stress_par.append(abs(par.mean()))
                err_Lambanot1[i,1] = abs(par.std())/np.sqrt(len(par))


                Lambdanot1_stress_shear1.append(abs(tau1.mean()))
                err_Lambanot1[i,2] = abs(tau1.std())/np.sqrt(len(tau1))


                Lambdanot1_stress_shear2.append(abs(tau2.mean()))
                err_Lambanot1[i,0] = abs(tau2.std())/np.sqrt(len(tau2))

        i+=1
        
    for Lambda in Lambdas:
        if Lambda == 1.0:
            plt.errorbar(time_steps,Lambda1_stress_perp, yerr=err_Lamba1[:,0], fmt="x", color=dfunc.colour_Lambda(Lambda))
            plt.errorbar(time_steps,Lambda1_stress_par, yerr=err_Lamba1[:,1], fmt="o",color=dfunc.colour_Lambda(Lambda))
            plt.errorbar(time_steps,Lambda1_stress_shear1, yerr=err_Lamba1[:,2], fmt="*", color=dfunc.colour_Lambda(Lambda))
            plt.errorbar(time_steps,Lambda1_stress_shear2, yerr=err_Lamba1[:,3], fmt="v", color=dfunc.colour_Lambda(Lambda))
        else:
            plt.errorbar(time_steps,Lambdanot1_stress_perp, yerr=err_Lambanot1[:,0], fmt="x", color=dfunc.colour_Lambda(Lambda))
            plt.errorbar(time_steps,Lambdanot1_stress_par, yerr=err_Lambanot1[:,1], fmt="o", color=dfunc.colour_Lambda(Lambda))
            plt.errorbar(time_steps,Lambdanot1_stress_shear1, yerr=err_Lambanot1[:,2], fmt="*", color=dfunc.colour_Lambda(Lambda))
            plt.errorbar(time_steps,Lambdanot1_stress_shear2, yerr=err_Lambanot1[:,3], fmt="v", color=dfunc.colour_Lambda(Lambda))
    return Lambdas

def plot_pressure_time(data_dir):
    files = get_file_paths(data_dir)

    time_steps = []

    Lambda1_alpha = [] #deviatoric active stress
    Lambda1_p = []

    Lambdanot1_alpha = []
    Lambdanot1_p = []


    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        df = pd.read_csv(file_path, sep="\t")
        time_steps.append(time_step * 0.1)
        dfunc.find_stress(df)
        Lambdas = dfunc.find_Lambdas(df)
        for Lambda in Lambdas:
            if Lambda == 1.0:
                cells = dfunc.find_Lambda_cells(df)
                p,alpha = dfunc.find_pressure(cells)
                Lambda1_p.append(p.mean())
                Lambda1_alpha.append(alpha.mean())
                
            else:
                cells = dfunc.find_Lambda_cells(df,Lambda=Lambda)
                p,alpha = dfunc.find_pressure(cells)
                Lambdanot1_p.append(p.mean())
                Lambdanot1_alpha.append(alpha.mean())

        
    for Lambda in Lambdas:
        if Lambda == 1.0:
            plt.plot(time_steps,Lambda1_p, "x", color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambda1_alpha, "o",color=dfunc.colour_Lambda(Lambda))
        else:
            plt.plot(time_steps,Lambdanot1_p, "x", color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambdanot1_alpha, "o",color=dfunc.colour_Lambda(Lambda))
    return Lambdas

def plot_av_length_time(data_dir):
    files = get_file_paths(data_dir)
    Lambdas = sorted(dfunc.find_Lambdas(pd.read_csv(files[0], sep="\t")))
    if len(Lambdas) > 1:
        Lambdas = sorted(Lambdas)
        Lambdas = Lambdas[1]/Lambdas[0]
    else:
        Lambdas = Lambdas[0]


    time_steps = []
    av_ltot=[]

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        time_steps.append(time_step*0.1)
    
        df = pd.read_csv(file_path, sep="\t")
        av_ltot.append(dfunc.average_length(df))
        
        plt.plot(time_steps,av_ltot,color=dfunc.colour_Lambda(Lambdas))
    return [Lambdas]

def orientation_time(data_dir,width=40):
    files = get_file_paths(data_dir)
    time_steps =[]

    theta_Lambda1 = []
    theta_Lambda2 = []

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        time_steps.append(time_step*0.1)
    
        df = pd.read_csv(file_path, sep="\t")
        dfunc.remove_out_channel(df,width)
        Lambdas = sorted(dfunc.find_Lambdas(df))

        df1 = dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
        theta_Lambda1.append(abs(dfunc.get_orientation_cells(df1).mean()))

        if len(Lambdas) != 1:

            df2 = dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
            theta_Lambda2.append(abs(dfunc.get_orientation_cells(df2).mean()))

    plt.plot(time_steps,theta_Lambda1,color=dfunc.colour_Lambda(Lambdas[0]))
    if len(Lambdas) != 1:
        plt.plot(time_steps,theta_Lambda2,color=dfunc.colour_Lambda(Lambdas[1]))

    return Lambdas



def plot_COM_interacting(data_dir):
    files = get_file_paths(data_dir)

    time_steps = []

    com1 = [] 
    com2 = []
    comtot=[]
    t_collision = dfunc.colonies_collided(files) #finds time when colonies collide
    i=0
    for file_path in files:
        if i==0: #registers initial COM to take away so that cells are centered
            df = pd.read_csv(file_path, sep="\t")
            Lambdas = sorted(dfunc.find_Lambdas(df))
            xtot,ytot,ztot= dfunc.centerBiofilm(df)
            x_com_initial = xtot
            y_com_initial = ytot
            if len(Lambdas) != 1:
                xi1,yi1,zi1= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[0]))
                xi2,yi2,zi2= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[1]))
                
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))

        if time_step*0.1 < t_collision:
            continue
        time_steps.append(i * 0.1) #only appends after colonies have collided
        i+=1
    
        df = pd.read_csv(file_path, sep="\t")
        Lambdas = sorted(dfunc.find_Lambdas(df))
        xtot,ytot,ztot= dfunc.centerBiofilm(df)
        
        comtot.append(np.sqrt((xtot-x_com_initial)**2 + (ytot - y_com_initial)**2))

        if len(Lambdas) != 1:
            df1=dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
            x1,y1,z1= dfunc.centerBiofilm(df1)
            com1.append(np.sqrt((x1-xi1)**2 + (y1-yi1)**2))

            df2=dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
            x2,y2,z2= dfunc.centerBiofilm(df2)
            com2.append(np.sqrt((x2-xi2)**2 + (y2-yi2)**2))
            
    if len(Lambdas) == 1:
        plt.plot(time_steps,comtot, color=dfunc.colour_Lambda(Lambdas[0]))
    else:
        #plt.plot(time_steps,xcomtot, "--", color="k")
        plt.plot(time_steps,com1, color=dfunc.colour_Lambda(Lambdas[0]))
        plt.plot(time_steps,com2,color=dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def plot_perimeter_area_time(data_dir):
    time_steps, perimeter1, perimeter2, area1,area2, Lambdas = dfunc.perimeter_area_time(data_dir)
    Lambdas = sorted(Lambdas)

    plt.plot(time_steps, perimeter1, '.', color=dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(time_steps, area1, 'x', color=dfunc.colour_Lambda(Lambdas[0]))

    if len(Lambdas) != 1.0:
        plt.plot(time_steps, perimeter2, '.', color=dfunc.colour_Lambda(Lambdas[1]))
        plt.plot(time_steps, area2, 'x', color=dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def plot_perimeter_area_N(data_dir,ax,plot=True):
    time_steps, perimeter1, perimeter2, ptot, area1, area2, atot, Lambdas = dfunc.perimeter_area_time(data_dir)
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)
    Lambdas = sorted(Lambdas)

    colonytot = np.asarray([np.asarray(Lambda1_counts)+np.asarray(Lambda2_counts),np.asarray(ptot),np.asarray(atot)])
    Ppartot,Perrtot = dfunc.find_scaling_law(colonytot[0,:],colonytot[1,:])
    Apartot,Aerrtot = dfunc.find_scaling_law(colonytot[0,:],colonytot[2,:])
    
    #print(f"Total, A(N) = ({Apartot[0]}+-{Aerrtot[0]})*N^({Apartot[1]}+-{Aerrtot[1]})")
    print(f"Total, P(N) = ({Ppartot[0]}+-{Perrtot[0]})*N^({Ppartot[1]}+-{Perrtot[1]})")

    plot = False

    if len(Lambdas) !=1:
        colony1 = np.asarray([np.asarray(Lambda1_counts),np.asarray(perimeter1),np.asarray(area1)])
        Ppar1,Perr1 = dfunc.find_scaling_law(colony1[0,:],colony1[1,:])
        Apar1,Aerr1 = dfunc.find_scaling_law(colony1[0,:],colony1[2,:])
        #colony1[1,:] = (colony1[1,:])/par1[0]
        #print(f"Lambda = {Lambdas[0]}, A(N) = ({Apar1[0]}+-{Aerr1[0]})*N^({Apar1[1]}+-{Aerr1[1]})")
        print(f"Lambda = {Lambdas[0]}, P(N) = ({Ppar1[0]}+-{Perr1[0]})*N^({Ppar1[1]}+-{Perr1[1]})")

        colony2 = np.asarray([np.asarray(Lambda2_counts),np.asarray(perimeter2),np.asarray(area2)])
        Ppar2,Perr2 = dfunc.find_scaling_law(colony2[0,:],colony2[1,:])
        Apar2,Aerr2 = dfunc.find_scaling_law(colony2[0,:],colony2[2,:])
        #colony2[1,:] = (colony2[1,:])/par2[0]
        #print(f"Lambda = {Lambdas[1]}, A(N) = ({Apar2[0]}+-{Aerr2[0]})*N^({Apar2[1]}+-{Aerr2[1]})")
        print(f"Lambda = {Lambdas[1]}, P(N) = ({Ppar2[0]}+-{Perr2[0]})*N^({Ppar2[1]}+-{Perr2[1]})")
        if plot:
            n1 = np.geomspace(colony1[0,:].min(),colony1[0,:].max(),1000)
            fit1 = dfunc.linear_exp(n1,*Ppar1)#/par1[0]
            ax[0].scatter(colony1[0,:],colony1[1,:],marker=".",color = dfunc.colour_Lambda(Lambdas[0]))
            ax[0].plot(n1,fit1,color = dfunc.colour_Lambda(Lambdas[0]))

            fit1 = dfunc.linear_exp(n1,*Apar1)#/par1[0]
            ax[1].scatter(colony1[0,:],colony1[2,:],marker=".",color = dfunc.colour_Lambda(Lambdas[0]))
            ax[1].plot(n1,fit1,color = dfunc.colour_Lambda(Lambdas[0]))


            n2 = np.geomspace(colony2[0,:].min(),colony2[0,:].max(),1000)
            fit2 = dfunc.linear_exp(n2,*Ppar2)#/par2[0]
            ax[0].scatter(colony2[0,:],colony2[1,:],marker=".",color = dfunc.colour_Lambda(Lambdas[1]))
            ax[0].plot(n2,fit2,color = dfunc.colour_Lambda(Lambdas[1]))

            fit2 = dfunc.linear_exp(n2,*Apar2)#/par2[0]
            ax[1].scatter(colony2[0,:],colony2[2,:],marker=".",color = dfunc.colour_Lambda(Lambdas[1]))
            ax[1].plot(n2,fit2,color = dfunc.colour_Lambda(Lambdas[1]))

            return Lambdas,[Ppar1[0],Ppar2[0]],[Ppar1[1],Ppar2[1]],[Apar1[0],Apar2[0]],[Apar1[1],Apar2[1]]
    
    return Lambdas,[Ppar1[0]],[Ppar1[1]],[Apar1[0]],[Apar1[1]]

def perim_area_compare_FreeXInteracting(data_dir,ax):
    time_steps, perimeter1, perimeter2, area1, area2, Lambdas = dfunc.perimeter_area_time(data_dir)

    area1 = np.asarray(area1)
    area2 = np.asarray(area2)
    perimeter1 = np.asarray(perimeter1)
    perimeter2 = np.asarray(perimeter2)

    iQ1 = (area1*4*np.pi)/(perimeter1**2)
    iQ2 = (area2*4*np.pi)/(perimeter2**2)

    ratio1 = (perimeter1**2)/area1
    ratio2 = (perimeter2**2)/area2

    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)

    Lambda_params = [1.0,1.5,5.0,10.0]
    
    Pk_params = [27,24,26,26]
    Ak_params = [7.5,7.5,10.2,9.6]

    Ax_params = [0.946,0.94,0.882,0.88]
    Px_params = [0.429,0.45,0.415,0.41]

    for i,Lambda in enumerate(Lambdas):
        Lambda1 = round(Lambda,1)
        idx = Lambda_params.index(Lambda1)
        Ppar = [Pk_params[idx],Px_params[idx]]
        Apar = [Ak_params[idx],Ax_params[idx]]
        if i == 0:
            pFree = np.asarray(dfunc.linear_exp(np.asarray(Lambda1_counts),*Ppar))
            aFree = np.asarray(dfunc.linear_exp(np.asarray(Lambda1_counts),*Apar))
            iQFree = (aFree*4*np.pi)/(pFree**2)
            ratioFree = (pFree**2)/aFree
            ax[0].plot(Lambda1_counts,perimeter1/pFree,color=dfunc.colour_Lambda(Lambda))
            ax[1].plot(Lambda1_counts,area1/aFree,color=dfunc.colour_Lambda(Lambda))
            ax[2].plot(Lambda1_counts,iQ1/iQFree,color=dfunc.colour_Lambda(Lambda))
            ax[3].plot(Lambda1_counts,ratio1/ratioFree,color=dfunc.colour_Lambda(Lambda))
        else:
            pFree = np.asarray(dfunc.linear_exp(np.asarray(Lambda2_counts),*Ppar))
            aFree = np.asarray(dfunc.linear_exp(np.asarray(Lambda2_counts),*Apar))
            iQFree = (aFree*4*np.pi)/(pFree**2)
            ratioFree = (pFree**2)/aFree
            ax[0].plot(Lambda1_counts,perimeter2/pFree,color=dfunc.colour_Lambda(Lambda))
            ax[1].plot(Lambda1_counts,area2/aFree,color=dfunc.colour_Lambda(Lambda))
            ax[2].plot(Lambda1_counts,iQ2/iQFree,color=dfunc.colour_Lambda(Lambda))
            ax[3].plot(Lambda1_counts,ratio2/ratioFree,color=dfunc.colour_Lambda(Lambda))
    return Lambdas

def plot_perimeter_total_colony(data_dir):
    time_steps, perimeter1, perimeter2, ptot, area1, area2,atot, Lambdas = dfunc.perimeter_area_time(data_dir)
    Lambdas = sorted(Lambdas)
    files = get_file_paths(data_dir)
    t = dfunc.colonies_collided(files)
    plt.scatter(t,ptot[int(t*10)],c="r")
    plt.plot(time_steps, ptot, color="k")
    #plt.plot(time_steps, atot, 'x', color=dfunc.colour_Lambda(Lambdas[0]))
    return Lambdas



def plot_IQ_time(data_dir):
    time_steps, perimeter1, perimeter2, area1,area2, Lambdas = dfunc.perimeter_area_time(data_dir)
    Lambdas = sorted(Lambdas)
    perimeter1 = np.asarray(perimeter1)
    area1 = np.asarray(area1)
    iq1 = 4*np.pi*area1/perimeter1**2
    plt.plot(time_steps, iq1, '.', color=dfunc.colour_Lambda(Lambdas[0]))


    if len(Lambdas) != 1.0:
        perimeter2 = np.asarray(perimeter2)
        area2 = np.asarray(area2)
        iq2 = 4*np.pi*area2/perimeter2**2
        plt.plot(time_steps, iq2, '.', color=dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def surface_fraction_time(data_dir):
    time_steps, frac1, frac2, Lambdas = dfunc.surface_fraction(data_dir)
    plt.plot(time_steps, frac1, color=dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(time_steps, frac2, color=dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def interfacexsurface_time(data_dir):
    time_steps, frac_interface, frac_external, Lambdas = dfunc.surface_fraction_high(data_dir)
    plt.plot(time_steps, frac_interface, color=dfunc.colour_Lambda(Lambdas[1]))
    return [Lambdas[1]]

def plot_average_interfacexsurface(data_dirs,ax=None):
    raw1 = []
    raw2 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, frac_interface, frac_external, Lambdas = dfunc.surface_fraction_high(data_dir)

        raw1.append(frac_interface)
        raw2.append(frac_external)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = np.asarray(time_steps)
            max_len = len(time_steps)
    Lam = np.unique(Lambdas)
    fraction1 = []
    fraction2 = []

    if len(Lambdas) != 1:

        # 2. Pad the shorter sequences with NaN up to max_len
        for item1, item2 in zip(raw1, raw2):
            padded1 = list(item1) + [np.nan] * (max_len - len(item1))
            fraction1.append(padded1)
            
            padded2 = list(item2) + [np.nan] * (max_len - len(item2))
            fraction2.append(padded2)


        fraction1 = np.asarray(fraction1)
        avg1 = np.nanmean(fraction1, axis=0)
        count1 = np.sum(~np.isnan(fraction1), axis=0)
        
        # Only retain positions with >= 2 datapoints
        mask1 = count1 >= 2
        avg1 = np.nanmean(fraction1[:, mask1], axis=0)
        std1 = np.nanstd(fraction1[:, mask1], axis=0)/np.sqrt(count1[mask1])

        if ax == None:
            plt.errorbar(time[mask1], avg1,yerr=std1,fmt="-", color=dfunc.colour_Lambda(Lam[1]))
        else:
            ax.errorbar(time[mask1], avg1,yerr=std1,fmt="-", color=dfunc.colour_Lambda(Lam[1]))

        """fraction2 = np.asarray(fraction2)
        avg2 = np.nanmean(fraction2, axis=0)
        std2 =np.nanstd(fraction2, axis=0)/(np.sqrt(fraction2.shape[0]))
        plt.errorbar(time, avg2,yerr=std2,fmt="-", color=dfunc.colour_Lambda(Lam[1]))"""
    save = True
    if save:
        data = np.asarray([avg1,std1])
        header = "avg\tstdev"
        np.savetxt("Inner_fraction_Lambda"+str(Lam[1]), np.transpose(data), delimiter="\t", fmt="%g",header=header)
    return [Lam[1]]

def plot_average_interfacexsurface_ratio(data_dirs):
    raw1 = []
    raw2 = []

    ratio = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, frac_interface, frac_external, Lambdas = dfunc.surface_fraction_high(data_dir)
        raw1.append(frac_interface)
        raw2.append(frac_external)

        if abs(1 - Lambdas[1]/Lambdas[0]) < 0.1:
            ratio.append(1)
        else:
            ratio.append(Lambdas[1]/Lambdas[0])
        if len(Lambdas) != 1:
            raw1.append(frac_interface)
            raw2.append(frac_external)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
        
        if len(set(ratio)) != 1:
            print("more than one ratio:", ratio)


    fraction1 = []
    fraction2 = []

        # 2. Pad the shorter sequences with NaN up to max_len
    for item1, item2 in zip(raw1, raw2):
        padded1 = list(item1) + [np.nan] * (max_len - len(item1))
        fraction1.append(padded1)
            
        padded2 = list(item2) + [np.nan] * (max_len - len(item2))
        fraction2.append(padded2)

    fraction1 = np.asarray(fraction1)
    avg1 = np.nanmean(fraction1, axis=0)
    std1 =np.nanstd(fraction1, axis=0)/(np.sqrt(fraction1.shape[0]))

    fraction2 = np.asarray(fraction2)
    avg2 = np.nanmean(fraction2, axis=0)
    std2 =np.nanstd(fraction2, axis=0)/(np.sqrt(fraction2.shape[0]))

    errorbar = False
    if errorbar:
        #plt.errorbar(time, avg2,yerr=std2,fmt="-", color=dfunc.colour_Lambda(ratio[0]))
        plt.errorbar(time, avg1,yerr=std1,fmt="--", color=dfunc.colour_Lambda(ratio[0]))
    else:
        #plt.errorbar(time, avg2,fmt="-", color=dfunc.colour_Lambda(ratio[0]))
        plt.errorbar(time, avg1,fmt="--", color=dfunc.colour_Lambda(ratio[0]))


def plot_average_surface_fraction(data_dirs):
    raw1 = []
    raw2 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, frac1, frac2, Lambdas = dfunc.surface_fraction(data_dir)

        raw1.append(frac1)
        raw2.append(frac2)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = np.asarray(time_steps)
            max_len = len(time_steps)
    Lam = np.unique(Lambdas)
    fraction1 = []
    fraction2 = []

    if len(Lambdas) != 1:

        # 2. Pad the shorter sequences with NaN up to max_len
        for item1, item2 in zip(raw1, raw2):
            padded1 = list(item1) + [np.nan] * (max_len - len(item1))
            fraction1.append(padded1)
            
            padded2 = list(item2) + [np.nan] * (max_len - len(item2))
            fraction2.append(padded2)


        fraction1 = np.asarray(fraction1)
        avg1 = np.nanmean(fraction1, axis=0)
        std1 =np.nanstd(fraction1, axis=0)/(np.sqrt(fraction1.shape[0]))
        #plt.errorbar(time, avg1,yerr=std1,fmt="-", color=dfunc.colour_Lambda(Lam[0]))

        fraction2 = np.asarray(fraction2)
        # Number of valid datapoints at each position
        count2 = np.sum(~np.isnan(fraction2), axis=0)
        # Only retain positions with >= 2 datapoints
        mask2 = count2 >= 2
        avg2 = np.nanmean(fraction2[:, mask2], axis=0)
        std2 = np.nanstd(fraction2[:, mask2], axis=0) / np.sqrt(count2[mask2])
        plt.errorbar(time[mask2], avg2,yerr=std2,fmt="-", color=dfunc.colour_Lambda(Lam[1]))

        
    return Lambdas

def plot_average_surface_fraction_ratio(data_dirs):
    raw1 = []
    raw2 = []

    ratio = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, frac1, frac2, Lambdas = dfunc.surface_fraction(data_dir)
        raw1.append(frac1)
        raw2.append(frac2)

        if abs(1 - Lambdas[1]/Lambdas[0]) < 0.1:
            ratio.append(1)
        else:
            ratio.append(Lambdas[1]/Lambdas[0])
        if len(Lambdas) != 1:
            raw1.append(frac1)
            raw2.append(frac2)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
        
        if len(set(ratio)) != 1:
            print("more than one ratio:", ratio)


    fraction1 = []
    fraction2 = []

        # 2. Pad the shorter sequences with NaN up to max_len
    for item1, item2 in zip(raw1, raw2):
        padded1 = list(item1) + [np.nan] * (max_len - len(item1))
        fraction1.append(padded1)
            
        padded2 = list(item2) + [np.nan] * (max_len - len(item2))
        fraction2.append(padded2)

    fraction1 = np.asarray(fraction1)
    avg1 = np.nanmean(fraction1, axis=0)
    std1 =np.nanstd(fraction1, axis=0)/(np.sqrt(fraction1.shape[0]))

    fraction2 = np.asarray(fraction2)
    avg2 = np.nanmean(fraction2, axis=0)
    std2 =np.nanstd(fraction2, axis=0)/(np.sqrt(fraction2.shape[0]))

    errorbar = False
    if errorbar:
        plt.errorbar(time, avg2,yerr=std2,fmt="-", color=dfunc.colour_Lambda(ratio[0]))
        plt.errorbar(time, avg1,yerr=std1,fmt="--", color=dfunc.colour_Lambda(ratio[0]))
    else:
        plt.errorbar(time, avg2,fmt="-", color=dfunc.colour_Lambda(ratio[0]))
        plt.errorbar(time, avg1,fmt="--", color=dfunc.colour_Lambda(ratio[0]))


def plot_average_orientation(data_dirs,width=40):
    raw1 = []

    ratio = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        files = get_file_paths(data_dir)
        time_steps =[]

        theta_Lambda1 = []
        for file_path in files:
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))
            time_steps.append(time_step*0.1)
        
            df = pd.read_csv(file_path, sep="\t")
            dfunc.remove_out_channel(df,width)
            Lambdas = sorted(dfunc.find_Lambdas(df))
            theta_Lambda1.append(abs(dfunc.get_orientation_cells(df).mean()))
        raw1.append(theta_Lambda1)
        ratio.append(Lambdas[1]/Lambdas[0])
        

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = np.asarray(time_steps)
            max_len = len(time_steps)
    ratio = list(set(ratio))

    orientation1 = []

        # 2. Pad the shorter sequences with NaN up to max_len
    for item1 in raw1:
        padded1 = list(item1) + [np.nan] * (max_len - len(item1))
        orientation1.append(padded1)

    orientation1 = np.asarray(orientation1)
    avg1 = np.nanmean(orientation1, axis=0)
    std1 =np.nanstd(orientation1, axis=0)/(np.sqrt(orientation1.shape[0]))


    errorbar = True
    if errorbar:
        plt.errorbar(time, avg1,yerr=std1, color=dfunc.colour_Lambda(ratio[0]))
    else:
        plt.errorbar(time, avg1, color=dfunc.colour_Lambda(ratio[0]))
    return ratio


def COM_N(data_dir,Nlong1=100, Nlong2=100,plot=True):
    files = get_file_paths(data_dir)

    com1 = [] 
    com2 = []
    counts1 =[]
    counts2 = []

    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)

    Npoints_fit = 5
    if len(Lambdas)>1:
        if Lambda1_counts[-Npoints_fit] < Nlong1  or Lambda2_counts[-Npoints_fit] < Nlong2: #at least 4 data points are higher than Nlong
            return [0],[0],[0]
    else:
        if Lambda1_counts[-Npoints_fit] < Nlong1: #at least 4 data points are higher than Nlong
            return [0],[0],[0]

    for i, file_path in enumerate(files):
        if i==0: #registers initial COM to take away so that cells are centered
            df = pd.read_csv(file_path, sep="\t")
            Lambdas = sorted(dfunc.find_Lambdas(df))
    
            xi1,yi1,zi1= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[0]))
            if len(Lambdas)>1:
                xi2,yi2,zi2= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[1]))
                    
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue

        if len(Lambdas)>1:
            if Lambda1_counts[i]<=Nlong1 or Lambda2_counts[i]<=Nlong2:
                continue
        else:
            if Lambda1_counts[i]<=Nlong1:
                continue

        df = pd.read_csv(file_path, sep="\t")

        df1=dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
        x1,y1,z1= dfunc.centerBiofilm(df1)
        com1.append(np.sqrt((x1-xi1)**2 + (y1-yi1)**2))
        counts1.append(Lambda1_counts[i])
        if len(Lambdas)>1:
            df2=dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
            x2,y2,z2= dfunc.centerBiofilm(df2)
            com2.append(np.sqrt((x2-xi2)**2 + (y2-yi2)**2))
            counts2.append(Lambda2_counts[i])

    
    
    Lambdas = sorted(Lambdas)

    Lambda1_counts = Lambda1_counts[len(com1):]
    colony1 = np.asarray([np.asarray(counts1),np.asarray(com1)])
    par1,err1 = dfunc.find_scaling_law(colony1[0,:],colony1[1,:])
    colony1[1,:] = (colony1[1,:]-par1[2])/par1[0]

    if len(Lambdas)>1:
        colony2 = np.asarray([np.asarray(counts2),np.asarray(com2)])
        Lambda2_counts = Lambda2_counts[len(com2):]
        par2,err2 = dfunc.find_scaling_law(colony2[0,:],colony2[1,:])
        colony2[1,:] = (colony2[1,:]-par2[2])/par2[0]

    

    if plot:
        n1 = np.geomspace(colony1[0,:].min(),colony1[0,:].max(),1000)
        fit1 = (dfunc.linear_exp(n1,*par1)-par1[2])/par1[0]
        #print(f"Lambda = {Lambdas[0]}, R = ({par1[0]}+-{err1[0]})*N^({par1[1]}+-{err1[1]})")
        plt.scatter(colony1[0,:],colony1[1,:],marker=".",color = dfunc.colour_Lambda(Lambdas[0]))
        #plt.plot(colony1[0,:],colony1[1,:],color = dfunc.colour_Lambda(Lambdas[0]),alpha=0.5)
        plt.plot(n1,fit1,color = dfunc.colour_Lambda(Lambdas[0]))


        if len(Lambdas)>1:
            n2 = np.geomspace(colony2[0,:].min(),colony2[0,:].max(),1000)
            fit2 = (dfunc.linear_exp(n2,*par2)-par2[2])/par2[0]
            #print(f"Lambda = {Lambdas[1]}, R = ({par2[0]}+-{err2[0]})*N^({par2[1]}+-{err2[1]})")
            plt.scatter(colony2[0,:],colony2[1,:],marker=".",color = dfunc.colour_Lambda(Lambdas[1]))
            #plt.plot(colony2[0,:],colony2[1,:],color = dfunc.colour_Lambda(Lambdas[1]))
            plt.plot(n2,fit2,color = dfunc.colour_Lambda(Lambdas[1]))

    if len(Lambdas)>1:
        return Lambdas,[par1[0],par2[0]],[par1[1],par2[1]]
    else:
        return Lambdas,[par1[0]],[par1[1]]
    




#plots quantity vs distance
def plot_stress_distance(filepath):
    df = pd.read_csv(filepath, sep="\t")
    cells = ut.getCells(filepath)
    Lambdas = dfunc.find_Lambdas(cells)

    for Lambda in Lambdas:
        cells_Lambda = dfunc.find_Lambda_cells(cells,Lambda)
        df_Lambda = dfunc.find_Lambda_cells(df,Lambda)

        par_stress, perp_stress, shear1, shear2 = dfunc.find_stress(df_Lambda)
        par_stress = np.asarray(par_stress)
        perp_stress = np.asarray(perp_stress)
        shear1 = np.asarray(shear1)
        shear2 = np.asarray(shear2)


        
        cells_Lambda = dfunc.centerCells(cells_Lambda)
        d_cells = np.asarray(dfunc.distance_from_origin(cells_Lambda))
        dmax = max(d_cells)
        dmin = min(d_cells)

        n_points = 10
        stress = np.zeros((5,n_points))
        stress[0,:] = np.linspace(dmin,dmax,n_points)
        std_err = np.zeros((4,n_points))

        for i in range(0,len(stress[0,:])-1):
            n_cells_sqrt = np.sqrt(np.sum(np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))))

            stress[1,i] = par_stress[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            std_err[0,i] = par_stress[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].std()/n_cells_sqrt

            stress[2,i] = perp_stress[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            std_err[1,i] = perp_stress[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].std()/n_cells_sqrt

            stress[3,i] = shear1[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            std_err[2,i] = shear1[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].std()/n_cells_sqrt

            stress[4,i] = shear2[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            std_err[3,i] = shear2[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].std()/n_cells_sqrt


        c_Lambda = dfunc.colour_Lambda(Lambda)
        plt.errorbar(stress[0,:], stress[2,:],yerr=std_err[1,i],fmt= "x",color = c_Lambda)
        plt.errorbar(stress[0,:], stress[1,:],yerr=std_err[0,i],fmt=  "o",color = c_Lambda)
        plt.errorbar(stress[0,:], stress[3,:],yerr=std_err[2,i],fmt=  "*",color = c_Lambda)
        plt.errorbar(stress[0,:], stress[4,:],yerr=std_err[3,i],fmt=  "v",color = c_Lambda)
    return Lambdas

def plot_pressure_distance(filepath):
    df = pd.read_csv(filepath, sep="\t")
    cells = ut.getCells(filepath)
    Lambdas = dfunc.find_Lambdas(cells)

    for Lambda in Lambdas:
        cells_Lambda = dfunc.find_Lambda_cells(cells,Lambda)
        df_Lambda = dfunc.find_Lambda_cells(df,Lambda)

        p,alpha = dfunc.find_pressure(df_Lambda)
        p = np.asarray(p)
        alpha = np.asarray(alpha)

        cells_Lambda = dfunc.centerCells(cells_Lambda)
        d_cells = np.asarray(dfunc.distance_from_origin(cells_Lambda))
        dmax = max(d_cells)
        dmin = min(d_cells)

        n_points = 10
        stress = np.zeros((3,n_points))
        stress[0,:] = np.linspace(dmin,dmax,n_points)
        std_err = np.zeros((2,n_points))
 

        for i in range(0,len(stress[0,:])-1):
            n_cells_sqrt = np.sqrt(len(d_cells))

            stress[1,i] = p[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            std_err[0,i] = p[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].std()/n_cells_sqrt

            stress[2,i] = alpha[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            std_err[1,i] = alpha[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].std()/n_cells_sqrt

        c_Lambda = dfunc.colour_Lambda(Lambda)
        plt.errorbar(stress[0,:], stress[1,:],yerr=std_err[1,i],fmt= "x",color = c_Lambda)
        plt.errorbar(stress[0,:], stress[2,:],yerr=std_err[0,i],fmt=  "o",color = c_Lambda)

    return Lambdas

def plot_fraction_distance(data_dir):
    files = get_file_paths(data_dir)
    d_normalised,fraction_higher, ratio= dfunc.find_fraction_distance(files[-1])
    plt.plot(d_normalised, fraction_higher,"--",color = dfunc.colour_ratio(ratio))
    return ratio

def average_fraction_distance(data_dirs,ax=None,com="total",width=False):
    #com = "total" or "higher" and it indicates which COM to be considered.
    #com = "total" is the total COM, the "higher" is the COM of the higher Lmabda colony

    fraction = []
    fraction_low =[]
    for data_dir in data_dirs:
        files = get_file_paths(data_dir)
        d_normalised,fraction_higher,fraction_lower, ratio= dfunc.find_fraction_distance(files[-1],n_points=50,com=com)
        fraction.append(fraction_higher)
        fraction_low.append(fraction_lower)
    fraction = np.asarray(fraction)
    fraction_low = np.asarray(fraction_low)

    avg_L1 = np.mean(fraction, axis=0)
    std_L1 =np.std(fraction, axis=0)/np.sqrt(len(data_dirs))
    avg_L2 = np.mean(fraction_low, axis=0)
    std_L2 =np.std(fraction_low, axis=0)/np.sqrt(len(data_dirs))


    if not width:
        if ax == None:
            if com == "higher":
                plt.errorbar(d_normalised, avg_L2,yerr=std_L2, color=dfunc.colour_Lambda(ratio))
            if com=="total":
                plt.errorbar(d_normalised, avg_L1,yerr=std_L1, color=dfunc.colour_Lambda(ratio))

            return ratio

        ax.errorbar(d_normalised, avg_L1,yerr=std_L1, color=dfunc.colour_Lambda(ratio))
        if com == "higher":
            ax.errorbar(d_normalised, avg_L2,"--",yerr=std_L2, color=dfunc.colour_Lambda(ratio))
        return ratio
    else:
        avg_L2 = avg_L2
        d_normalised = d_normalised-0.5
        params,error = dfunc.find_tanh_fraction_dist(d_normalised,avg_L2)
        print(f"ratio={ratio}: f2(rho) = ({round(params[1],2)}+-{round(error[1],2)})*(0.5+tanh((rho-0.5)/({round(params[0],2)}+-{round(error[0],2)})))")
        return params,error,ratio




#plots averages
def plot_average_growth(data_dirs,ax=None):
    raw_L1 = []
    raw_Ln1 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
        
        raw_L1.append(Lambda1_counts)
        raw_Ln1.append(Lambdanot1_counts)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    L1_count = []
    Ln1_count = []

    # 2. Pad the shorter sequences with NaN up to max_len
    for item1, itemnot1 in zip(raw_L1, raw_Ln1):
        # Pad Lambda 1 counts
        padded_L1 = list(item1) + [np.nan] * (max_len - len(item1))
        L1_count.append(padded_L1)
        
        # Pad Lambda not 1 counts
        padded_Ln1 = list(itemnot1) + [np.nan] * (max_len - len(itemnot1))
        Ln1_count.append(padded_Ln1)
        
    L1_count = np.asarray(L1_count)
    Ln1_count = np.asarray(Ln1_count)
    
    # 4. Calculate the averages ignoring the NaNs
    if ax == None:
        avg_L1 = np.nanmean(L1_count, axis=0)
        std_L1 =np.nanstd(L1_count, axis=0)/(avg_L1*np.log(2)*np.sqrt(L1_count.shape[0]))
        params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(avg_L1, time_step=0.1)
        print(f"Behaviour for Lambda={Lam[0]}: log_2 N(t) =({round(params_Lambda1[1],3)}+-{round(err_Lambda1[1],3)}) + t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
        if len(avg_L1) == len(time):
            plt.errorbar(time, np.log2(avg_L1),fmt= '.',yerr=std_L1, markersize=6, color=dfunc.colour_Lambda(Lam[0]))
            plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda1), '-', color=dfunc.colour_Lambda(Lam[0]),label=f"{Lam[0]}: ${round(params_Lambda1[0],3)}\pm{round(err_Lambda1[0],3)}$")
        if len(Lam) != 1:
            avg_L2 = np.nanmean(Ln1_count, axis=0)
            std_L2 =np.nanstd(Ln1_count, axis=0)/(avg_L2*np.log(2)*np.sqrt(Ln1_count.shape[0]))
            params_Lambda2, err_Lambda2 = dfunc.estimate_growth_rate(avg_L2, time_step=0.1)
            print(f"Behaviour for Lambda={Lam[1]}: log_2 N(t) =({round(params_Lambda2[1],3)}+-{round(err_Lambda2[1],3)}) + t/({round(params_Lambda2[0],3)}+-{round(err_Lambda2[0],3)})")
            if len(avg_L2) == len(time):
                plt.errorbar(time, np.log2(avg_L2),fmt= '.',yerr=std_L2, markersize=6, color=dfunc.colour_Lambda(Lam[1]))
                plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda2), '-', color=dfunc.colour_Lambda(Lam[1]),label=f"{Lam[1]}: ${round(params_Lambda2[0],3)}\pm{round(err_Lambda2[0],3)}$")
        return Lam
    
    avg_L1 = np.nanmean(L1_count, axis=0)
    std_L1 =np.nanstd(L1_count, axis=0)/(avg_L1*np.log(2)*np.sqrt(L1_count.shape[0]))
    params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(avg_L1, time_step=0.1)
    print(f"Behaviour for Lambda={Lam[0]}: log_2 N(t) =({round(params_Lambda1[1],3)}+-{round(err_Lambda1[1],3)}) + t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
    if len(avg_L1) == len(time):
        plt.errorbar(time, np.log2(avg_L1),fmt= '.',yerr=std_L1, markersize=6, color=dfunc.colour_Lambda(Lam[0]))
        plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda1), '-', color=dfunc.colour_Lambda(Lam[0]),label=f"{Lam[0]}: ${round(params_Lambda1[0],3)}\pm{round(err_Lambda1[0],3)}$")
    if len(Lam) != 1:
        avg_L2 = np.nanmean(Ln1_count, axis=0)
        std_L2 =np.nanstd(Ln1_count, axis=0)/(avg_L2*np.log(2)*np.sqrt(Ln1_count.shape[0]))
        params_Lambda2, err_Lambda2 = dfunc.estimate_growth_rate(avg_L2, time_step=0.1)
        print(f"Behaviour for Lambda={Lam[1]}: log_2 N(t) =({round(params_Lambda2[1],3)}+-{round(err_Lambda2[1],3)}) + t/({round(params_Lambda2[0],3)}+-{round(err_Lambda2[0],3)})")
        if len(avg_L2) == len(time):
            plt.errorbar(time, np.log2(avg_L2),fmt= '.',yerr=std_L2, markersize=6, color=dfunc.colour_Lambda(Lam[1]))
            plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda2), '-', color=dfunc.colour_Lambda(Lam[1]),label=f"{Lam[1]}: ${round(params_Lambda2[0],3)}\pm{round(err_Lambda2[0],3)}$")

    return Lam, params_Lambda1[1], err_Lambda1[1]

def plot_average_rg(data_dirs,ax=None):
    raw_L1 = []
    raw_Ln1 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:


        time_steps, Rg_tot, Lambdas,  Lambda1_counts, Lambdanot1_counts= dfunc.RgLambda_time(data_dir)
        if len(Lambdas) == 1:
            Lambda1_counts = Rg_tot
        
        raw_L1.append(Lambda1_counts)
        raw_Ln1.append(Lambdanot1_counts)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    L1_count = []
    Ln1_count = []
    # 2. Pad the shorter sequences with NaN up to max_len
    if Lam != 1:
        for item1, itemnot1 in zip(raw_L1, raw_Ln1):
            # Pad Lambda 1 counts
            padded_L1 = list(item1) + [np.nan] * (max_len - len(item1))
            L1_count.append(padded_L1)
            
            # Pad Lambda not 1 counts
            padded_Ln1 = list(itemnot1) + [np.nan] * (max_len - len(itemnot1))
            Ln1_count.append(padded_Ln1)
    else:
        for item1 in raw_L1:
            # Pad Lambda 1 counts
            padded_L1 = list(item1) + [np.nan] * (max_len - len(item1))
            L1_count.append(padded_L1)
        
    L1_count = np.asarray(L1_count)
    Ln1_count = np.asarray(Ln1_count)
    time = np.asarray(time)
    
    # 4. Calculate the averages ignoring the NaNs
    if ax == None:
        avg_L1 = np.nanmean(L1_count, axis=0)
        mask = np.where(avg_L1 != 0)
        std_L1 =np.nanstd(L1_count, axis=0)/(avg_L1*np.log(2)*np.sqrt(L1_count.shape[0]))
        params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(avg_L1[mask[0]], time_step=0.1)
        print(f"Behaviour for Lambda={Lam[0]}: log_2 Rg(t) =({round(params_Lambda1[1],3)}+-{round(err_Lambda1[1],3)}) + t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
        if len(avg_L1) == len(time):
            plt.errorbar(time[mask[0]], np.log2(avg_L1[mask[0]]),fmt= '.',yerr=std_L1[mask[0]], markersize=6, color=dfunc.colour_Lambda(Lam[0]))
            plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda1), '-', color=dfunc.colour_Lambda(Lam[0]),label=f"{Lam[0]}: ${round(params_Lambda1[0],3)}\pm{round(err_Lambda1[0],3)}$")
        if len(Lam) != 1:
            avg_L2 = np.nanmean(Ln1_count, axis=0)
            std_L2 =np.nanstd(Ln1_count, axis=0)/(avg_L2*np.log(2)*np.sqrt(Ln1_count.shape[0]))
            params_Lambda2, err_Lambda2 = dfunc.estimate_growth_rate(avg_L2, time_step=0.1)
            print(f"Behaviour for Lambda={Lam[1]}: log_2 Rg(t) =({round(params_Lambda2[1],3)}+-{round(err_Lambda2[1],3)}) + t/({round(params_Lambda2[0],3)}+-{round(err_Lambda2[0],3)})")
            if len(avg_L2) == len(time):
                plt.errorbar(time, np.log2(avg_L2),fmt= '.',yerr=std_L2, markersize=6, color=dfunc.colour_Lambda(Lam[1]))
                plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda2), '-', color=dfunc.colour_Lambda(Lam[1]),label=f"{Lam[1]}: ${round(params_Lambda2[0],3)}\pm{round(err_Lambda2[0],3)}$")
        return Lam
    
    avg_L1 = np.nanmean(L1_count, axis=0)
    std_L1 =np.nanstd(L1_count, axis=0)/(avg_L1*np.log(2)*np.sqrt(L1_count.shape[0]))
    params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(avg_L1, time_step=0.1)
    print(f"Behaviour for Lambda={Lam[0]}: log_2 Rg(t) =({round(params_Lambda1[1],3)}+-{round(err_Lambda1[1],3)}) + t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
    if len(avg_L1) == len(time):
        plt.errorbar(time, np.log2(avg_L1),fmt= '.',yerr=std_L1, markersize=6, color=dfunc.colour_Lambda(Lam[0]))
        plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda1), '-', color=dfunc.colour_Lambda(Lam[0]),label=f"{Lam[0]}: ${round(params_Lambda1[0],3)}\pm{round(err_Lambda1[0],3)}$")
    if len(Lam) != 1:
        avg_L2 = np.nanmean(Ln1_count, axis=0)
        std_L2 =np.nanstd(Ln1_count, axis=0)/(avg_L2*np.log(2)*np.sqrt(Ln1_count.shape[0]))
        params_Lambda2, err_Lambda2 = dfunc.estimate_growth_rate(avg_L2, time_step=0.1)
        print(f"Behaviour for Lambda={Lam[1]}: log_2 Rg(t) =({round(params_Lambda2[1],3)}+-{round(err_Lambda2[1],3)}) + t/({round(params_Lambda2[0],3)}+-{round(err_Lambda2[0],3)})")
        if len(avg_L2) == len(time):
            plt.errorbar(time, np.log2(avg_L2),fmt= '.',yerr=std_L2, markersize=6, color=dfunc.colour_Lambda(Lam[1]))
            plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda2), '-', color=dfunc.colour_Lambda(Lam[1]),label=f"{Lam[1]}: ${round(params_Lambda2[0],3)}\pm{round(err_Lambda2[0],3)}$")

    return Lam, params_Lambda1[1], err_Lambda1[1]



def plot_average_asphericity(data_dirs):
    raw_L1 = []
    raw_L2 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, Lambdas, asph1, asph2 = dfunc.calc_asphericity(data_dir)
        
        raw_L1.append(asph1)
        raw_L2.append(asph2)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    L1_Asph = []
    L2_Asph = []

    # 2. Pad the shorter sequences with NaN up to max_len
    for item1, itemnot1 in zip(raw_L1, raw_L2):
        # Pad Lambda 1 asphericity
        padded_L1 = list(item1) + [np.nan] * (max_len - len(item1))
        L1_Asph.append(padded_L1)
        
        # Pad Lambda 2 asphericity
        if len(Lam) != 1 :
            padded_L2 = list(itemnot1) + [np.nan] * (max_len - len(itemnot1))
            L2_Asph.append(padded_L2)
        
    L1_Asph = np.asarray(L1_Asph)
    L2_Asph = np.asarray(L2_Asph)

    avg_Asph1 = np.nanmean(L1_Asph, axis=0)
    err_Asph1 = np.nanstd(L1_Asph, axis=0)/np.sqrt(L1_Asph.shape[0])
    plt.errorbar(time, avg_Asph1,fmt= 'o',yerr=err_Asph1, markersize=6, color=dfunc.colour_Lambda(Lam[0]))

    if len(Lam) != 1:
        avg_Asph2 = np.nanmean(L2_Asph, axis=0)
        err_Asph2 = np.nanstd(L2_Asph, axis=0)/np.sqrt(L2_Asph.shape[0])
        plt.errorbar(time, avg_Asph2,fmt= 'o',yerr=err_Asph2, markersize=6, color=dfunc.colour_Lambda(Lam[1]))

    return Lam

def plot_average_dasph(data_dirs):
    raw_L1 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, Lambdas, asph1, asph2 = dfunc.calc_asphericity(data_dir)
        
        raw_L1.append(np.asarray(asph1)-np.asarray(asph2))
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    L1_dAsph = []


    # 2. Pad the shorter sequences with NaN up to max_len
    for item in raw_L1:
        # Pad Lambda 1 asphericity
        padded_L1 = list(item) + [np.nan] * (max_len - len(item))
        L1_dAsph.append(padded_L1)

    L1_dAsph = np.asarray(L1_dAsph)


    avg_dAsph = np.nanmean(L1_dAsph, axis=0)
    err_dAsph = np.nanstd(L1_dAsph, axis=0)/np.sqrt(L1_dAsph.shape[0])
    plt.errorbar(time, avg_dAsph,fmt= 'o',yerr=err_dAsph, markersize=6, color=dfunc.colour_Lambda(Lam[1]))

    return Lam

def plot_average_COM(data_dirs):
    raw1 = []
    raw2 = []
    rawtot = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        files = get_file_paths(data_dir)

        time_steps = []

        com1 = [] 
        com2 = []
        comtot=[]
        t_collision = dfunc.colonies_collided(files) #finds time when colonies collide

        i=0
        for file_path in files:
            if i==0: #registers initial COM to take away so that cells are centered
                df = pd.read_csv(file_path, sep="\t")
                Lambdas = sorted(dfunc.find_Lambdas(df))
                xtot,ytot,ztot= dfunc.centerBiofilm(df)
                x_com_initial = xtot
                y_com_initial = ytot
                if len(Lambdas) != 1:
                    xi1,yi1,zi1= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[0]))
                    xi2,yi2,zi2= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[1]))
                    
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))

            if time_step*0.1 < t_collision:
                continue
            time_steps.append(i * 0.1) #only appends after colonies have collided
            i+=1

            df = pd.read_csv(file_path, sep="\t")
            Lambdas = sorted(dfunc.find_Lambdas(df))

            comtot.append(np.sqrt((xtot-x_com_initial)**2 + (ytot - y_com_initial)**2))

            if len(Lambdas) != 1:
                df1=dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
                x1,y1,z1= dfunc.centerBiofilm(df1)
                com1.append(np.sqrt((x1-xi1)**2 + (y1-yi1)**2))

                df2=dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
                x2,y2,z2= dfunc.centerBiofilm(df2)
                com2.append(np.sqrt((x2-xi2)**2 + (y2-yi2)**2))

        rawtot.append(comtot)
        Lam.append(Lambdas)
        if len(Lambdas) != 1:
            raw1.append(com1)
            raw2.append(com2)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = np.asarray(time_steps)
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    com_1 = []
    com_2 = []
    com_tot = []

    if len(Lambdas) != 1:

        # 2. Pad the shorter sequences with NaN up to max_len
        for item1, item2, itemtot in zip(raw1, raw2, rawtot):
            padded1 = list(item1) + [np.nan] * (max_len - len(item1))
            com_1.append(padded1)
            
            padded2 = list(item2) + [np.nan] * (max_len - len(item2))
            com_2.append(padded2)

            paddedtot = list(itemtot) + [np.nan] * (max_len - len(itemtot))
            com_tot.append(paddedtot)
        
        com_tot = np.asarray(com_tot)
        avgtot = np.nanmean(com_tot, axis=0)
        stdtot =np.nanstd(com_tot, axis=0)/(np.sqrt(com_tot.shape[0]))
        #plt.errorbar(time, avgtot,yerr=stdtot,fmt="--", color="k")

        com_1 = np.asarray(com_1)
        avg1 = np.nanmean(com_1, axis=0)
        std1 =np.nanstd(com_1, axis=0)/(np.sqrt(com_1.shape[0]))
        plt.errorbar(time, avg1,yerr=std1,fmt=".", color=dfunc.colour_Lambda(Lam[0]))
        par,err = dfunc.find_scaling_law(avg1)
        fit1 = dfunc.linear_exp(time,*par)
        print(f"Lambda = {Lam[0]}, R =({par[0]}+-{err[0]})+ ({par[1]}+-{err[1]})*t^({par[2]}+-{err[2]})")
        #print(f"Lambda = {Lam[0]}, R = e^t/({par}+-{err})")
        plt.plot(time,fit1,color=dfunc.colour_Lambda(Lam[0]))

        com_2 = np.asarray(com_2)
        avg2 = np.nanmean(com_2, axis=0)
        std2 =np.nanstd(com_2, axis=0)/(np.sqrt(com_2.shape[0]))
        plt.errorbar(time, avg2,yerr=std2,fmt=".", color=dfunc.colour_Lambda(Lam[1]))
        par,err = dfunc.find_scaling_law(avg2)
        fit2 = dfunc.linear_exp(time,*par)
        plt.plot(time,fit2,color=dfunc.colour_Lambda(Lam[1]))
        print(f"Lambda = {Lam[1]}, R =({par[0]}+-{err[0]})+ ({par[1]}+-{err[1]})*t^({par[2]}+-{err[2]})")
        #print(f"Lambda = {Lam[1]}, R = e^t/({par}+-{err})")
        
    else:
        for item in rawtot:
            paddedtot = list(item) + [np.nan] * (max_len - len(item))
            com_tot.append(paddedtot)

        com_tot = np.asarray(com_tot)

        avgtot = np.nanmean(com_tot, axis=0)
        stdtot =np.nanstd(com_tot, axis=0)/(np.sqrt(com_tot.shape[0]))
        plt.errorbar(time, avgtot,yerr=stdtot, color=dfunc.colour_Lambda(Lam[0]))
            
    return Lambdas

def plot_average_COM_ratio(data_dirs):
    raw1 = []
    raw2 = []

    ratio = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        files = get_file_paths(data_dir)

        time_steps = []

        com1 = [] 
        com2 = []
        t_collision = dfunc.colonies_collided(files) #finds time when colonies collide

        i=0
        df = pd.read_csv(files[0], sep="\t") #registers initial COM to take away so that cells are centered
        Lambdas = sorted(dfunc.find_Lambdas(df))
        xi1,yi1,zi1= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[0]))
        xi2,yi2,zi2= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[1]))

        for file_path in files:    
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))

            if time_step*0.1 < t_collision:
                continue
            time_steps.append(i * 0.1) #only appends after colonies have collided
            i+=1

            df = pd.read_csv(file_path, sep="\t")
            Lambdas = sorted(dfunc.find_Lambdas(df))

            if len(Lambdas) != 1:
                df1=dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
                x1,y1,z1= dfunc.centerBiofilm(df1)
                com1.append(np.sqrt((x1-xi1)**2 + (y1-yi1)**2))

                df2=dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
                x2,y2,z2= dfunc.centerBiofilm(df2)
                com2.append(np.sqrt((x2-xi2)**2 + (y2-yi2)**2))

        if abs(1 - Lambdas[1]/Lambdas[0]) < 0.1:
            ratio.append(1)
        else:
            ratio.append(Lambdas[1]/Lambdas[0])
        if len(Lambdas) != 1:
            raw1.append(com1)
            raw2.append(com2)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
        
        if len(set(ratio)) != 1:
            print("more than one ratio:", ratio)


    com_1 = []
    com_2 = []

        # 2. Pad the shorter sequences with NaN up to max_len
    for item1, item2 in zip(raw1, raw2):
        padded1 = list(item1) + [np.nan] * (max_len - len(item1))
        com_1.append(padded1)
            
        padded2 = list(item2) + [np.nan] * (max_len - len(item2))
        com_2.append(padded2)

    com_1 = np.asarray(com_1)
    avg1 = np.nanmean(com_1, axis=0)
    std1 =np.nanstd(com_1, axis=0)/(np.sqrt(com_1.shape[0]))

    com_2 = np.asarray(com_2)
    avg2 = np.nanmean(com_2, axis=0)
    std2 =np.nanstd(com_2, axis=0)/(np.sqrt(com_2.shape[0]))

    errorbar = False
    if errorbar:
        plt.errorbar(time, avg2,yerr=std2,fmt="-", color=dfunc.colour_Lambda(ratio[0]))
        plt.errorbar(time, avg1,yerr=std1,fmt="--", color=dfunc.colour_Lambda(ratio[0]))
    else:
        plt.errorbar(time, avg2,fmt="-", color=dfunc.colour_Lambda(ratio[0]))
        plt.errorbar(time, avg1,fmt="--", color=dfunc.colour_Lambda(ratio[0]))

def plot_average_IQ_time(data_dirs):
    raw_L1 = []
    raw_L2 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, perimeter1, perimeter2, area1,area2, Lambdas = dfunc.perimeter_area_time(data_dir)
        Lambdas = sorted(Lambdas)
        
        raw_L1.append(asph1)
        raw_L2.append(asph2)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    L1_Asph = []
    L2_Asph = []

    # 2. Pad the shorter sequences with NaN up to max_len
    for item1, itemnot1 in zip(raw_L1, raw_L2):
        # Pad Lambda 1 asphericity
        padded_L1 = list(item1) + [np.nan] * (max_len - len(item1))
        L1_Asph.append(padded_L1)
        
        # Pad Lambda 2 asphericity
        if len(Lam) != 1 :
            padded_L2 = list(itemnot1) + [np.nan] * (max_len - len(itemnot1))
            L2_Asph.append(padded_L2)
        
    L1_Asph = np.asarray(L1_Asph)
    L2_Asph = np.asarray(L2_Asph)

    avg_Asph1 = np.nanmean(L1_Asph, axis=0)
    err_Asph1 = np.nanstd(L1_Asph, axis=0)/np.sqrt(L1_Asph.shape[0])
    plt.errorbar(time, avg_Asph1,fmt= 'o',yerr=err_Asph1, markersize=6, color=dfunc.colour_Lambda(Lam[0]))

    if len(Lam) != 1:
        avg_Asph2 = np.nanmean(L2_Asph, axis=0)
        err_Asph2 = np.nanstd(L2_Asph, axis=0)/np.sqrt(L2_Asph.shape[0])
        plt.errorbar(time, avg_Asph2,fmt= 'o',yerr=err_Asph2, markersize=6, color=dfunc.colour_Lambda(Lam[1]))

    return Lam

def average_COM_N(data_dirs):
    n_colony1 = []
    n_colony2 = []
    com_colony1 = []
    com_colony2 = []

    Lam = []

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        files = get_file_paths(data_dir)

        com1 = [] 
        com2 = []
        t_collision = 1 #dfunc.colonies_collided(files) #finds time when colonies collide

        i=0
        for file_path in files:
            if i==0: #registers initial COM to take away so that cells are centered
                df = pd.read_csv(file_path, sep="\t")
                Lambdas = sorted(dfunc.find_Lambdas(df))
    
                if len(Lambdas) != 1:
                    xi1,yi1,zi1= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[0]))
                    xi2,yi2,zi2= dfunc.centerBiofilm(dfunc.find_Lambda_cells(df,Lambdas[1]))
                    
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))

            if time_step < int(10*t_collision):
                continue
            i+=1

            df = pd.read_csv(file_path, sep="\t")

            df1=dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
            x1,y1,z1= dfunc.centerBiofilm(df1)
            com1.append(np.sqrt((x1-xi1)**2 + (y1-yi1)**2))

            df2=dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
            x2,y2,z2= dfunc.centerBiofilm(df2)
            com2.append(np.sqrt((x2-xi2)**2 + (y2-yi2)**2))


        time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)
        n_colony1 = n_colony1+Lambda1_counts[int(t_collision*10):]
        n_colony2 = n_colony2+Lambda2_counts[int(t_collision*10):]
        Lambdas = sorted(Lambdas)

        com_colony1 = com_colony1 + com1
        com_colony2 = com_colony2 + com2

        Lam.append(Lambdas)
    colony1 = np.asarray([np.asarray(n_colony1),np.asarray(com_colony1)])
    colony2 = np.asarray([np.asarray(n_colony2),np.asarray(com_colony2)])

    Lam = np.unique(Lam)

    n1 = np.geomspace(100,colony1[0,:].max(),1000)
    par,err = dfunc.find_scaling_law(colony1[0,:],colony1[1,:])
    fit1 = dfunc.linear_exp(n1,*par)
    print(f"Lambda = {Lam[0]}, R =({par[0]}+-{err[0]})*N^({par[1]}+-{err[1]})")
    plt.scatter(colony1[0,:],colony1[1,:],marker=".",color = dfunc.colour_Lambda(Lam[0]))
    #plt.plot(colony1[0,:],colony1[1,:],color = dfunc.colour_Lambda(Lam[0]),alpha=0.5)
    plt.plot(n1,fit1,color = dfunc.colour_Lambda(Lam[0]))

    n2 = np.geomspace(100,colony1[0,:].max(),1000)
    par,err = dfunc.find_scaling_law(colony2[0,:],colony2[1,:])
    fit2 = dfunc.linear_exp(n2,*par)
    print(f"Lambda = {Lam[1]}, R =({par[0]}+-{err[0]})*N^({par[1]}+-{err[1]})")
    plt.scatter(colony2[0,:],colony2[1,:],marker=".",color = dfunc.colour_Lambda(Lam[1]))
    #plt.plot(colony2[0,:],colony2[1,:],color = dfunc.colour_Lambda(Lam[1]))
    plt.plot(n2,fit2,color = dfunc.colour_Lambda(Lam[1]))

    return Lambdas


#fancier plots
def t_doub_Lambda(data_dirs):

    tdoubl1=[]
    tdoubl2=[]

    for data_dir in data_dirs:
        time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)
        params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(Lambda1_counts, time_step=0.1)
        tdoubl1.append(params_Lambda1[0])

        if len(Lambdas)!=1:
            params_Lambda2, err_Lambda2 = dfunc.estimate_growth_rate(Lambda2_counts, time_step=0.1)
            tdoubl2.append(params_Lambda2[0])
    
    plt.errorbar(Lambdas[0],np.mean(tdoubl1),yerr=np.std(tdoubl1)/np.sqrt(len(tdoubl1)),fmt="x", color="k")
    plt.scatter(np.ones(len(tdoubl1))*Lambdas[0],tdoubl1,color="k",alpha=0.6,s=3)
                 
    if len(Lambdas)!=1:
        plt.errorbar(Lambdas[1],np.mean(tdoubl2),yerr=np.std(tdoubl2)/np.sqrt(len(tdoubl2)),fmt="x", color="k")
        plt.scatter(np.ones(len(tdoubl2))*Lambdas[1],tdoubl2,color="k",alpha=0.6,s=3)
    
    return Lambdas

def plot_yfraction(filepath,width=0):
    df = pd.read_csv(filepath, sep="\t")
    Lambdas = sorted(dfunc.find_Lambdas(df))

    df = dfunc.centerCells(df)

    ycells = np.asarray(df["pos_y"])
    
    df1 = dfunc.find_Lambda_cells(df,Lambdas[0])
    df2 = dfunc.find_Lambda_cells(df,Lambdas[1])
    ycells1 = np.asarray(abs(df1["pos_y"])) #reflect cells 
    ycells2 = np.asarray(abs(df2["pos_y"]))
    
    ymax = width/2
    if width == 0:
        ymax = max(ycells)

    n_points = 20
    fraction = np.zeros((3,n_points))
    fraction[0,:] = np.linspace(0,ymax,n_points)

    for i in range(1,len(fraction[0,:])):
        n1 = len(ycells1[np.where((ycells1>=fraction[0,i-1]) & (ycells1<=fraction[0,i]))])
        n2 = len(ycells2[np.where((ycells2>=fraction[0,i-1]) & (ycells2<=fraction[0,i]))])
        frac1,frac2 = dfunc.calc_fraction(n1,n2)
        fraction[1,i] = frac1
        fraction[2,i] = frac2
    
    plt.plot(fraction[0,1:], fraction[1,1:],color = dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(fraction[0,1:], fraction[2,1:],color = dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def yfraction_repeats(data_dirs,width=0):
    filepaths = []
    Lambdas = []
    ymax = width/2

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        files = get_file_paths(data_dir)
        filepaths.append(files[-1])

        df = pd.read_csv(files[-1], sep="\t")
        ycells = np.asarray(abs(df["pos_y"]))
        Lambdas = sorted(dfunc.find_Lambdas(df))
        if max(ycells) > ymax and width ==0:
            ymax = max(ycells)

    n_points = 20
    fraction = np.zeros((3,n_points))
    fraction[0,:] = np.linspace(0,ymax,n_points)
    n1=np.zeros((n_points))
    n2=np.zeros((n_points))
    for filepath in filepaths:
        df = pd.read_csv(filepath, sep="\t")
        df1 = dfunc.find_Lambda_cells(df,Lambdas[0])
        df2 = dfunc.find_Lambda_cells(df,Lambdas[1])
        ycells1 = np.asarray(abs(df1["pos_y"]))
        ycells2 = np.asarray(abs(df2["pos_y"]))

        for i in range(1,len(fraction[0,:])):
            n1[i] += len(ycells1[np.where((ycells1>=fraction[0,i-1]) & (ycells1<=fraction[0,i]))])
            n2[i] += len(ycells2[np.where((ycells2>=fraction[0,i-1]) & (ycells2<=fraction[0,i]))])
    frac1,frac2 = dfunc.calc_fraction(n1,n2)
    fraction[1,:] = frac1
    fraction[2,:] = frac2       
    plt.plot(fraction[0,1:], fraction[1,1:],"--",color = dfunc.colour_Lambda(Lambdas[1]))
    plt.plot(fraction[0,1:], fraction[2,1:],color = dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas 

def plot_xfraction(filepath,width=0,marker="."):
    df = pd.read_csv(filepath, sep="\t")
    Lambdas = sorted(dfunc.find_Lambdas(df))

    df = dfunc.centerCells(df)

    xcells = np.asarray(df["pos_x"])
    
    df1 = dfunc.find_Lambda_cells(df,Lambdas[0])
    df2 = dfunc.find_Lambda_cells(df,Lambdas[1])
    xcells1 = np.asarray(abs(df1["pos_x"])) #reflect cells 
    xcells2 = np.asarray(abs(df2["pos_x"]))
    
    xmax = width*1.5/2
    if width == 0:
        xmax = max(xcells)

    n_points = 20
    fraction = np.zeros((3,n_points))
    fraction[0,:] = np.linspace(0,xmax,n_points)

    for i in range(1,len(fraction[0,:])):
        n1 = len(xcells1[np.where((xcells1>=fraction[0,i-1]) & (xcells1<=fraction[0,i]))])
        n2 = len(xcells2[np.where((xcells2>=fraction[0,i-1]) & (xcells2<=fraction[0,i]))])
        frac1,frac2 = dfunc.calc_fraction(n1,n2)
        fraction[1,i] = frac1
        fraction[2,i] = frac2
    
    plt.plot(fraction[0,1:], fraction[1,1:],marker,color = dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(fraction[0,1:], fraction[2,1:],marker,color = dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas


#channel plots
def plot_leaving_position(ax,data_dir,plot_tot=False,width=60):
    df = get_exit_data(data_dir)
    if plot_tot:
        ytot = df["y"]
        ax.hist(ytot,color="k", bins=int(width/5),range=(-width/2,width/2),histtype='step',label="total")
    Lambdas = sorted(dfunc.find_Lambdas(df))
    for Lambda in Lambdas:
        df1 = dfunc.find_Lambda_cells(df,Lambda)
        y1 = df1["y"]
        ax.hist(y1,color=dfunc.colour_Lambda(Lambda), bins=int(width/5),range=(-width/2,width/2),histtype='step',label="$\Lambda=$"+str(Lambda))

    ax.legend()
    ax.set_xlabel("y pos (microns)")
    ax.set_xlim(-(width/2 +1),(width/2 +1))

    xticks= np.linspace(-(width/2),(width/2), 7)
    ax.set_xticks(xticks)

def plot_leaving_cells(data_dir):
    files = get_file_paths(data_dir)
    df = get_exit_data(data_dir)
    time = np.arange(2.9,len(files)*0.1+0.05,0.1)

    Lambdas = sorted(dfunc.find_Lambdas(df))

    if len(Lambdas) != 2:
        print("Cells haven't left")
        return 0

    df1 = dfunc.find_Lambda_cells(df,Lambdas[0])
    df2 = dfunc.find_Lambda_cells(df,Lambdas[1])

    time1 = np.asarray(df1["time"])
    time2 = np.asarray(df2["time"])
    fraction = np.zeros((3,len(time)-1))
    fraction[0,:] = time[1:]

    for i in range(0,len(fraction[0,:])):
        n1 = len(time1[np.where((time1>=time[i]) & (time1<=time[i+1]))])
        n2 = len(time2[np.where((time2>=time[i]) & (time2<=time[i+1]))])
        frac1,frac2 = dfunc.calc_fraction(n1,n2)
        fraction[1,i] = frac1
        fraction[2,i] = frac2

    ratio = Lambdas[1]/Lambdas[0]
    
    plt.plot(fraction[0,:], fraction[2,:],color = dfunc.colour_ratio(ratio))
    return [ratio]


def plot_fraction_time(data_dir, width=60):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir,channels=True, width=width)
    Lambdas = sorted(Lambdas)
    frac1,frac2 = dfunc.calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))

    plt.plot(time_steps,frac1,color=dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(time_steps,frac2,color=dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def plot_fraction_time_ratio(data_dir, width=40):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir,channels=True, width=width)
    Lambdas = sorted(Lambdas)
    ratio = Lambdas[1]/Lambdas[0]
    frac1,frac2 = dfunc.calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))

    #plt.plot(time_steps,frac1,"--",color=dfunc.colour_Lambda(ratio))
    plt.plot(np.asarray(time_steps)-3,frac2,"-",color=dfunc.colour_ratio(ratio))
    return ratio

def plot_average_fraction(data_dirs,width=60):
    raw_L1 = []
    raw_L2 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir,channels=True, width=width)
        frac1,frac2 = dfunc.calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))
        
        raw_L1.append(frac1)
        raw_L2.append(frac2)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    L1_frac = []
    L2_frac = []

    # 2. Pad the shorter sequences with NaN up to max_len
    for item1, itemnot1 in zip(raw_L1, raw_L2):
        # Pad Lambda 1 asphericity
        padded_L1 = list(item1) + [np.nan] * (max_len - len(item1))
        L1_frac.append(padded_L1)
        
        # Pad Lambda 2 asphericity
        if len(Lam) != 1 :
            padded_L2 = list(itemnot1) + [np.nan] * (max_len - len(itemnot1))
            L2_frac.append(padded_L2)
        
    L1_frac = np.asarray(L1_frac)
    L2_frac = np.asarray(L2_frac)

    avg_frac1 = np.nanmean(L1_frac, axis=0)
    err_frac1 = np.nanstd(L1_frac, axis=0)/np.sqrt(L1_frac.shape[0])
    plt.errorbar(time, avg_frac1,fmt= 'o',yerr=err_frac1, markersize=6, color=dfunc.colour_Lambda(Lam[0]))

    if len(Lam) != 1:
        avg_frac2 = np.nanmean(L2_frac, axis=0)
        err_frac2 = np.nanstd(L2_frac, axis=0)/np.sqrt(L2_frac.shape[0])
        plt.errorbar(time, avg_frac2,fmt= 'o',yerr=err_frac2, markersize=6, color=dfunc.colour_Lambda(Lam[1]))

    return Lam

def plot_initial_final(data_dir,width=40):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir,channels=True, width=width)
    Lambdas = np.asarray(sorted(Lambdas))
    frac1,frac2 = dfunc.calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))
    
    plt.scatter(frac1[-1],frac1[0],marker=".",color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]))
    plt.scatter(frac2[-1],frac2[0],marker="x",color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]))
    return [Lambdas[Lambdas != 1.0][0]]

def plot_whisker_initial_final(data_dirs,width=40,ax=None):
    initial_higher = []
    final_lower= []
    final_higher = []
    Lambdas = []
    for data_dir in data_dirs:
        time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir,channels=True, width=width)
        Lambdas += sorted(Lambdas)
        frac1,frac2 = dfunc.calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))
        initial_higher.append(frac2[0])
        final_lower.append(frac1[-1])
        final_higher.append(frac2[-1])

    initial_ratios = np.sort(np.unique(initial_higher))
    final_lower = np.asarray(final_lower)
    final_higher = np.asarray(final_higher)
    initial_higher = np.asarray(initial_higher)
    Lambdas = np.unique(np.asarray(Lambdas))
    for initial in initial_ratios:
        mask = np.where(initial_higher == initial)[0]
        if ax == None:
            bplot = plt.boxplot([final_lower[mask],final_higher[mask]],positions=np.ones(2)*initial,patch_artist=True,widths=np.ones(2)*0.08)
            
        else:
            bplot = ax.boxplot([final_lower[mask],final_higher[mask]],positions=np.ones(2)*initial,patch_artist=True,widths=np.ones(2)*0.08)

        for i,patch in enumerate(bplot['boxes']):
            patch.set_facecolor(dfunc.colour_Lambda(Lambdas[i]))
    plt.xticks(initial_ratios, [f"{x:.2f}" for x in initial_ratios])
    plt.xlim(initial_ratios[0]-0.05, initial_ratios[-1]+0.05)
    return Lambdas

def fraction_center(data_dir,width=40):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    
    time_steps = []
    Lambda1_frac= []
    Lambda2_frac = []

    Lambdas = sorted(dfunc.find_Lambdas(pd.read_csv(files[0], sep="\t")))
    ratio = Lambdas[1]/Lambdas[0]

    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        time_steps.append(time_step * 0.1)
        df = pd.read_csv(file_path, sep="\t")
        Lambdas = sorted(dfunc.find_Lambdas(df))

        df = dfunc.centerCells(df)
        
        df1 = dfunc.find_Lambda_cells(df,Lambdas[0])
        df2 = dfunc.find_Lambda_cells(df,Lambdas[1])
        xcells1 = np.asarray(abs(df1["pos_x"]))
        xcells2 = np.asarray(abs(df2["pos_x"]))

        xmax = width*0.1

        n1 = len(xcells1[np.where((xcells1>=-xmax) & (xcells1<=xmax))])
        n2 = len(xcells2[np.where((xcells2>=-xmax) & (xcells2<=xmax))])
        frac1,frac2 = dfunc.calc_fraction(n1,n2)

        Lambda1_frac.append(frac1)
        Lambda2_frac.append(frac2)
    
    plt.plot(time_steps, Lambda1_frac,color = dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(time_steps, Lambda2_frac,color = dfunc.colour_Lambda(Lambdas[1]))
    return ratio

def cells_change_fraction(data_dir, width=40,norm=False,t_star=False,par2=False):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir,channels=True, width=width)
    Lambdas = sorted(Lambdas)
    frac1,frac2 = dfunc.calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))

    initial_frac = frac2[0]

    time_steps = np.asarray(time_steps)-3.0
    frac1 = abs((frac1 - frac1[0])/(1-frac1[0]))
    frac2 = abs((frac2 - initial_frac)/(1-initial_frac))

    ratio = Lambdas[1]/Lambdas[0]

    par = [1,1]
    norm_timesteps = time_steps
    if norm:
        par,err = dfunc.find_tanh_1param(time_steps[30:],frac2[30:])
        if par2:
            try:
                par,err = dfunc.find_tanh_2param(time_steps[30:],frac2[30:],par[0])
                print(f"f0={initial_frac}: ({round(par[1],2)}+-{round(err[1],2)})*tanh(t/({round(par[0],2)}+-{round(err[0],2)}))")
            except:
                print(f"Could not find values for f0={initial_frac}, ratio={ratio}")
                return 0
            frac2 = frac2/par[1]
        else:
            print(f"f0={initial_frac}: tanh(t/({round(par[0],2)}+-{round(err[0],2)}))")
        norm_timesteps = time_steps/par[0]
    
    #plt.plot(time_steps[30:],frac1[30:],color=dfunc.colour_Lambda(Lambdas[0]))
    if not t_star:
        match round(initial_frac,2):
            case 0.50:
                plt.scatter(norm_timesteps,frac2,marker="x",s=5,color=dfunc.colour_ratio(ratio),alpha=0.5,linewidth=1)
            case 0.33:
                plt.scatter(norm_timesteps,frac2,marker="o",s=5,color=dfunc.colour_ratio(ratio),alpha=0.5,linewidth=1)
            case 0.2:
                plt.scatter(norm_timesteps,frac2,marker="*",s=5,color=dfunc.colour_ratio(ratio),alpha=0.5,linewidth=1)
            case 0.1:
                plt.scatter(norm_timesteps,frac2,marker="v",s=5,color=dfunc.colour_ratio(ratio),alpha=0.5,linewidth=1)
            case _:
                print(initial_frac)
                plt.plot(norm_timesteps,frac2,color=dfunc.colour_ratio(ratio),alpha=0.5)
        return ratio,norm_timesteps
    if t_star:
        return ratio,par,initial_frac

def surface_frac_analytical_param(data_dirs,c):
    raw1 = []
    raw2 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, frac_interface, frac_external, Lambdas = dfunc.surface_fraction_high(data_dir)

        raw1.append(frac_interface)
        raw2.append(frac_external)
        Lam.append(Lambdas)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = np.asarray(time_steps)
            max_len = len(time_steps)
    Lam = np.unique(Lambdas)
    fraction1 = []
    fraction2 = []

    if len(Lambdas) != 1:

        # Pad the shorter sequences with NaN up to max_len
        for item1, item2 in zip(raw1, raw2):
            padded1 = list(item1) + [np.nan] * (max_len - len(item1))
            fraction1.append(padded1)

            padded2 = list(item2) + [np.nan] * (max_len - len(item2))
            fraction2.append(padded2)

        fraction1 = np.asarray(fraction1)
        fraction2 = np.asarray(fraction2)

        # Number of valid datapoints at each position
        count1 = np.sum(~np.isnan(fraction1), axis=0)

        # Only retain positions with >= 2 datapoints
        mask1 = count1 >= 2

        avg1 = np.nanmean(fraction1[:, mask1], axis=0)
        std1 = np.nanstd(fraction1[:, mask1], axis=0)/np.sqrt(count1[mask1])

        plt.errorbar(time[mask1], avg1,yerr=std1,fmt="x", color=c)

    save = False
    if save:
        data = np.asarray([avg1,std1])
        header = "avg\tstdev"
        np.savetxt("Inner_fraction_Lambda"+str(Lam[1]), np.transpose(data), delimiter="\t", fmt="%g",header=header)
    return [Lam[1]]







#Plotting utilities
def get_file_paths(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    return files

def repeat_files(filepath,n_repeats):
    repeats = []
    for i in range(0,n_repeats):
        repeats.append(filepath+"\\repeat"+str(i))
    return repeats

def get_exit_data(data_dir):
    exit_file = data_dir + "\\exited_cells.dat"
    df = pd.read_csv(exit_file, sep="\t")
    return df


"""test_file_single= "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\FreeGrow\\Lambda1\\repeat4"
test_file_double= "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\Interacting_colonies\\stress\\Lambda1AND10\\repeat4"
test_issue = repeat_files("C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\FreeGrow\\Lambda5",5)
test_exit = "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\InfiniteChannel\\channels_test\\repeat7"

plt.figure(figsize=(5, 3.5))

plot_leaving_cells(test_exit)

plt.show()"""