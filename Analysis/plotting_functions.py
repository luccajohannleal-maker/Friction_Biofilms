import matplotlib.pyplot as plt
import numpy as np
import drag_functions as dfunc
import utilities as ut
import pandas as pd
import os
import re
import glob

#plots quantity vs time
def plot_count(data_dir):
    time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
    
    for Lambda in Lambdas:
        if Lambda == 1.0 and len(Lambda1_counts) == len(time_steps):
            plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(Lambda))

        if Lambda != 1.0 and len(Lambdanot1_counts) == len(time_steps):
            plt.plot(time_steps, Lambdanot1_counts, 'o', color=dfunc.colour_Lambda(Lambda))
    return Lambdas

def plot_rg_linear(data_dir):
    time_steps, Rg_tot, Lambdas, Rg1, Rg2 = dfunc.RgLambda_time(data_dir)

    if len(Lambdas) != 1:
        Rg_combined = np.asarray([Rg1,Rg2])
        for i in range(0,len(Lambdas)):
            params,err = dfunc.RgLambda_time_est(Rg_combined[:,i])
            print(f"Behaviour for $\Lambda={Lambdas[Lambdas != 1.0][0]}$: $log_2(R_g) = ({round(params[0],3)}+-{round(err[0],3)}) + t/({round(params[1],3)}+-{round(err[1],3)})")
            if len(Rg_combined[:,i]) == len(time_steps):
                plt.plot(time_steps, np.log2(Rg_combined[:,i]), '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[i]))
                plt.plot(time_steps, dfunc.Rg_Lambda_growth(np.array(time_steps), *params), '-', color=dfunc.colour_Lambda(Lambdas[i]))  
 
    else:
        params_tot,err_tot = dfunc.RgLambda_time_est(Rg_tot)
        print(f"Behaviour for $\Lambda={Lambdas[0]}$: $log_2(R_g) =({round(params_tot[0],3)}+-{round(err_tot[0],3)}) + t/({round(params_tot[1],3)}+-{round(err_tot[1],3)})")
        plt.plot(time_steps, np.log2(Rg_tot), '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[0]))#, label=" Counts $\Lambda = 1$")
        plt.plot(time_steps, dfunc.Rg_Lambda_growth(np.array(time_steps), *params_tot), '-', color=dfunc.colour_Lambda(Lambdas[0]))
    return Lambdas

def plot_GR_linear(data_dir):
    time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
    
    if len(Lambda1_counts)>0:
        params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(Lambda1_counts, time_step=0.1)
        print(f"Behaviour for $\Lambda=1$: $log_2 N(t) =({round(params_Lambda1[1],3)}+-{round(err_Lambda1[1],3)}) + t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
        if len(Lambda1_counts) == len(time_steps):
            plt.plot(time_steps, np.log2(Lambda1_counts), '.', markersize=6, color=dfunc.colour_Lambda(1.0))#, label=" Counts $\Lambda = 1$")
            plt.plot(time_steps, dfunc.doubling_linear_growth(np.array(time_steps), *params_Lambda1), '-', color=dfunc.colour_Lambda(1.0))#, label="Fit: $log_2(N)="+str(round(params_Lambda1[0], 2)) +" + t/"+str(round(params_Lambda1[1], 2)) +"$")
            
    if len(Lambdanot1_counts)>0:
        params_Lambdanot1, err_Lambdanot1 = dfunc.estimate_growth_rate(Lambdanot1_counts, time_step=0.1)
        print(f"Behaviour for $\Lambda={Lambdas[Lambdas != 1.0][0]}$: $log_2 N(t) =({round(params_Lambdanot1[1],3)}+-{round(err_Lambdanot1[1],3)}) + t/({round(params_Lambdanot1[0],3)}+-{round(err_Lambdanot1[0],3)})")
        if len(Lambdanot1_counts) == len(time_steps):
            plt.plot(time_steps, np.log2(Lambdanot1_counts), '.', markersize=6, color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]))#, label="$Counts \Lambda =" + str(Lambdas[Lambdas != 1.0][0])+"$")
            plt.plot(time_steps, dfunc.doubling_linear_growth(np.array(time_steps), *params_Lambdanot1), '-', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]))#, label="Fit: $log_2(N)="+str(round(params_Lambdanot1[0], 2)) +" + t/"+str(round(params_Lambdanot1[1], 2)) +"$")
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


def plot_COM_interacting(data_dir):
    files = get_file_paths(data_dir)

    time_steps = []

    xcom1 = [] 
    xcom2 = []
    xcomtot=[]


    for file_path in files:
        match = re.search(r'(\d+)', os.path.basename(file_path))
        if not match:
            continue
        time_step = int(match.group(1))
        df = pd.read_csv(file_path, sep="\t")
        time_steps.append(time_step * 0.1)
        Lambdas = sorted(dfunc.find_Lambdas(df))

        xtot,ytot,ztot= dfunc.centerBiofilm(df)
        xcomtot.append(abs(xtot))

        if len(Lambdas) != 1:
            df1=dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
            x1,y1,z1= dfunc.centerBiofilm(df1)
            xcom1.append(abs(x1))

            df2=dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
            x2,y1,z1= dfunc.centerBiofilm(df2)
            xcom2.append(abs(x2))
            
    if len(Lambdas) == 1:
        plt.plot(time_steps,xcomtot, ".", color=dfunc.colour_Lambda(Lambdas[0]))
    else:
        #plt.plot(time_steps,xcomtot, "--", color="k")
        plt.plot(time_steps,xcom1, color=dfunc.colour_Lambda(Lambdas[0]))
        plt.plot(time_steps,xcom2,color=dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas


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

        for i in range(0,len(stress[0,:])-1):
            stress[1,i] = par_stress[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            stress[2,i] = perp_stress[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            stress[3,i] = shear1[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            stress[4,i] = shear2[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()


        c_Lambda = dfunc.colour_Lambda(Lambda)
        plt.plot(stress[0,:], stress[2,:], "x",color = c_Lambda)
        plt.plot(stress[0,:], stress[1,:], "o",color = c_Lambda)
        plt.plot(stress[0,:], stress[3,:], "*",color = c_Lambda)
        plt.plot(stress[0,:], stress[4,:], "v",color = c_Lambda)
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

        for i in range(0,len(stress[0,:])-1):
            stress[1,i] = p[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()
            stress[2,i] = alpha[np.where((d_cells>=stress[0,i]) & (d_cells<=stress[0,i+1]))].mean()

        c_Lambda = dfunc.colour_Lambda(Lambda)
        plt.plot(stress[0,:], stress[1,:], "x",color = c_Lambda)
        plt.plot(stress[0,:], stress[2,:], "o",color = c_Lambda)

    return Lambdas



#plots averages
def plot_average_growth(data_dirs):
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
        if np.any(Lam == 1):
            padded_L1 = list(item1) + [np.nan] * (max_len - len(item1))
            L1_count.append(padded_L1)
        
        # Pad Lambda not 1 counts
        elif np.any(Lam != 1):
            padded_Ln1 = list(itemnot1) + [np.nan] * (max_len - len(itemnot1))
            Ln1_count.append(padded_Ln1)
        
    L1_count = np.asarray(L1_count)
    Ln1_count = np.asarray(Ln1_count)
    
    # 4. Calculate the averages ignoring the NaNs
    for Lambda in Lam:
        if Lambda == 1:
            avg_L1 = np.nanmean(L1_count, axis=0)
            std_L1 =np.nanstd(L1_count, axis=0)/(avg_L1*np.log(2)*np.sqrt(L1_count.shape[0]))
            params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(avg_L1, time_step=0.1)
            print(f"Behaviour for $\Lambda={Lambda}$: $log_2 N(t) =({round(params_Lambda1[1],3)}+-{round(err_Lambda1[1],3)}) + t/({round(params_Lambda1[0],3)}+-{round(err_Lambda1[0],3)})")
            if len(avg_L1) == len(time):
                plt.errorbar(time, np.log2(avg_L1),fmt= '.',yerr=std_L1, markersize=6, color=dfunc.colour_Lambda(Lambda))
                plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambda1), '-', color=dfunc.colour_Lambda(Lambda))
            
        else:
            avg_Ln1 = np.nanmean(Ln1_count, axis=0)
            std_Ln1 = np.nanstd(Ln1_count, axis=0)/(avg_Ln1*np.log(2)*np.sqrt(Ln1_count.shape[0]))
            params_Lambdanot1, err_Lambdanot1 = dfunc.estimate_growth_rate(avg_Ln1, time_step=0.1)
            print(f"Behaviour for $\Lambda={Lambda}$: $log_2 N(t) =({round(params_Lambdanot1[1],3)}+-{round(err_Lambdanot1[1],3)}) + t/({round(params_Lambdanot1[0],3)}+-{round(err_Lambdanot1[0],3)})")
            if len(avg_Ln1) == len(time):
                plt.errorbar(time, np.log2(avg_Ln1),fmt= '.',yerr=std_Ln1, markersize=6, color=dfunc.colour_Lambda(Lambda))
                plt.plot(time, dfunc.doubling_linear_growth(np.array(time), *params_Lambdanot1), '-', color=dfunc.colour_Lambda(Lambda))
    return Lambdas

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

        xcom1 = [] 
        xcom2 = []
        xcomtot=[]

        for file_path in files:
            match = re.search(r'(\d+)', os.path.basename(file_path))
            if not match:
                continue
            time_step = int(match.group(1))
            df = pd.read_csv(file_path, sep="\t")
            time_steps.append(time_step * 0.1)
            Lambdas = sorted(dfunc.find_Lambdas(df))

            xtot,ytot,ztot= dfunc.centerBiofilm(df)
            xcomtot.append(abs(xtot))

            if len(Lambdas) != 1:
                df1=dfunc.find_Lambda_cells(df,Lambda=Lambdas[0])
                x1,y1,z1= dfunc.centerBiofilm(df1)
                xcom1.append(abs(x1))

                df2=dfunc.find_Lambda_cells(df,Lambda=Lambdas[1])
                x2,y1,z1= dfunc.centerBiofilm(df2)
                xcom2.append(abs(x2))

        rawtot.append(xcomtot)
        Lam.append(Lambdas)
        if len(Lambdas) != 1:
            raw1.append(xcom1)
            raw2.append(xcom2)

        # Track the longest time sequence
        if len(time_steps) > max_len:
            time = time_steps
            max_len = len(time_steps)
    Lam = np.unique(Lam)
    com1 = []
    com2 = []
    comtot = []

    if len(Lambdas) != 1:

        # 2. Pad the shorter sequences with NaN up to max_len
        for item1, item2, itemtot in zip(raw1, raw2, rawtot):
            padded1 = list(item1) + [np.nan] * (max_len - len(item1))
            com1.append(padded1)
            
            padded2 = list(item2) + [np.nan] * (max_len - len(item2))
            com2.append(padded2)

            paddedtot = list(itemtot) + [np.nan] * (max_len - len(itemtot))
            comtot.append(paddedtot)
        
        comtot = np.asarray(comtot)
        avgtot = np.nanmean(comtot, axis=0)
        stdtot =np.nanstd(comtot, axis=0)/(np.sqrt(comtot.shape[0]))
        #plt.errorbar(time, avgtot,yerr=stdtot,fmt="--", color="k")

        com1 = np.asarray(com1)
        avg1 = np.nanmean(com1, axis=0)
        std1 =np.nanstd(com1, axis=0)/(np.sqrt(com1.shape[0]))
        plt.errorbar(time, avg1,yerr=std1, color=dfunc.colour_Lambda(Lam[0]))

        com2 = np.asarray(com2)
        avg2 = np.nanmean(com2, axis=0)
        std2 =np.nanstd(com2, axis=0)/(np.sqrt(comtot.shape[0]))
        plt.errorbar(time, avg2,yerr=std2, color=dfunc.colour_Lambda(Lam[1]))
        
    else:
        for item in rawtot:
            paddedtot = list(item) + [np.nan] * (max_len - len(item))
            comtot.append(paddedtot)

        comtot = np.asarray(comtot)

        avgtot = np.nanmean(comtot, axis=0)
        stdtot =np.nanstd(comtot, axis=0)/(np.sqrt(comtot.shape[0]))
        plt.errorbar(time, avgtot,yerr=stdtot, color=dfunc.colour_Lambda(Lam[0]))
            
    return Lambdas

#fancier plots
def t_doub_Lambda(data_dirs):

    tdoubl1=[]
    tdoubl2=[]

    for data_dir in data_dirs:
        time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
        if len(Lambda1_counts)>0:
            params_Lambda1, err_Lambda1 = dfunc.estimate_growth_rate(Lambda1_counts, time_step=0.1)
            tdoubl1.append(params_Lambda1[0])

                 
        if len(Lambdanot1_counts)>0:
            Lambda = Lambdas[Lambdas != 1.0][0]
            params_Lambdanot1, err_Lambdanot1 = dfunc.estimate_growth_rate(Lambdanot1_counts, time_step=0.1)
            tdoubl2.append(params_Lambdanot1[0])
    
    if len(tdoubl1)>0:
        plt.errorbar(1,np.mean(tdoubl1),yerr=np.std(tdoubl1)/np.sqrt(len(tdoubl1)),fmt="o", color=dfunc.colour_Lambda(1))
        plt.scatter(np.ones(len(tdoubl1)),tdoubl1,color=dfunc.colour_Lambda(1),alpha=0.4,s=1)
                 
    if len(tdoubl2)>0:
        plt.errorbar(Lambda,np.mean(tdoubl2),yerr=np.std(tdoubl2)/np.sqrt(len(tdoubl2)),fmt="o", color=dfunc.colour_Lambda(Lambda))
        plt.scatter(np.ones(len(tdoubl2))*Lambda,tdoubl2,color=dfunc.colour_Lambda(Lambda),alpha=0.4,s=2)
    
    return Lambdas


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




plt.figure(figsize=(5, 3.5))
test_file_single= "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\FreeGrow\\Lambda1\\repeat4"

test_file_double= "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\Interacting_colonies\\stress\\Lambda1AND10\\repeat4"

test_issue = repeat_files("C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\FreeGrow\\Lambda5",5)
plot_COM_interacting(test_file_double)

plt.legend()
plt.show()