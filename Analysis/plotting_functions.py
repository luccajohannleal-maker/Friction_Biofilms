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
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)
    Lambdas = sorted(Lambdas)

    plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(Lambdas[0]))

    if len(Lambdas) != 1.0:
        plt.plot(time_steps, Lambda2_counts, 'o', color=dfunc.colour_Lambda(Lambdas[1]))
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

        time_steps = []

        com1 = [] 
        com2 = []
        t_collision = dfunc.colonies_collided(files) #finds time when colonies collide

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

    n1 = np.geomspace(40,colony1[0,:].max(),1000)
    par,err = dfunc.find_scaling_law(colony1[0,:],colony1[1,:])
    fit1 = dfunc.linear_exp(n1,*par)
    print(f"Lambda = {Lam[0]}, R =N^({par}+-{err})")
    plt.scatter(colony1[0,:],colony1[1,:],marker=".",color = dfunc.colour_Lambda(Lam[0]),alpha=0.5)
    plt.plot(n1,fit1,color = dfunc.colour_Lambda(Lam[0]))

    n2 = np.geomspace(40,colony1[0,:].max(),1000)
    par,err = dfunc.find_scaling_law(colony2[0,:],colony2[1,:])
    fit2 = dfunc.linear_exp(n2,*par)
    print(f"Lambda = {Lam[1]}, R =N^({par}+-{err})")
    plt.scatter(colony2[0,:],colony2[1,:],marker=".",color = dfunc.colour_Lambda(Lam[1]))
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
    
    plt.errorbar(Lambdas[0],np.mean(tdoubl1),yerr=np.std(tdoubl1)/np.sqrt(len(tdoubl1)),fmt="o", color=dfunc.colour_Lambda(Lambdas[0]))
    plt.scatter(np.ones(len(tdoubl1))*Lambdas[0],tdoubl1,color=dfunc.colour_Lambda(Lambdas[0]),alpha=0.4,s=1)
                 
    if len(Lambdas)!=1:
        plt.errorbar(Lambdas[1],np.mean(tdoubl2),yerr=np.std(tdoubl2)/np.sqrt(len(tdoubl2)),fmt="o", color=dfunc.colour_Lambda(Lambdas[1]))
        plt.scatter(np.ones(len(tdoubl2))*Lambdas[1],tdoubl2,color=dfunc.colour_Lambda(Lambdas[1]),alpha=0.4,s=2)
    
    return Lambdas

def plot_yfraction(filepath):
    df = pd.read_csv(filepath, sep="\t")
    Lambdas = sorted(dfunc.find_Lambdas(df))

    df = dfunc.centerCells(df)

    ycells = np.asarray(df["pos_y"])
    
    df1 = dfunc.find_Lambda_cells(df,Lambdas[0])
    df2 = dfunc.find_Lambda_cells(df,Lambdas[1])
    ycells1 = np.asarray(df1["pos_y"])
    ycells2 = np.asarray(df2["pos_y"])
    

    ymax = max(ycells)
    ymin = min(ycells)

    n_points = 20
    fraction = np.zeros((3,n_points))
    fraction[0,:] = np.linspace(ymin,ymax,n_points)

    for i in range(1,len(fraction[0,:])):
        n1 = len(ycells1[np.where((ycells1>=fraction[0,i-1]) & (ycells1<=fraction[0,i]))])
        n2 = len(ycells2[np.where((ycells2>=fraction[0,i-1]) & (ycells2<=fraction[0,i]))])
        frac1,frac2 = dfunc.calc_fraction(n1,n2)
        fraction[1,i] = frac1
        fraction[2,i] = frac2
    
    plt.plot(fraction[0,:], fraction[1,:],color = dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(fraction[0,:], fraction[2,:],color = dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def yfraction_repeats(data_dirs):
    filepaths = []
    Lambdas = []
    ymax=0
    ymin=0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        files = get_file_paths(data_dir)
        filepaths.append(files[-1])

        df = pd.read_csv(files[-1], sep="\t")
        ycells = np.asarray(df["pos_y"])
        Lambdas = sorted(dfunc.find_Lambdas(df))
        if max(ycells) > ymax:
            ymax = max(ycells)
        if min(ycells) <ymin:
            ymin = min(ycells)

    n_points = 20
    fraction = np.zeros((3,n_points))
    fraction[0,:] = np.linspace(ymin,ymax,n_points)
    n1=np.zeros((n_points))
    n2=np.zeros((n_points))
    for filepath in filepaths:
        df = pd.read_csv(filepath, sep="\t")
        df1 = dfunc.find_Lambda_cells(df,Lambdas[0])
        df2 = dfunc.find_Lambda_cells(df,Lambdas[1])
        ycells1 = np.asarray(df1["pos_y"])
        ycells2 = np.asarray(df2["pos_y"])

        for i in range(1,len(fraction[0,:])):
            n1[i] += len(ycells1[np.where((ycells1>=fraction[0,i-1]) & (ycells1<=fraction[0,i]))])
            n2[i] += len(ycells2[np.where((ycells2>=fraction[0,i-1]) & (ycells2<=fraction[0,i]))])
    frac1,frac2 = dfunc.calc_fraction(n1,n2)
    fraction[1,:] = frac1
    fraction[2,:] = frac2       
    plt.plot(fraction[0,:], fraction[1,:],"--",color = dfunc.colour_Lambda(Lambdas[1]))
    plt.plot(fraction[0,:], fraction[2,:],color = dfunc.colour_Lambda(Lambdas[1]))
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

def plot_fraction_time(data_dir):
    time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)
    Lambdas = sorted(Lambdas)
    frac1,frac2 = dfunc.calc_fraction(np.asarray(Lambda1_counts), np.asarray(Lambda2_counts))

    plt.plot(time_steps,frac1,color=dfunc.colour_Lambda(Lambdas[0]))
    plt.plot(time_steps,frac2,color=dfunc.colour_Lambda(Lambdas[1]))
    return Lambdas

def plot_average_fraction(data_dirs):
    raw_L1 = []
    raw_L2 = []
    Lam = []
    time = []
    max_len = 0

    # 1. Collect raw data and find the longest timeline
    for data_dir in data_dirs:
        time_steps, Lambda1_counts, Lambda2_counts, Lambdas = dfunc.counts(data_dir)
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

plot_COM_interacting(test_file_double)

plt.show()"""
