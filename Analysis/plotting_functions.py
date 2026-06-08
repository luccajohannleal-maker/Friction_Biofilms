import matplotlib.pyplot as plt
import numpy as np
import drag_functions as dfunc
import pandas as pd
import os
import re
import glob

def plot_count(data_dir):
    time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
    
    for Lambda in Lambdas:
        if Lambda == 1.0 and len(Lambda1_counts) == len(time_steps):
            plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(Lambda), label="$\Lambda = 1$")

        if Lambda != 1.0 and len(Lambdanot1_counts) == len(time_steps):
            plt.plot(time_steps, Lambdanot1_counts, 'o', color=dfunc.colour_Lambda(Lambda), label="$\Lambda = "+str(Lambda)+"$")

def plot_rg(data_dir):
    time_steps, Rg_tot, Lambdas, Rg_Lambda1, Rg_Lambdanot1 = dfunc.RgLambda_time(data_dir)
    if len(Lambdas) != 1:
        if len(Rg_Lambdanot1) == len(time_steps):
            plt.plot(time_steps, Rg_Lambdanot1, 'o', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]), label="$\Lambda = " + str(Lambdas[Lambdas != 1.0][0])+"$")    

        if len(Rg_Lambda1) == len(time_steps):
            plt.plot(time_steps, Rg_Lambda1, 'o', color=dfunc.colour_Lambda(1.0), label="$ \Lambda = 1$")

    else:
            plt.plot(time_steps, Rg_tot, 'o', color=dfunc.colour_Lambda(Lambdas[0]), label="$ \Lambda = " + str(Lambdas[0])+"$")

def plot_GR(data_dir):
    time_steps, Lambda1_counts, Lambdanot1_counts, Lambdas = dfunc.counts(data_dir)
    
    if len(Lambda1_counts)>0:
        growth_rate_Lambda1,err_Lambda1 = dfunc.estimate_growth_rate(Lambda1_counts, time_step=0.1)
        print(f"Estimated doubling time for Lambda=1: {growth_rate_Lambda1}+-{err_Lambda1}")
        if len(Lambda1_counts) == len(time_steps):
            plt.plot(time_steps, Lambda1_counts, 'o', color=dfunc.colour_Lambda(1.0))#, label=" Counts $\Lambda = 1$")
            plt.plot(time_steps, 2**(np.array(time_steps)/growth_rate_Lambda1), '-', color=dfunc.colour_Lambda(1.0), label="Fit $\Lambda = 1$: $N(t)= 2^{(t/"+str(round(growth_rate_Lambda1, 2)) +")}$")
            
    if len(Lambdanot1_counts)>0:
        growth_rate_Lambdanot1, err_Lambdanot1 = dfunc.estimate_growth_rate(Lambdanot1_counts, time_step=0.1)
        print(f"Estimated doubling time for Lambda={Lambdas[Lambdas != 1.0][0]}: {growth_rate_Lambdanot1}+-{err_Lambdanot1}")
        if len(Lambdanot1_counts) == len(time_steps):
            plt.plot(time_steps, Lambdanot1_counts, 'o', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]))#, label="$Counts \Lambda =" + str(Lambdas[Lambdas != 1.0][0])+"$")
            plt.plot(time_steps, 2**(np.array(time_steps)/growth_rate_Lambdanot1), '-', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]), label="Fit $ \Lambda = " + str(Lambdas[Lambdas != 1.0][0])+"$" + ": $ N(t)= 2^{(t/"+str(round(growth_rate_Lambdanot1, 2)) +")}$")

def plot_shape_anis_time(data_dir):

    """
    Relative shape anisotropy in this case is defined:
    For eig1, eig2 being the 2 eigenvalues of the gyration tensor,
    diff = abs( (eig1-eig2)/(eig1+eig2) )

    """
    time_steps, Rg_tot_eig, Lambdas, Rg_Lambda1_eig, Rg_Lambdanot1_eig = dfunc.Gyr_eig_time(data_dir)


    shape_tot1 = []

    if len(Lambdas) == 1:
        for i in range(0,len(time_steps)):
            if Rg_tot_eig[i,0] and Rg_tot_eig[i,1] == 0: #one cell
                shape_tot1.append(1) 
            else:
                shape_tot1.append(1.5*(Rg_tot_eig[i,0]**4+Rg_tot_eig[i,1]**4)/((Rg_tot_eig[i,0]**2+Rg_tot_eig[i,1]**2)**2) - 0.5)
        plt.plot(time_steps,shape_tot1, 'o', color=dfunc.colour_Lambda(Lambdas[0]), label="$ \Lambda = " + str(Lambdas[0])+"$")


    else: #If only one Lambda, only Rg_tot is required, so skip the rest
        shape_tot2 = []
        for Lambda in Lambdas:
            if Lambda == 1.0:
                for i in range(0,len(time_steps)):
                    if Rg_Lambda1_eig[i,0] + Rg_Lambda1_eig[i,1] == 0: #one cell
                        shape_tot1.append(1) 
                    else:
                        shape_tot1.append(1.5*(Rg_Lambda1_eig[i,0]**4+Rg_Lambda1_eig[i,1]**4)/((Rg_Lambda1_eig[i,0]**2+Rg_Lambda1_eig[i,1]**2)**2) - 0.5)
                plt.plot(time_steps,shape_tot1, 'o', color=dfunc.colour_Lambda(1.0), label="$ \Lambda = 1$")
  
            else:
                for i in range(0,len(time_steps)):
                    if Rg_Lambdanot1_eig[i,0] + Rg_Lambdanot1_eig[i,1] == 0: #one cell
                        shape_tot2.append(1) 
                    else:
                        shape_tot2.append(1.5*(Rg_Lambdanot1_eig[i,0]**4+Rg_Lambdanot1_eig[i,1]**4)/((Rg_Lambdanot1_eig[i,0]**2+Rg_Lambdanot1_eig[i,1]**2)**2) - 0.5)
                plt.plot(time_steps,shape_tot2, 'o', color=dfunc.colour_Lambda(Lambdas[Lambdas != 1.0][0]), label="$\Lambda = " + str(Lambdas[Lambdas != 1.0][0])+"$")
                    


def plot_stress_time(data_dir):
    file_pattern = os.path.join(data_dir, "biofilm_*.dat")
    files = sorted(glob.glob(file_pattern), key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    time_steps = []

    Lambda1_stress_perp = []
    Lambda1_stress_par = []
    Lambda1_stress_shear1 = []
    Lambda1_stress_shear2 = []

    Lambdanot1_stress_shear1 = []
    Lambdanot1_stress_shear2 = []
    Lambdanot1_stress_par = []
    Lambdanot1_stress_perp = []


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
                Lambda1_stress_par.append(abs(par.mean()))
                Lambda1_stress_shear1.append(abs(tau1.mean()))
                Lambda1_stress_shear2.append(abs(tau2.mean()))

            else:

                cells = dfunc.find_Lambda_cells(df,Lambda=Lambda)
                par,perp,tau1,tau2 = dfunc.find_stress(cells)
                Lambdanot1_stress_perp.append(abs(perp.mean()))
                Lambdanot1_stress_par.append(abs(par.mean()))
                Lambdanot1_stress_shear1.append(abs(tau1.mean()))
                Lambdanot1_stress_shear2.append(abs(tau2.mean()))

        
    for Lambda in Lambdas:
        if Lambda == 1.0:
            plt.plot(time_steps,Lambda1_stress_perp, "x",label= "$\sigma_{\perp}$", color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambda1_stress_par, "o",label= "$\sigma_{\parallel}$",color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambda1_stress_shear1, "*",label= "$shear 1$", color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambda1_stress_shear2, ".",label= "$shear 2$", color=dfunc.colour_Lambda(Lambda))
        else:
            plt.plot(time_steps,Lambdanot1_stress_perp, "x",label= "$\sigma_{\perp}$", color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambdanot1_stress_par, "o",label= "$\sigma_{\parallel}$",color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambdanot1_stress_shear1, "*",label= "$shear 1$", color=dfunc.colour_Lambda(Lambda))
            plt.plot(time_steps,Lambdanot1_stress_shear2, ".",label= "$shear 2$", color=dfunc.colour_Lambda(Lambda))

    
    



#Plotting utilities
def repeat_files(filepath,n_repeats):
    repeats = []
    for i in range(0,n_repeats):
        repeats.append(filepath+"\\repeat"+str(i))
    return repeats


plt.figure(figsize=(5, 3.5))
test_file_single= "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\FreeGrow\\Lambda1\\repeat4"

test_file_double= "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\data\\Interacting_colonies\\Lambda1AND5\\repeat0"

test_issue = "C:\\Users\\lucca\\Desktop\\GeneratedOutput\\SimOutput\\test\\Stress\\repeat1"

plot_stress_time(test_issue)
plt.legend()
plt.show()
