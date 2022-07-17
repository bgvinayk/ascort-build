#ASCORT Calculations - Vinay
from reaktoro import *
import numpy as np
import pandas as pd
# import pickle
import streamlit as st
# import csv
from scipy.optimize import fsolve

##### Decision Tree functions------
def corr(Sc_1,Sc_3,pH,mu_flow,evap_rate,a_op,a_oph):
    a_corr='No Product'
    Avg_corr='_'
    if Sc_1==0:
        a=st.radio('Is heavy metal allowed? (Zn/Mo)',['Yes','No'])
        if a=='No':
            a_corr='No Product Available'
        else:
            if Sc_3=='No':
                b=st.radio('Is phosphate in product allowed?',['Yes','No'])
                if b=='Yes':
                    if pH<8.5:
                        a_corr='Impackt 501C, Impackt 502C, Impackt 504C'
                        Avg_corr='10, 10, 10'
                        a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a2=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a3=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        Avg_dose_corr=str(round(a1,1))+', '+str(round(a2,1))+', '+str(round(a3,1))


                        
                    else:
                        a_corr='Impackt 501C, Impackt 502C '
                        Avg_corr='10,10'
                        a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a2=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        Avg_dose_corr=str(round(a1,1))+','+str(round(a2,1))
                        
                else:
                    a_corr='No Product Available'

            else:
                c=st.radio('Is phosphate in product allowed?',['Yes','No'])
                if c=='No':
                    a_corr='No Product Available'
                else:
                    if pH<8.5:
                        a_corr='Impackt 504C, Impackt 103C'
                        Avg_corr='10, 5'
                        a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a2=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        Avg_dose_corr=str(round(a1,1))+', '+str(round(a2,1))
                        
                    else:
                        a_corr='Impackt 502C, Impackt 103C'
                        Avg_corr='10, 5'
                        a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a2=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        Avg_dose_corr=str(round(a1,1))+', '+str(round(a2,1))

    else:
        d=st.radio('Is heavy metal allowed? (Zn/Mo)',['Yes','No'])
        if d=='No':
            a_corr='No Product Available'
        else:
            if Sc_3=='No':
                
                    f=st.radio('Is phosphate in product allowed?',['Yes','No'])
                    if f=='Yes':
                        if pH<8.5:
                            a_corr='Impackt 504C, Impackt 501C, Impackt 100C'
                            Avg_corr='10, 10, 5'
                            a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                            a2=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                            a3=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                            Avg_dose_corr=str(round(a1,1))+', '+str(round(a2,1))+', '+str(round(a3,1))



                        else:
                            a_corr='Impackt 501C, Impackt 502C, Impackt 100C'
                            Avg_corr='10, 10, 5'
                            a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                            a2=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                            a3=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                            Avg_dose_corr=str(round(a1,1))+', '+str(round(a2,1))+', '+str(round(a3,1))


                    else:
                        a_corr='No Product Available'
            else:
                g=st.radio('Is phosphate in product allowed?',['Yes','No'])
                if g=='Yes':
                    if pH<8.5:
                        a_corr='Impackt 504C, Impackt 103C, Impackt 100C'
                        Avg_corr='10, 5, 5'
                        a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a2=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a3=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        Avg_dose_corr=str(round(a1,1))+', '+str(round(a2,1))+', '+str(round(a3,1))

                    else:
                        a_corr='Impackt 502C, Impackt 103C ,Impackt 100C'
                        Avg_corr='10, 5, 5'
                        a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a2=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        a3=(5*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        Avg_dose_corr=str(round(a1,1))+', '+str(round(a2,1))+', '+str(round(a3,1))


                else:
                    a_corr='No Product Available'
    return a_corr, Avg_corr, Avg_dose_corr
def scale(LSI,Phosphate_SI,CaSo4_SI,Silica_SI,Total_Fe,mu_flow,evap_rate,a_op,a_oph):
    a_scale='No Product'
    Avg_Scaling='_'
    Avg_dose_Scaling='_'
    if LSI>1:
        if LSI>2.2:
           if Phosphate_SI>100:
              a1=st.radio('Is free halogen content greater than 0.5ppm?',['Yes','No'])   
              if a1=='No':
                a_scale='No Product Available'
              else:
                a2=st.radio('Is Phospahte in product allowed?',['Yes','No'])
                if a2=='Yes':
                    a_scale='Impackt 527S'
                    Avg_Scaling='50'
                    a1=(50*(mu_flow - evap_rate)/1000)*a_oph*a_op
                    Avg_dose_Scaling=str(round(a1,1))
                else:
                    a_scale='No Product Available'
           else:
                z1=st.radio('Is free halogen content greater than 0.5ppm?',['Yes','No'])
                if z1=='Yes':
                    a_scale='No Product Available'
                else:
                    if CaSO4_SI>1:
                        a_scale='No Product Available'
                    else:
                        z2=st.radio('Is Phosphorus in product allowed?',['Yes','No'])
                        if z2=='Yes':
                            z3=st.radio('Is Phosphate in product allowed?',['Yes','No'])
                            if z3=='Yes':
                                if Silica_SI>1:
                                    a_scale='No Product Available'
                                else:
                                    if Total_Fe>3:
                                        a_scale='No Product Available'
                                    else:
                                        a_scale='Impackt 521S'
                            else:
                                a_scale='No Product Available'
                        else:
                            a_scale='No Product Available'
                
        else:
            if Phosphate_SI>100:
                a5=st.radio('Is free halogen content greater than 0.5 ppm?',['Yes','No'])
                if a5=='No':
                    a_scale='No Product Available'
                else:
                    if CaSo4_SI>1:
                       a6=st.radio('Is phosphorus in product allowed?',['Yes','No'])
                       if a6=='No':
                        a_scale='No Product Available'
                       else:
                        if Silica_SI>1:
                            a_scale='No Product Available'
                        else:
                            a_scale='Impackt 527S, Impackt 528S'
                            Avg_Scaling='50, 50'
                            a1=(50*(mu_flow - evap_rate)/1000)*a_oph*a_op
                            a2=(50*(mu_flow - evap_rate)/1000)*a_oph*a_op  
                            Avg_dose_Scaling=str(round(a1,1))+', '+str(round(a2,1))                 
                    else:
                        a7=st.radio('Is Phosphorus in product allowed?',['Yes','No'])
                        if a7=='Yes':
                            a_scale='Impackt 521S'
                            Avg_Scaling='10'
                            a1=(10*(mu_flow - evap_rate)/1000)*a_oph*a_op
                              
                            Avg_dose_Scaling=str(round(a1,1))

                    
                        else:
                            a_scale='No Product Available'
            else:
                a3=st.radio('Is free halogen content greater than 0.5 ppm?',['Yes','No'])
                if a3=='Yes':
                    a_scale='No Product Available'
                else:
                    a4=st.radio('Is Phosphorus in product allowed?',['Yes','No'])
                    if a4=='No':
                        a_scale='No Product Available'
                    else:
                        a_scale='Impackt 527S'
                        Avg_Scaling='50'
                        a1=(50*(mu_flow - evap_rate)/1000)*a_oph*a_op
                        Avg_dose_Scaling=str(round(a1,1))
    else:
        if LSI>2.2:
           a_scale='No Product Available'
        else:
            if Phosphate_SI>100:
                if Total_Fe>3:
                    a_scale='No Product Available'
                else:
                    a_scale='Impackt 527S'
                    Avg_Scaling='50'
                    a1=(50*(mu_flow - evap_rate)/1000)*a_oph*a_op
                    Avg_dose_Scaling=str(round(a1,1))
            else:
                a_scale='No Product Available'
    return a_scale,Avg_Scaling,Avg_dose_Scaling
def multi(Sc_3,LSI,mu_flow,evap_rate,a_op,a_oph):
    a_multi='No Product'
    Avg_Multi='_'
    Avg_dose_multi='_'
    if LSI>1:
        if LSI>2.2:
            a_multi='Impackt 547M'
            Avg_Multi='100'
            a1=(100*(mu_flow - evap_rate)/1000)*a_oph*a_op
            Avg_dose_multi=str(round(a1,1))
        else:
            a=st.radio('Is heavy metal allowed? (Zn/Mo)',['Yes','No'])
            if a=='Yes':
                if Sc_3=='Yes':
                    a_multi='Impackt 541M'
                    Avg_Multi='175'
                    a1=(175*(mu_flow - evap_rate)/1000)*a_oph*a_op
                    Avg_dose_multi=str(round(a1,1))
                else:
                    a_multi='Impackt 540M'
                    Avg_Multi='150'
                    a1=(150*(mu_flow - evap_rate)/1000)*a_oph*a_op
                    Avg_dose_multi=str(round(a1,1))
            else:
                a_multi='Impackt 542M'
                Avg_Multi='100'
                a1=(100*(mu_flow - evap_rate)/1000)*a_oph*a_op
                Avg_dose_multi=str(round(a1,1))
    else:
        a_multi='Impackt 545M'
        Avg_Multi='100'
        a1=(100*(mu_flow - evap_rate)/1000)
        Avg_dose_multi=str(round(a1,1))
    return a_multi,Avg_Multi,Avg_dose_multi
def ox(pH):
    a_ox='No Product'
    a=st.radio('Is there a high oxidant demand from organics?',['Yes','No'])
    if a=='No':
       b=st.radio('Is there a high oxidant demand from ammonia?',['Yes','No'])
       if b=='Yes':
        a_ox='Bulab 6044 and Bulab 6040/Bulab 6041 or Bulab 6038'
       else:
        a_ox='Bulab 6044'
    else:
        c=st.radio('Is there a high oxidant demand from ammonia?',['Yes','No'])
        if c=='Yes':
            a_ox='Bulab 6044 and Bulab 6040/Bulab 6041'
        else:
            if pH<8.5:
                a_ox='Oxamine Program'
            else:
                a_ox='Oxamine Program or Bulab 6044 and Bulab 6040/Bulab 6041 or Bulab 6038'
    return a_ox
def Non(Non_1,pH):
 a_Non='No Product'
 if Non_1=='Algae':
    a1=st.radio('Is it only Algae?',['Yes','No'])
    if a1=='No':
        a2=st.radio('Is Algae a primary concern?',['Yes','No'])
        if a2=='No':
            a3=st.radio('Is copper free product required?',['Yes','No'])
            if a3=='No':
                a_Non='Bulab 6057'
            else:a_Non='Busan 1078CF'
        else:
            a_Non='Bulab 6162'
    else:
        a_Non='Bulab 6024'
 elif Non_1=='Bacteria':
    if pH<8.3:
        a1=st.radio('Is Zinc in use at the site?',['Yes','No'])
        if a1=='Yes':
            a2=st.radio('Is Fungi a concern?',['Yes','No'])
            if a2=='Yes':
                a3=st.radio('Is Foaming a concern?',['Yes','No'])
                if a3=='Yes':
                    a4=st.radio('Is additional efficacy against SRB required?',['Yes','No'])
                    if a4=='Yes':
                        a_Non='Bulab 6158'
                    else:
                        a5=st.radio('Is copper free product required?',['Yes','No'])
                        if a5=='No':
                            a_Non='Bulab 6057'
                        else:
                            a_Non='Busan 1078CF'
                else:
                    a_Non='Bulab 6010'
            else:
                a6=st.radio('Is contact time greater than 4 hours required?',['Yes','No'])
                if a6=='No':
                    a_Non='Bulab 6042'
                else:
                    a7=st.radio('Is additional efficacy against SRB required?',['Yes','No'])
                    if a7=='No':
                        st.write('Is copper free product required?',['Yes','No'])
                        if a7=='No':
                            a_Non='Bulab 6057'
                        else:
                            a_Non='Busan 1078F'
                    else:
                        a_Non='Bulab 6158'
        else:
            a_Non='Bulab 6019'            
    else:
        a8=st.radio('Does blowdown need to meet loacl regulatory requirements for discharge to open waterways?',['Yes','No'])
        if a8=='Yes':
            a_Non='Bulab 6142'
        else:
            a9=st.radio('Is zinc in use at the site?',['Yes','No'])
            if a9=='No':
                a_Non='Bulab 6019'
            else:
                a10=st.radio('Is Fungi a concern?',['Yes','No'])
                if a10=='No':
                    a11=st.write('Is additional efficacy against SRB required?',['Yes','No'])
                    if a11=='No':
                        a12=st.radio('Is copper free product required?',['Yes','No'])
                        if a12=='No':
                            a_Non='Bulab 6057'
                        else:
                            a_Non='Busan 1078CF'
                    else:
                        a_Non='Bulab 6158'
                else:
                    a12=st.radio('Is foaming a concern?',['Yes','No'])
                    if a12=='No':
                        a_Non='Bulab 6010'
                    else:
                        a13=st.radio('Is additional efficacy against SRB required?',['Yes','No'])
                        if a13=='No':
                            a14=st.radio('Is copper free product required?',['Yes','No'])
                            if a14=='No':
                                a_Non='Bulab 6057'
                            else:
                                a_Non='Busan 1078CF'
                        else:
                            a_Non='Bulab 6158'
 else:
    a15=st.radio('Is foaming a concern?',['Yes','No'])
    if a15=='Yes':
        a16=st.radio('Is additional efficacy against SRB required?',['Yes','No'])
        if a16=='Yes':
            a_Non='Bulab 6158'
        else:
            a17=st.radio('Is copper free product required?',['Yes','No'])
            if a17=='No':
                a_Non='Bulab 6057'
            else:
                a_Non='Busan 1078CF'
    else:
        a_Non='Bulab 6010'
    return a_Non    
def disp(Disp_5):
    a_disp='No Product'
    a=st.radio('Does the system contain rubber?',['Yes','No'])
    if a=='Yes':
        b=st.radio('Is the system heavily fouled with organic material?',['Yes','No'])
        if b=='Yes':
            a_disp='Bulab 8065'
        else:
            a_disp='Bulab 8152'
    else:
        c=st.radio('Does the deposit contain both organic and mineral components?',['Yes','No'])
        if c=='yes':
            a_disp='Bulab 8002F'
        else:
            d=st.radio('Is the system heavily fouled with biofilm?',['Yes','No'])
            if d=='No':
                e=st.radio('Are heavy oils present?',['Yes','No'])
                if e=='Yes':
                    a_disp='Bulab 8065'
                else:
                    a_disp='Bulab 8152'
            else:
                if Disp_5=='Yes':
                    f=st.radio('Is aromatic free product required?',['Yes','No'])
                    if f=='Yes':
                        a_disp='Bulab 8012'
                    else:
                        a_disp='Bulab 8148'
                else:
                    a_disp='Bulab 8065'
    return a_disp

#####ASCORT functions---------
def effInh(sProdDos,flag):

  HEDP = sProdDos * (formProducts['PHOS 6 (60% HEDP)'][ind] * 0.6 + formProducts['HEDP (100% HEDP)'][ind])
  Acumer = sProdDos * (formProducts['Acumer 2000 (40%)'][ind] * 0.4 + formProducts['Flosperse DP/KY 2736 (40%)'][ind] * 0.4)
  PBTC = sProdDos * (formProducts['PHOS 9 (50% PBTC)'][ind] * 0.5 + formProducts['PBTC'][ind])
  Na = sProdDos * (formProducts['Sodium Hydroxide (50%)'][ind] * 0.5 * 0.574786 + formProducts['Tetrasodium Ethylene Diamine Tetraacetate (EDTA)'][ind]* 0.241889 + formProducts['Sodium Tolyltriazole'][ind]* 0.148194 
                   + formProducts['Sodium Molybdate (35%)'][ind]*0.35 * 0.223270 + formProducts['Sodium Tolyltriazole, 50%'][ind]*0.5 * 0.148194)
  PCA = sProdDos * (formProducts['BSI 361 (50%)'][ind] * 0.5)
  K = sProdDos * (formProducts['Potassium Hydroxide 45%'][ind] * 0.45 * 0.696869)
  OH = sProdDos * (formProducts['Potassium Hydroxide 45%'][ind] * 0.45 *0.303131 + formProducts['Sodium Hydroxide (50%)'][ind] * 0.5 * 0.425214)
  BSI = sProdDos * (formProducts['BSI 361 (50%)'][ind] * 0.5 + formProducts['BSI 97 (50%)'][ind] * 0.5)
  Disp = sProdDos * (formProducts['Busperse 39'][ind])
  Zn = sProdDos * (formProducts['BL9050 Zinc Chloride 49%'][ind] * 0.49 * 0.479726)
  Cl = sProdDos * (formProducts['BL9050 Zinc Chloride 49%'][ind] * 0.49 * 0.520274)
  SO4 = sProdDos * (formProducts['Sulfuric acid (93%)'][ind] * 0.93 * 0.979446)
  H_ion = sProdDos * (formProducts['Sulfuric acid (93%)'][ind] * 0.93 * 0.020554)
  H3PO4 = sProdDos * (formProducts['Phosphoric acid (75%)'][ind] * 0.75)
  Mo = sProdDos * (formProducts['Sodium Molybdate (35%)'][ind]*0.35 * 0.465967)

  #N = 2.0
  Cai = Ca_ion['mg/L0']
  PO4i = PO4_ion['mg/L0']
  print("SIs: and flag: ", SI_CaCO3, SI_TCP, flag, HEDP, PBTC, PCA)
  if (SI_CaCO3>1.2 and flag==1):

    if (float(HEDP)):
      k_HEDP = np.log(10.568/float(HEDP))/21.116
      I_HEDP = 100/(k_HEDP*(Cai - Cae)/169.56+1)
    else:
      #k_HEDP = 10
      I_HEDP = 0

    if (float(PBTC)):
      k_PBTC = np.log(9.1678/float(PBTC))/16.19
      I_PBTC = 100/(k_PBTC*(Cai - Cae)/169.56+1)
    else:
      #k_PBTC = 10
      I_PBTC = 0

    if (float(PCA)):
      k_PCA = np.log(7.7979/float(PCA))/0.479
      I_PCA = 100/(k_PCA*(Cai - Cae)/169.56+1)
    else:
      #k_PCA = 10
      I_PCA = 0

    iEff = 100*(1-(1-I_HEDP/100)*(1-I_PCA/100)*(1-I_PBTC/100))

  elif (SI_TCP >1000 and flag==2):

    if (float(Acumer)):
      #k_Acumer = np.log(10.31827/float(Acumer))/65422676.1512
      k_Acumer = (2.746/float(Acumer))**(1/0.199)
      # if k_Acumer < 1e-10:
      #   k_Acumer = 1e-10
      I_Acumer = 100/(k_Acumer*float(PO4i - Pae)**4/1.947e-5+1)**0.25
      print("I_Acumer: ", Acumer, k_Acumer, PO4i, Pae, I_Acumer)
    else:
      #k_Acumer = 10
      I_Acumer = 90
    
    iEff = 100*(1-(1-I_Acumer/100))

  else: iEff = 999999

  print("Efficiencies: ",iEff)
  return iEff - targetInh

def pH_est_const(pH, Malk):
  global pH_m, pH_c
  logMalk = np.log10(Malk)
  c17 = 4.5652
  c16 = 1.4301
  dc17 = -0.1212
  dc16 = 0.1642
  x = (pH - c17 - c16*logMalk)/(dc17+dc16*logMalk)
  pH_m = c16 + x * dc16
  pH_c = c17 + x * dc17
def newtonRaphsonDos(x0,e,N,flg):
    print('\n\n*** NEWTON RAPHSON pH METHOD IMPLEMENTATION ***')
    step = 1
    flag = 1
    h= 0.01
    condition = True
    f_x = effInh(x0,flg)
    while condition:
        f_xh = effInh(x0+h,flg)
        print("F_x+h, F_x: ", f_xh, f_x)
        g_x = (f_xh-f_x)/h
        if g_x == 0.0:
            print('Divide by zero error!')
            return 0.1
        x1 = x0 - f_x/g_x
        if (x1 < 0.1):
            return 0.1
        elif (x1 > 1000):
            x1 = 1000
        print('pH Iteration-%d, x1 = %0.6f, x0 = %0.6f' % (step, x1, x0))
        x0 = x1
        step = step + 1
        if step > N:
          flag = 0
          return x0
          #break
        f_x = effInh(x1,flg)
        condition = abs(float(f_x)) > e
    
    if flag==1:
        return x0
    else:
        print('\nNot Convergent.')

def Ion_balance():

    M_alk_molar = 2*M_alk_input /(CaCO3_MW*1000)
    Ca_ion["mg/L0"] = Ca_ion_input
    Mg_ion["mg/L0"] = Mg_ion_input
    Na_ion["mg/L0"] = Na_ion_input        
    K_ion["mg/L0"] = K_ion_input
    Ba_ion["mg/L0"] = Ba_ion_input
    Sr_ion["mg/L0"] = Sr_ion_input 
    Zn_ion["mg/L0"] = Zn_ion_input
    Fe_ion["mg/L0"] = Fe_ion_input
    Al_ion["mg/L0"] = Al_ion_input
    Cl_ion["mg/L0"] = Cl_ion_input 
    F_ion["mg/L0"] = F_ion_input
    SO4_ion["mg/L0"] = SO4_ion_input
    SiO2_ion["mg/L0"] = SiO2_ion_input 
    PO4_ion["mg/L0"] = PO4_ion_input
    NO3_ion["mg/L0"] = NO3_ion_input

    for ions in Input_ions_list: ions["Molar0"] = ions["mg/L0"]/(ions["MW"]*1000)
    
    pH_est = pH_c + pH_m * np.log10(float(M_alk_input))

    H_ion["Molar0"] = 10**(-pH_est)/H_ion["Gamma"]
    H_ion["Act"] = H_ion["Molar0"] * H_ion["Gamma"]
    OH_ion["Act"] = Kw/H_ion["Act"]
    OH_ion["Molar0"] = OH_ion["Act"] / OH_ion["Gamma"]

    PO4_ion["Act"] = PO4_ion["Molar0"]/(1/PO4_ion["Gamma"]+H_ion["Act"]/(K3_H3PO4*HPO4_ion["Gamma"])+H_ion["Act"]**2/(K2_H3PO4*K3_H3PO4*H2PO4_ion["Gamma"])+ H_ion["Act"]**3/K_H3PO4)
    PO4_ion["Molar0"] = PO4_ion["Act"]/PO4_ion["Gamma"]

    HPO4_ion["Act"] = (H_ion["Act"]*PO4_ion["Act"])/K3_H3PO4
    HPO4_ion["Molar0"] = HPO4_ion["Act"]/HPO4_ion["Gamma"]
              
    H2PO4_ion["Act"] = (H_ion["Act"]*HPO4_ion["Act"])/K2_H3PO4
    H2PO4_ion["Molar0"] = H2PO4_ion["Act"]/H2PO4_ion["Gamma"]
              
    H3PO4_ion["Act"] = (H_ion["Act"]*H2PO4_ion["Act"])/K1_H3PO4
    H3PO4_ion["Molar0"] = H3PO4_ion["Act"]/H3PO4_ion["Gamma"]

    HCO3_ion["Molar0"] = (M_alk_molar + H_ion["Molar0"] - OH_ion["Molar0"])/((1+2*K2*HCO3_ion["Gamma"]/(H_ion["Act"]*CO3_ion["Gamma"])))
    HCO3_ion["Act"] = HCO3_ion["Molar0"] * HCO3_ion["Gamma"]
    CO3_ion["Molar0"] = K2*HCO3_ion["Act"]/(H_ion["Act"]*CO3_ion["Gamma"])
    CO3_ion["Act"] = CO3_ion["Molar0"] * CO3_ion["Gamma"]

    charge_sum = 0

    for ai in Ion_bal_ions: charge_sum += ai["Z"]*ai["Molar0"]
          
    if (charge_sum < 0):
        Na_ion["Molar0"] = Na_ion["Molar0"] - charge_sum / Na_ion["Z"]
        Na_ion_ppm = Na_ion["Molar0"] * 1000 * Na_ion["MW"] 
        Cl_ion_ppm = Cl_ion["mg/L0"]
    
    else:
        Cl_ion["Molar0"] = Cl_ion["Molar0"] - charge_sum / Cl_ion["Z"]
        Cl_ion_ppm = Cl_ion["Molar0"] * 1000 * Cl_ion["MW"] 
        Na_ion_ppm = Na_ion["mg/L0"]
        
    return Na_ion_ppm, Cl_ion_ppm
def fcoc(ConcR, pH_control):

  global Acid_ppm, Base_ppm, M_alk_molar, acidalk, M_alk_pH, opt_coc, opt_pH, aprops, HTI
  
  if (pH_control_option!='None'):
    pH_target = pH_control 
  else:
    pH_target = pH_c + pH_m * np.log10(float(M_alk_input*ConcR))


  opt_coc = ConcR
  opt_pH = pH_target
  points_coc = 1 # No. of points for generating plots
  Coc = [ConcR]

  for m in range(points_coc):
    Ca_ion["mg/L0"] = Ca_ion_input * Coc[m]
    Mg_ion["mg/L0"] = Mg_ion_input* Coc[m]
    Na_ion["mg/L0"] = Na_ion_input * Coc[m]       
    K_ion["mg/L0"] = K_ion_input  * Coc[m]
    Ba_ion["mg/L0"] = Ba_ion_input * Coc[m]
    Sr_ion["mg/L0"] = Sr_ion_input * Coc[m]
    Zn_ion["mg/L0"] = Zn_ion_input * Coc[m]
    Fe_ion["mg/L0"] = Fe_ion_input * Coc[m]
    Al_ion["mg/L0"] = Al_ion_input * Coc[m]
    Cl_ion["mg/L0"] = Cl_ion_input * Coc[m]
    F_ion["mg/L0"] = F_ion_input * Coc[m]
    SO4_ion["mg/L0"] = SO4_ion_input * Coc[m]
    SiO2_ion["mg/L0"] = SiO2_ion_input * Coc[m]
    PO4_ion["mg/L0"] = PO4_ion_input * Coc[m]
    NO3_ion["mg/L0"] = NO3_ion_input * Coc[m]
    M_alk = M_alk_input * Coc[m]      #mg/l as CaCO3
    P_alk = P_alk_input * Coc[m]       #mg/l as CaCO3
    # Determine pH

    M_alk_pH = 10**((float(pH_target)-pH_c)/pH_m)
    M_alk_molar = 2*M_alk_pH /(CaCO3_MW*1000)

    for ions in Input_ions_list: 
      ions["Molar0"] = ions["mg/L0"]/(ions["MW"]*1000)

    Acid_ppm = Acid_mmolar = 0
    Base_ppm = Base_mmolar = 0

    if (pH_control_acid_option != 'None'):
      delta_M_alk = M_alk - M_alk_pH
      if (delta_M_alk>0):
        Acid_ppm = delta_M_alk/CaCO3_MW*pH_control_acid[pH_control_acid_option]["MW"]*2/pH_control_acid[pH_control_acid_option]["Conc"]/pH_control_acid[pH_control_acid_option]["H+"]
        Acid_mmolar = Acid_ppm*pH_control_acid[pH_control_acid_option]["Conc"] / pH_control_acid[pH_control_acid_option]["MW"]
        if(pH_control_acid_option == '98% H2SO4'):
          acidalk = '98% H2SO4'
          SO4_ion["Molar0"] += Acid_mmolar/1000

        elif (pH_control_acid_option == '35% HCl'):
          acidalk = '35% HCl'
          Cl_ion["Molar0"] += Acid_mmolar/1000
    
    #print("Acid Requirement: ", Acid_ppm, " ppm")

    if (pH_control_base_option != 'None'):
      delta_M_alk = M_alk_pH - M_alk
      if (delta_M_alk>0):
        Base_ppm = delta_M_alk/CaCO3_MW*pH_control_base[pH_control_base_option]["MW"]*2/pH_control_base[pH_control_base_option]["Conc"]/pH_control_base[pH_control_base_option]["OH-"]
        Base_mmolar = Base_ppm*pH_control_base[pH_control_base_option]["Conc"] / pH_control_base[pH_control_base_option]["MW"]
        if(pH_control_base_option == '47% NaOH'):
          acidalk = '47% NaOH'
          Na_ion["Molar0"] += Base_mmolar/1000
        elif (pH_control_acid_option == '45% KOH'):
          acidalk = '45% KOH'
          K_ion["Molar0"] += Base_mmolar/1000

  #print("Base Requirement: ", Base_ppm, " ppm")
  global SI_TCP, SI_CaCO3, SI_Arag, SI_Gyp, SI_Anh, LSI, RSI, SI_BaSO4, SI_SrCO3, SI_SiO2, SI_MgSi, SI_CaF2, SI_FeCO3, SI_FeOH3, SI_FePO4, SI_Arag, SI_MgOH, SI_MgCO3, SI_MgHPO4, SI_ZnCO3, SI_ZnOH, SI_ZnO, SI_ZnPO4, SI_AlOH
  aprops = catanion(pH_target, M_alk_molar)
#   Ion_S = float(aprops.ionicStrengthStoichiometric())
#   Ion_S = ChemicalProperty.ionicStrength(system)
  Ion_S = ionic_strength(aprops.properties()).val

  for ai in Ion_bal_ions: ai["Gamma"] = np.exp(-A*ai["Z"]**2*(Ion_S**0.5/(1+Ion_S**0.5)-0.3*Ion_S))
  
  pCa = -np.log10(float(Ca_ion["Molar0"])*float(Ca_ion["Gamma"]))
  #pCa = -np.log10(Ca_ion["Molar0"])
  pKsp = -np.log10(Ksp_CaCO3)
  pK2 =  -np.log10(K2)
  Alk_gamma = M_alk_molar*HCO3_ion["Gamma"]
  #Alk_gamma = 2*CO3_ion["Act"]+HCO3_ion["Act"]+OH_ion["Act"]-H_ion["Act"]
  pAlk =  -np.log10(float(Alk_gamma))
  pHs = pK2 - pKsp + pCa + pAlk
  LSI = pH_target - pHs
  RSI = 2*pHs - pH_target


  SI_CaCO3 = aprops.saturationIndex("Calcite")

  SI_TCP = aprops.saturationIndex("Ca3(PO4)2(beta)")

  SI_Gyp = aprops.saturationIndex("Gypsum")
  SI_Anh = aprops.saturationIndex("Anhydrite")
  SI_BaSO4 = aprops.saturationIndex("Barite") 
  SI_SrCO3 = aprops.saturationIndex("Strontianite") 
  SI_SiO2 = aprops.saturationIndex("SiO2(am-gel)")
  SI_MgSi = aprops.saturationIndex("Sepiolite(A)")
  SI_CaF2 = aprops.saturationIndex("Fluorite")/20
  SI_FeCO3 = aprops.saturationIndex("Siderite")
  # SI_FeOH3 = aprops.saturationIndex("Ferrihydrite")
  SI_FePO4 = aprops.saturationIndex("Strengite")
  SI_ZnOH = aprops.saturationIndex("Zn(OH)2(epsilon)")
  SI_ZnPO4 = aprops.saturationIndex("Zn3(PO4)2:4H2O")
  SI_Arag = aprops.saturationIndex("Aragonite")
  SI_MgOH = aprops.saturationIndex("Brucite")
  SI_MgCO3 = aprops.saturationIndex("Magnesite")
  SI_MgHPO4 = aprops.saturationIndex("MgHPO4:3H2O")
  SI_ZnCO3 = aprops.saturationIndex("Smithsonite")
  SI_ZnOH = aprops.saturationIndex("Zn(OH)2(epsilon)")
  SI_ZnO = aprops.saturationIndex("Zincite")
  SI_ZnPO4 = aprops.saturationIndex("Zn3(PO4)2:4H2O")
  SI_AlOH = aprops.saturationIndex("Gibbsite")
  #SI_CaF2 = 0
  #SI_FeCO3 = 0
  SI_FeOH3 = 0
  #SI_FePO4 = 0
  HTI = 0.69314718056*sys_vol*(opt_coc-1)/evap_rate

  Max_SI = max(SI_TCP/TCP_SI_Level, SI_CaCO3/Calcite_SI_Level, SI_Gyp/Gyp_SI_Level, SI_BaSO4/Bari_SI_Level, SI_SrCO3/Stro_SI_Level, SI_SiO2/SiO2_SI_Level,
               SI_MgSi/MgSi_SI_Level, Cl_ion["mg/L0"]/Cl_Level, HTI/HTI_Level)

  print("SIs: ", SI_TCP/TCP_SI_Level, SI_CaCO3/Calcite_SI_Level, SI_Gyp/Gyp_SI_Level, SI_BaSO4/Bari_SI_Level, SI_SrCO3/Stro_SI_Level, SI_SiO2/SiO2_SI_Level,
               SI_MgSi/MgSi_SI_Level, SI_CaF2/CaF2_SI_Level, SI_FeCO3/Side_SI_Level, SI_FeOH3/FeOH_SI_Level, SI_FePO4/Stre_SI_Level, Cl_ion["mg/L0"]/Cl_Level, HTI/HTI_Level, Max_SI)

  return Max_SI - 1
def newtonRaphson(x0, ConcR, e,N):
    print('\n\n*** NEWTON RAPHSON pH METHOD IMPLEMENTATION ***')
    global step_pH
    step_pH = 1
    flag = 1
    h= 0.01
    condition = True
    f_x = fpH(x0, ConcR)
    while condition:
        f_xh  = fpH(x0+h, ConcR)
        print("F_x+h, F_x: ", f_xh, f_x)
        g_x = (f_xh-f_x)/h
        if g_x == 0.0:
            print('Divide by zero error!')
            break
        x1 = x0 - f_x/g_x
        if (x1 < 6):
            x1 = 6
        elif (x1 > 10):
            x1 = 10
        print('pH Iteration-%d, x1 = %0.6f, x0 = %0.6f' % (step_pH, x1, x0))
        x0 = x1
        step_pH = step_pH + 1
        if step_pH > N:
          flag = 0
          return x1
        f_x = fpH(x1, ConcR)
        condition = abs(float(f_x)) > e
    
    if flag==1:
        return x0
    else:
        print('\nNot Convergent.')
def newtonRaphson_coc(x0,pH_target,e,N):
    print('\n\n*** NEWTON RAPHSON COC METHOD IMPLEMENTATION ***')
    step = 1
    flag = 1
    h= 0.01
    condition = True
    fcoc_x = fcoc(x0, pH_target)

    while condition:
        fcoc_xh = fcoc(x0+h, pH_target)
        gcoc_x = (fcoc_xh-fcoc_x)/h
        if gcoc_x == 0.0:
            print('Divide by zero error!')
            break
        x1 = x0 - fcoc_x/gcoc_x
        if (x1 < 1.01):
            x1 = 1.01
        elif (x1 > 30):
            x1 = 30
        print('COC Iteration-%d, x1 = %0.6f, x0 = %0.6f' % (step, x1, x0))
        x0 = x1
        step = step + 1
        if step > N:
            flag = 0
            return x1
        fcoc_x = fcoc(x1, pH_target)
        condition = abs(float(fcoc_x)) > e
    
    if flag==1:
        return x0
    else:
        print('\nNot Convergent.')
def fpH(pH_control, ConcR):

  global corrR
  print("pH_iter , Coc : ", pH_control, ConcR)

  if (Coc_control_type != 'Optimize'):
      SI_tar = fcoc(ConcR, pH_control)
  else :
      pH_target = newtonRaphson_coc(ConcR, pH_control, Tol_err, Max_iter)

  corrR = corr_R(M_alk_pH)

  Max_fpH = max(corrR/Corrosion_target, SI_BaSO4/Bari_SI_Level, SI_Gyp/Gyp_SI_Level, SI_Anh/Anh_SI_Level)

  return (Max_fpH - 1)
def corr_R(M_alk_pH):

    global Dos, smInh_choice, raw_data, x1
    HCO3_i = aprops.speciesMolality("HCO3-")*HCO3_ion["MW"]*1000
    CO3_i = aprops.speciesMolality("CO3-2")*CO3_ion["MW"]*1000
    Ca_i = aprops.speciesMolality("Ca+2")*Ca_ion["MW"]*1000
    Cl_i = aprops.speciesMolality("Cl-")*Cl_ion["MW"]*1000
    F_i = aprops.speciesMolality("F-")*F_ion["MW"]*1000   
    Mg_i = aprops.speciesMolality("Mg+2")*Mg_ion["MW"]*1000 
    Al_i = aprops.speciesMolality("Al+3")*Al_ion["MW"]*1000 
    K_i = aprops.speciesMolality("K+")*K_ion["MW"]*1000    
    PO4_i = aprops.speciesMolality("PO4-3")*PO4_ion["MW"]*1000+aprops.speciesMolality("HPO4-2")*HPO4_ion["MW"]*1000+aprops.speciesMolality("H2PO4-")*H2PO4_ion["MW"]*1000
    Na_i = aprops.speciesMolality("Na+")*Na_ion["MW"]*1000    
    SO4_i = aprops.speciesMolality("SO4-2")*SO4_ion["MW"]*1000
    NO3_i = aprops.speciesMolality("NO3-")*NO3_ion["MW"]*1000 

    cond = HCO3_i*0.715 + Ca_i*2.6 + CO3_i*2.82 + Cl_i*2.14 + F_i*2.86 + Mg_i*3.82 + Al_i*3.44 + K_i*1.84 + PO4_i*2.18 + Na_i*2.13 + SO4_i*1.54 + NO3_i*1.15
    cond_T = cond * 1.02**(Tc - 25)
    
    raw_data = {
        "pH":opt_pH,
        "Ca":float(Ca_ion_input * opt_coc),
        "T_Alk":M_alk_pH,
        "Cl": float(Cl_ion_input * opt_coc),
        "SO4": float(SO4_ion_input * opt_coc),
        "Mg": float(Mg_ion_input * opt_coc),
        "Avg_conductivity": float(cond_T),
        "LSI": LSI,
    }

    print("Avg Cond: ", raw_data["Avg_conductivity"])

    if (inh_type == 'None' or (smInh_choice == 'None' and scInh_choice == 'None')):
      x1 = 6.0863 * raw_data["Cl"]**0.1 * raw_data["SO4"]**0.1 * raw_data["T_Alk"]**0.2 /(raw_data["Ca"]**0.01 * raw_data["Mg"]**0.0 *10**(raw_data["LSI"]*0.11))
      print("x1: ", x1, raw_data["Cl"], raw_data["SO4"], raw_data["T_Alk"], raw_data["Ca"], raw_data["Mg"], raw_data["LSI"])
    else:
      x1 = Corrosion_df[Corrosion_level]
    
    return x1
def catanion(pHest, Malk_molar):
    
#     db = Database("supcrt98.xml")
    editor = ChemicalEditor()
    editor.addAqueousPhaseWithElements("H O Ca Mg Na K Zn Fe C Si P S Cl F N Ba Sr Al")\
    .setChemicalModelDebyeHuckel()

#     solution = AqueousPhase(speciate("H O Ca Mg Na K Zn Fe C Si P S Cl F N Ba Sr Al"))


#     solution = AqueousPhaseWithElements("H O Ca Mg Na K Zn Fe C Si P S Cl F N Ba Sr Al")
#     solution.setActivityModel(ActivityModelDebyeHuckelPHREEQC())
#     system = ChemicalSystem(db)
    system = ChemicalSystem(editor)
#     system = ChemicalSystem(db, solution)

    state = EquilibriumInverseProblem(system)
    state.setTemperature(Tc, "celsius")
    state.setPressure(P, "bar")
    state.pH(pHest, "HCl", "NaOH")
    state.add("H2O", water_bmass, "kg")
#     specs.add("CaCO3", 1, "g")
#     specs.temperature()
#     specs.pressure()
#     specs.pH()

#     solver = EquilibriumSolver(specs)

#     state = ChemicalState(system)
#     state.pressure(P, "atm")
#     state.temperature(Tc, "celsius")
#     state.set("H2O"    , water_bmass, "kg")
    state.add("Na+"    , Na_ion["Molar0"]*water_bmass, "mol")
    state.add("K+"    , K_ion["Molar0"]*water_bmass, "mol")
    state.add("Ca+2"    , Ca_ion["mg/L0"], "g")
    state.add("Mg+2"    , Mg_ion["mg/L0"], "g")
    state.add("Ba+2"    , Ba_ion["mg/L0"], "g")
    state.add("Sr+2"    , Sr_ion["mg/L0"], "g")
    state.add("Fe+2"    , Fe_ion["mg/L0"], "g")
    state.add("Zn+2"    , Zn_ion["mg/L0"], "g")
    state.add("Al+3"    , Al_ion["mg/L0"], "g")
    state.add("CO3-2"    ,  Malk_molar*water_bmass, "mol")
    state.add("H4SiO4"    ,  SiO2_ion["Molar0"]*water_bmass, "mol")
    state.add("SO4-2"   , SO4_ion["Molar0"]*water_bmass, "mol")
    state.add("PO4-3"   , PO4_ion["mg/L0"], "g")
    state.add("Cl-"    , Cl_ion["Molar0"]*water_bmass, "mol")
    state.add("NO3-"    , NO3_ion["mg/L0"], "g")
    state.add("F-"    , F_ion["mg/L0"], "g")
    state = equilibrate(state)
#     conditions = EquilibriumConditions(specs)
#     conditions.temperature(Tc, "celsius")
#     conditions.pressure(P, "atm")
#     conditions.pH(pHest)
#     result = solver.solve(state, conditions)
    return state

def deltax(xe,Ca,PO4,SI):
  xe = float(xe)
  return (-108)*(xe)**5+ (108*Ca + 108*PO4)*(xe)**4- (108*Ca*PO4 +27*float(PO4)**2 + 36*float(Ca)**2)*(xe)**3 + (36*float(Ca)**2*PO4 + 27*Ca*float(PO4)**2 + 4*float(Ca)**3)*(xe)**2- (4*PO4*float(Ca)**3 + 9*float(Ca)**2*float(PO4)**2)*float(xe)+ float(Ca)**3*float(PO4)**2*(1-(1000/SI))

###-------------------------------------MAIN CODE---------
# Sc_1,Non_1,Disp_5,Sc_3,a_rr,a_sys,a_coc,a_rpH,a_max_temp,a_dt,a_op,a_pH_control,a_Na,a_K,a_Ca,a_Mg,a_Ba,a_Sr,a_Zn,a_Fe,a_Al,a_Cl,a_F,a_Alk,a_NO3,a_SO4,a_PO4,a_SiO2,controlled_by,pH_control,a_pH_sp
Sc_1='No'
st.title("Cooling Water Program Recommendation")
st.markdown("## Enter cooling tower parameters")
st.subheader("MECHANICAL CONDITIONS")
a=st.multiselect("What kind of Metallurgy is present?",['Galvanized Steel','Mild Steel','Copper'])
for i in range(len(a)):
    if a[i]=='Copper':Sc_1='Yes'
st.subheader("History of Microbiological Problems")
st.radio('Is Sulfate reducing bacteria present?',['Yes','No'])
Non_1=st.radio('Which microbe is predominant?',['Algae','Bacteria','Fungi'])
# st.subheader('History of Deposit Problems')
# st.radio('Is Calcium Carbonate present?',['Yes','No'])
# st.radio('Is Calcium Phosphate present?',['Yes','No'])
# st.radio('Is Silica present?',['Yes','No'])
# st.radio('Is Calcium Sulfate present?',['Yes','No'])
st.subheader('History of Corrosion problems')
Disp_5=st.radio('Microbially influenced corrosion?',['Yes','No'])
Sc_3=st.radio('Is Pitting Corrosion a concern?',['Yes','No'])
st.subheader('OPERATIONAL CONDITIONS')
a_rr=st.number_input('Recirculation Rate (m3/hr)', value=5000)
#a_rr=st.number_input('Reciculation Rate (m3/hr)',value=5000.0)
a_sys=st.number_input('System Volume (m3)',value=4000)
a_coc=st.number_input('Cycles of Concentration',value=3.0)
arp=st.radio('Is Recirculation pH available?',['Yes', 'No'])##Estimate or Enter ---
if arp=='Yes':
    a_str='Enter'
    a_rpH=st.number_input('',value=8.0)
else:
    a_str='Estimate'
    a_rpH= 8.0
a_max_temp=st.number_input('System Temperature (Celsius)',value=30.0)
a_dt=st.number_input('Delta Temperature (Celsius)',value=10.0)
a_op=st.number_input('Operation time per year (days)',value=365)
a_oph= st.number_input('Operation hours per day',value=24)
pH_control=st.radio('pH Control?',['No','Yes'])
if pH_control=="Yes":
    controlled_by=st.radio('Controlled by',['Sulfuric Acid','Hydrochloric Acid'])
    a_pH_sp=st.number_input('pH setpoint',value=8.0)
else: 
    a_pH_sp = 8.0
    controlled_by='None'
st.subheader('MAKE-UP CHEMICAL CONDITIONS')
st.markdown("### Cations") 
a_Ca=st.number_input('Calcium (ppm)',value=50.0)
a_Mg=st.number_input('Magnesium (ppm)',value=10.0)
a_Ba=st.number_input('Barium (ppm)',value=0.0)
a_Sr=st.number_input('Strontium (ppm)',value=0.0)
a_Na=st.number_input('Sodium (ppm)',value=0.0)
a_K=st.number_input('Potassium (ppm)',value=0.0)
a_Fe=Total_Fe=st.number_input('Iron (ppm)',value=0.0)
a_Al=st.number_input('Aluminium (ppm)',value=0.0)
a_Zn=st.number_input('Zinc (ppm)',value=0.0)

st.markdown("### Anions")
a_Alk=st.number_input('M Alkalinity (ppm as CaCO3)',value=100.0)
a_Cl=st.number_input('Chloride (ppm)',value=50.0)
a_SO4=st.number_input('Sulfate (ppm)',value=20.0)
a_PO4=st.number_input('Total Phosphate (ppm)',value=1.0)
a_F=st.number_input('Fluoride (ppm)',value=0.0)
a_NO3=st.number_input('Nitrate (ppm)',value=0.0)
a_SiO2=st.number_input('Silica (ppm)',value=0.0)

#Initialize
LSI=2.5
pH=8.5
Phosphate_SI=500
CaSO4_SI=0.1
Silica_SI=0

if st.button('Submit'):

    placeholder = st.empty()
    type_ = 'warning'
    msg = 'Please Wait...'
    styledMsg = f'\
        <div class="element-container" style="width: 693px;">\
            <div class="alert alert-{type_} stAlert" style="width: 693px;">\
                <div class="markdown-text-container">\
                    <p>{msg}</p></div></div></div>\
    '
    placeholder.markdown(styledMsg, unsafe_allow_html=True)
    
    global CR_x, formProducts, costProducts
    global pH_c, pH_m, Tc, T, TCP_SI_Level, Calcite_SI_Level_Base, Calcite_SI_Level, Gyp_SI_Level, Anh_SI_Level, Stro_SI_Level, Bari_SI_Level, SiO2_SI_Level, CaF2_SI_Level, Side_SI_Level, FeOH_SI_Level, MgSi_SI_Level, Stre_SI_Level
    CR_x=[50,75,100,125,150]

    #Ions data
    Ca_ion = {'Name': 'Calcium', 'Z': 2, 'MW': 40.078}
    Mg_ion = {'Name': 'Magnesium', 'Z': 2, 'MW': 24.305}
    Na_ion = {'Name': 'Sodium', 'Z': 1, 'MW': 22.989769}
    K_ion = {'Name': 'Potassium', 'Z': 1, 'MW': 39.0983}
    Ba_ion = {'Name': 'Barium', 'Z': 2, 'MW': 137.327}
    Sr_ion = {'Name': 'Sr ion', 'Z': 2, 'MW': 87.62}
    Fe_ion = {'Name': 'Fe ion', 'Z': 2, 'MW': 55.845}
    Zn_ion = {'Name': 'Zinc ion', 'Z': 2, 'MW': 65.38}
    Al_ion = {'Name': 'Al ion', 'Z': 3, 'MW': 26.981539}
    CO3_ion = {'Name': 'Carbonate', 'Z': -2, 'MW': 60.0089}
    HCO3_ion = {'Name': 'Bicarbonate', 'Z': -1, 'MW': 61.0168, 'Gamma': 1}
    SO4_ion = {'Name': 'Sulfate', 'Z': -2, 'MW': 96.06}
    Cl_ion = {'Name': 'Chloride', 'Z': -1, 'MW': 35.453, 'Gamma': 1}
    F_ion = {'Name': 'Fluoride', 'Z': -1, 'MW': 18.99840320}
    NO3_ion = {'Name': 'Nitrate', 'Z': -1, 'MW': 62.0049}
    PO4_ion = {'Name': 'o-Phosphate', 'Z': -3, 'MW': 94.9714}
    HPO4_ion = {'Name': 'HPO4', 'Z': -2, 'MW': 95.979, 'Gamma': 1}
    H2PO4_ion = {'Name': 'H2PO4', 'Z': -1, 'MW': 96.987, 'Gamma': 1}
    H3PO4_ion = {'Name': 'H3PO4', 'Z': 0, 'MW': 97.994, 'Gamma': 1}
    SiO2_ion = {'Name': 'SiO2', 'Z': 0, 'MW': 60.08430, 'Gamma': 1}
    H_ion = {'Name': 'Hydrogen', 'Z': 1, 'MW': 1.00784, 'Gamma': 1}
    OH_ion = {'Name': 'Hydroxide', 'Z': -1, 'MW': 17.008, 'Gamma': 1}
    CaCO3_ion = {'Name': 'CaCO3', 'Z': 0, 'MW': 100.0869, 'Gamma': 1}

    All_ions_list = [Ca_ion, Mg_ion, Na_ion, K_ion, Ba_ion, Sr_ion, Fe_ion, Al_ion, H_ion, OH_ion, Zn_ion, HCO3_ion, CO3_ion, Cl_ion, F_ion, SO4_ion, PO4_ion, HPO4_ion, H2PO4_ion, SiO2_ion, NO3_ion]

    Ion_bal_ions = [Na_ion, K_ion, Ca_ion, Mg_ion, Ba_ion, Sr_ion, Fe_ion, Al_ion, H_ion, OH_ion, Zn_ion, Cl_ion, F_ion, SO4_ion, PO4_ion, NO3_ion, HPO4_ion, H2PO4_ion, HCO3_ion, CO3_ion]

    Input_ions_list = [Na_ion, K_ion, Ca_ion, Mg_ion, Ba_ion, Sr_ion, Zn_ion, Fe_ion, Al_ion, Cl_ion, F_ion, SO4_ion, PO4_ion, NO3_ion, SiO2_ion]

    CaCO3_MW = CaCO3_ion["MW"]

    R = 8.314       #SI units

    Corrosion_df = {'Moderate': 10, 'Good': 5, 'Very Good': 3}
    Corrosion_level = 'Moderate'
    Corrosion_target = Corrosion_df[Corrosion_level]

    costProducts = pd.read_csv('Model_Product_Costs.csv')
    formProducts = pd.read_csv('Model_Product_Formulations.csv')
    # costProducts = pd.read_csv('/app/Model_Product_Costs.csv')
    # formProducts = pd.read_csv('/app/Model_Product_Formulations.csv')

    singleScale = {'521S': {'lim': 150, 'limPO4': 2000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][0])},
                '522S': {'lim': 1.2, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][1])},
                '524S': {'lim': 150, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][2])},
                '527S': {'lim': 150, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][3])},
                '528S': {'lim': 1.2, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][4])},
                'Select Best': {'lim': 150, 'limPO4': 15000, },
                'None': {'lim': 1.2, 'limPO4': 2000, 'Cost': 0}}

    singleCorr = {'501C': {'lim': Corrosion_target, 'Cost': float(costProducts['Cur.Std.Cost/KG'][5])},
                '502C': {'lim': Corrosion_target, 'Cost': float(costProducts['Cur.Std.Cost/KG'][6])},
                '504C': {'lim': Corrosion_target, 'Cost': float(costProducts['Cur.Std.Cost/KG'][7])},
                'Select Best': {'lim': Corrosion_target},
                'None': {'Cost': 0}}

    corrProdMap = {'501C': '541M',
                '502C': '542M',
                '504C': '540M',
                }

    multiFunc = {'540M': {'lim': 150, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][8])},
                '541M': {'lim': 200, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][9])},
                '542M': {'lim': 150, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][10])},
                '543M': {'lim': 200, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][11])},
                '544M': {'lim': 150, 'limPO4': 15000, 'Cost': float(costProducts['Cur.Std.Cost/KG'][12])},
                'Select Best': {'lim': 150, 'limPO4': 15000},
                'None': {'lim': 1.2, 'limPO4': 2000, 'Cost': 0}}

    CaCO3MultiScaleProds = ['540M', '541M', '542M', '543M', '544M']
    TCPMultiScaleProds = ['540M', '541M', '542M', '543M', '544M']

    pH_control_acid = {
            '98% H2SO4': {"MW": 98.079, "Conc": 0.98, "H+": 2, "Cost": 0.5},
            '35% HCl': {"MW": 36.458, "Conc": 0.35, "H+": 1, "Cost": 0.1},
            'None': {"Cost": 0},
    }

    pH_control_base = {
            '47% NaOH': {"MW": 39.997, "Conc": 0.47, "OH-": 1, "Cost": 0.6},
            '45% KOH': {"MW": 56.1056, "Conc": 0.45, "OH-": 1, "Cost": 0.6},
            'None': {"Cost": 0},
    }  

    Tol_err = 0.0001
    Max_iter = 50
    Initial_coc = 3
    Initial_pH = 7.6

    Ion_S = 0.01

    db = Phreeqc("minteq.v4.dat")
    # db = Phreeqc("C:/Docs/Acku Predictor/Practice/minteq.v4.dat")

    # db = open("C:/Docs/Acku Predictor/Practice/minteq.v4.dat", "r")

    water_bmass = 1000 # Taken 1000 kg of water as basis

    P = 1 #atm

    inputs = {}
    inputs['sVol'] = a_sys # System Vol m3
    inputs['rCirc'] = a_rr #Q Recirc rate m3/hr
    inputs['driftRate'] = 0.005
    inputs['leakRate'] = 0.01
    inputs['temp'] = a_max_temp #Q
    inputs['deltaT'] = a_dt #Q deg C
    inputs['coc'] = a_coc #Q base case coc
    inputs['Na'] = a_Na #Q
    inputs['K'] = a_K #Q
    inputs['Ca'] = a_Ca #Q ppm
    inputs['Mg'] = a_Mg #Q ppm
    inputs['Ba'] = a_Ba #Q
    inputs['Sr'] = a_Sr #Q
    inputs['Zn'] = a_Zn #Q
    inputs['Fe'] = a_Fe #Q
    inputs['Al'] = a_Al #Q
    inputs['Cl'] = a_Cl #Q
    inputs['F'] = a_F #Q
    inputs['MAlk'] = a_Alk    #QUESTION 21)- M Alkalinity
    inputs['NO3'] = a_NO3 #Q
    inputs['SO4'] = a_SO4 #Q
    inputs['PO4'] = a_PO4 #Q 
    inputs['SiO2'] = a_SiO2 

    # Setting submitted case conditions #
    inputs['pHBcaseOption'] = a_str   #'Enter' if recirc water pH is entered
    inputs['pHBcase'] = a_rpH #Recirc water pH

    if pH_control=='Yes':
        pH_control='Target'
    else:
        pH_control='None'
    inputs['pHOption'] = pH_control #'Target' if pH setpoint is entered
    inputs['pHTarget'] = a_pH_sp #pH setpoint
    if controlled_by=='Sulfuric Acid':
        controlled_by='98% H2SO4'
    else:controlled_by='35% HCl'
    inputs['pHAcid'] = controlled_by  #pH Acid selection
    inputs['pHBase'] = '47% NaOH'  #pH Base selection
    inputs['cocOption'] = 'Optimize' #'Target' if coc target is selected, 'None' if same coc as base case
    inputs['cocTarget'] = 5 #If coc 'Target' is selected
    inputs['inhType'] = 'Multi'
    inputs['stValue'] = 'None'
    inputs['ctValue'] = 'None'
    inputs['mtValue'] = '542M'

    TCP_SI_Level = TCP_SI_Level_Base = 1000
    Calcite_SI_Level_Base = Calcite_SI_Level = 1.2
    Anh_SI_Level = 1
    Gyp_SI_Level = 1
    Bari_SI_Level = 1
    Stro_SI_Level = 1
    SiO2_SI_Level = 1.1
    MgSi_SI_Level = 1
    CaF2_SI_Level = 1
    Side_SI_Level = 1.2
    FeOH_SI_Level = 1000
    Stre_SI_Level = 1
    Cl_Level = 500
    HEDP_pH_Level = 8.9
    HEDP_Ca_Level = 240
    HTI_Level = 70

    Tc = inputs['temp']     # degree C
    T = Tc + 273.15  #degree K
    coc_base = inputs['coc'] # base coc

    global Na_ion_input, K_ion_input, Ca_ion_input, Mg_ion_input, Ba_ion_input, Sr_ion_input, Zn_ion_input, Fe_ion_input, Al_ion_input, Cl_ion_input, F_ion_input, NO3_ion_input, SO4_ion_input, PO4_ion_input, SiO2_ion_input, M_alk_input, P_alk_input
    # Water chemistry inputs (ppm)
    Na_ion_input = inputs['Na']  
    K_ion_input  = inputs['K']
    Ca_ion_input = inputs['Ca']
    Mg_ion_input = inputs['Mg']
    Ba_ion_input = inputs['Ba']
    Sr_ion_input = inputs['Sr']
    if (inputs['Zn']<1e-5): Zn_ion_input = 1e-5
    else: Zn_ion_input = inputs['Zn']
    Al_ion_input = inputs['Al']
    if (inputs['Fe']<1e-5): Fe_ion_input = 1e-5
    else: Fe_ion_input = inputs['Fe']
    Cl_ion_input = inputs['Cl']
    if (inputs['F']<1e-5): F_ion_input = 1e-5
    else: F_ion_input = inputs['F']
    NO3_ion_input = inputs['NO3']
    SO4_ion_input = inputs['SO4']
    PO4_ion_input = inputs['PO4']
    SiO2_ion_input = inputs['SiO2']
    M_alk_input = inputs['MAlk']        #ppm as CaCO3
    P_alk_input = 0              #ppm as CaCO3

    #System Parameters
    sys_vol = inputs['sVol']        #m3
    recirc_rate = inputs['rCirc']     #m3/hr
    deltaTc = inputs['deltaT']
    deltaTf = deltaTc * 9/5
    evap_rate_constant = 1   #evaporation % per every 10 deg F
    evap_rate = recirc_rate * evap_rate_constant * deltaTf / 1000
    drift_rate_constant = inputs['driftRate']
    drift_rate = recirc_rate * drift_rate_constant/100
    #leak_rate_constant = 0.05
    leak_rate_constant = inputs['leakRate']
    leak_rate = recirc_rate * leak_rate_constant/100

    if inputs['pHBcaseOption']=='Estimate':
        pH_c = 4.4
        pH_m = 1.6
    else:
        pH_est_const(inputs['pHBcase'], M_alk_input*coc_base)
    
    global Kw, K1, K2, K2_H2SO4, K1_H3PO4, K2_H3PO4, K3_H3PO4, K_H3PO4, Ksp_CaCO3, epsilon, A #, Ksp_CaSO4_gy, Ksp_CaSO4_an, Ksp_beta_CaPO4, Ksp_TCP

    K = np.array([[-283.971000000,	-356.309400000,	-107.887100000], [-0.050698420,	-0.060919640,	-0.032528490]
                    , [13323.000000000,	21834.370000000,	5151.790000000], [102.244470000,	126.833900000,	38.925610000]
                    , [-1119669.000000000,	-1684915.000000000,	-563713.900000000]])

    [Kw, K1, K2] = np.power(10,(K[0,:] + K[1,:]*T + K[2,:]/T + K[3,:]*np.log10(T) + K[4,:]/ (T**2)))


    # K2 of H2SO4
    K2_H2SO4 = 10**(577.214 - 246.01 * np.log10(T) - 12717/T + 0.283133 * T - 0.000137566 *T**2)

    # K1, K2, K3 of H3PO4
    K1_H3PO4 = 10**(-799.31/T+4.5535-0.013486*T)
        
    K2_H3PO4 = np.exp(-16.82131+0.013846*Tc-0.000163*Tc**2)

    K3_H3PO4 = 10**(-0.0000397727*Tc**2 + 0.016344697*Tc - 12.7188333333)

    K_H3PO4 = K1_H3PO4*K2_H3PO4*K3_H3PO4


    # Ksp of CaCO3
    Ksp_CaCO3_cf = 0.27
    Ksp_CaCO3 = 10**(-171.9065 +Ksp_CaCO3_cf -0.077993*T+2839.319/T+71.595*np.log10(T))

    epsilon = 87.74 - 0.4008*Tc + 0.0009398*Tc**2 - 0.00000141*Tc**3
    A = 1824000*(epsilon*T)**-1.5

    for ai in Ion_bal_ions: ai["Gamma"] = np.exp(-A*ai["Z"]**2*(Ion_S**0.5/(1+Ion_S**0.5)-0.3*Ion_S))

    Na_ppm, Cl_ppm = Ion_balance()
    Na_ion_input = Na_ppm
    Cl_ion_input = Cl_ppm

    global SI_TCP, SI_CaCO3, SI_Gyp, SI_Anh, corrR, SI_BaSO4, SI_SrCO3, SI_SiO2, SI_MgSi, SI_CaF2, SI_FeCO3, SI_FeOH3, SI_FePO4, SI_Arag, SI_MgOH, SI_MgCO3, SI_MgHPO4, SI_ZnCO3, SI_ZnOH, SI_ZnO, SI_ZnPO4, SI_AlOH
    global pH_control_type, pH_control_acid_option, pH_control_base_option, pH_control_option, Coc_control_type, acidalk, opt_coc, opt_pH, corrR, M_alk_pH, bcInh_choice#, step_pH
    global inh_type, ssInh_choice, scInh_choice, smInh_choice, Model, smInhFiles, ind, indCost, targetInh, Cae

    #Inhibitor Info

    bInhType = 'None'
    bsInh_choice = 'None'
    bcInh_choice = 'None'
    bmInh_choice = 'None'
    bsInh_conc = 0
    bcInh_conc = 0
    bmInh_conc = 0

    sInhType = inputs['inhType']
    ssInh_choice = 'None'
    scInh_choice = 'None'
    smInh_choice = 'None'

    if(sInhType != 'Single'):
        smInh_choice = inputs['mtValue']
        if (smInh_choice == 'Select Best' or sInhType == 'Best'):
            smInh_choice != 'Select Best'
            smInhFiles = ['541M', '543M', '544M'] #Decision Trees Product List - Multifunctional
        else:
            smInhFiles = [smInh_choice]

    if(sInhType != 'Multi'):
        ssInh_choice = inputs['stValue']
        scInh_choice = inputs['ctValue']
        if (scInh_choice == 'Select Best' or sInhType == 'Best'):
            scInh_choice = 'Select Best'
            scInhFiles = ['501C', '502C', '504C'] #Decision Trees Product List - Single Corrosion
        else:
            scInhFiles = [scInh_choice]

        if (ssInh_choice == 'Select Best' or sInhType == 'Best'):
            ssInh_choice = 'Select Best'
            ssInhFiles = ['521S', '522S', '524S', '527S', '528S'] #Decision Trees Product List - Single Corrosion
        else:
            ssInhFiles = [ssInh_choice]

    print("inhType", sInhType, ssInh_choice, scInh_choice, smInh_choice)
    ssInh_conc = 0
    scInh_conc = 0
    smInh_conc = 0

    Model = []

    acidalk = 'None'
    pH_control_type = 'None'
    Coc_control_type = 'None'
    pH_control_option = 'None'
    pH_control_acid_option = 'None'
    pH_control_base_option = 'None'
    # step_pH = 1

    inh_type = bInhType
    pH_base = pH_c + pH_m * np.log10(float(M_alk_input*coc_base))
    fpH(pH_base, coc_base)   #Base case
    bCorrR = corrR
    bCorrC = bcInh_conc
    # bScaleC = bsInh_conc
    # bMultiC = bmInh_conc
    bLSI = LSI
    bRSI = RSI
    bSI_CaCO3 = SI_CaCO3
    bSI_TCP = SI_TCP
    bSI_Gyp = SI_Gyp
    bSI_BaSO4 = SI_BaSO4
    bSI_SrCO3 = SI_SrCO3
    bSI_SiO2 = SI_SiO2
    bSI_MgSi = SI_MgSi
    bSI_Anh = SI_Anh
    bSI_CaF2 = SI_CaF2
    bSI_FeCO3 = SI_FeCO3
    bSI_FeOH3 = SI_FeOH3
    bSI_FePO4 = SI_FePO4
    bSI_Arag = SI_Arag
    bSI_MgOH = SI_MgOH
    bSI_MgCO3 = SI_MgCO3
    bSI_MgHPO4 = SI_MgHPO4
    bSI_ZnCO3 = SI_ZnCO3
    bSI_ZnOH = SI_ZnOH
    bSI_ZnO = SI_ZnO
    bSI_ZnPO4 = SI_ZnPO4
    bSI_AlOH = SI_AlOH
    bMalk = M_alk_pH
    bCa = Ca_ion['mg/L0']
    bBa = Ba_ion['mg/L0']
    bSr = Sr_ion['mg/L0']
    bMg = Mg_ion['mg/L0']
    bZn = Zn_ion['mg/L0']
    bFe = Fe_ion['mg/L0']
    bAl = Al_ion['mg/L0']
    bNa = Na_ion["Molar0"]*Na_ion["MW"]*1000
    bK = K_ion["Molar0"]*K_ion["MW"]*1000
    bSO4 = SO4_ion["Molar0"]*SO4_ion["MW"]*1000
    bPO4 = PO4_ion['mg/L0']
    bCl = Cl_ion["Molar0"]*Cl_ion["MW"]*1000

    #inh_choice = 'HEDP'
    if (smInh_choice != 'None'):
        Calcite_SI_Level = multiFunc[smInh_choice]['lim'] 
        TCP_SI_Level = multiFunc[smInh_choice]['limPO4'] 

    elif (ssInh_choice != 'None'):
        Calcite_SI_Level = singleScale[ssInh_choice]['lim'] 
        TCP_SI_Level = singleScale[ssInh_choice]['limPO4'] 
    
    inh_type = sInhType
    pH_control_option = inputs['pHOption']
    #pH_control_option = 'Yes'
    
    if (pH_control_option != 'None'):
        #pH_control_type = 'Target'
        pH_control_type = inputs['pHOption']
        pH_control_acid_option = inputs['pHAcid'] #Atleast 1 should be selected
        pH_control_base_option = inputs['pHBase']
    
        if (pH_control_type == 'Target'):
            pH_target = inputs['pHTarget']
  
    Coc_control_type = inputs['cocOption']

    
    if (Coc_control_type == 'Target'):
        Coc_target = inputs['cocTarget']
        
    elif (Coc_control_type == 'None'):
        Coc_target = inputs['coc']

    if (pH_control_option == 'None'):

        if (Coc_control_type == 'Optimize'):
            CR_tar = fpH(Initial_pH, Initial_coc)          # Case 5
        else:
            CR_tar = fpH(Initial_pH, Coc_target)                           # Case 4
    
    elif (pH_control_type == 'Optimize'):
        if (Coc_control_type == 'Optimize'):
            pH = pH_c + pH_m * np.log10(float(M_alk_input*Initial_coc))
            opt_pH = newtonRaphson(pH, Initial_coc, Tol_err,Max_iter)           # Case 2
            if (opt_pH > 6.5):
                print("Optimum pH: ", opt_pH, ", Optimum Coc: ", opt_coc)
            else:
                print("Use a corrosion inhibitor or Call case 3")
        else:
            pH = pH_c + pH_m * np.log10(float(M_alk_input*Coc_target))
            opt_pH = newtonRaphson(pH, Coc_target, Tol_err,Max_iter)    # Case 6
            if (corrR > Corrosion_target):
                print("Unable to run at this Coc, try with a corrosion inhibitor")
    
    else:
        if (Coc_control_type == 'Optimize'):
            CR_tar = fpH(pH_target, Initial_coc)          # Case 3
            print("Target pH: ", pH_target, ", Optimum Coc: ", opt_coc)
            if (corrR > Corrosion_target):
                print("Cannot run CT at this Coc and pH, try with a corrosion inhibitor or change input conditions")
        else:
            CR_tar = fpH(pH_target, Coc_target)         # Case 1
            print("pH: ", opt_pH, ", Coc: ", opt_coc)
            if (corrR > Corrosion_target):
                print("Cannot run CT at this Coc and pH, try with a corrosion inhibitor or change input conditions")
  
    scInh_conc = 0
    ssInh_conc = 0
    smInh_conc = 0
    sMalk = M_alk_pH
    sCa = Ca_ion['mg/L0']
    sBa = Ba_ion['mg/L0']
    sSr = Sr_ion['mg/L0']
    sMg = Mg_ion['mg/L0']
    sZn = Zn_ion['mg/L0']
    sFe = Fe_ion['mg/L0']
    sAl = Al_ion['mg/L0']
    sNa = Na_ion["Molar0"]*Na_ion["MW"]*1000
    sK = K_ion["Molar0"]*K_ion["MW"]*1000
    sSO4 = SO4_ion["Molar0"]*SO4_ion["MW"]*1000
    sPO4 = PO4_ion['mg/L0']
    sCl = Cl_ion["Molar0"]*Cl_ion["MW"]*1000

    # print("SIs: ", SI_TCP, SI_CaCO3, SI_Gyp, "CR: ", corrR, "Na: ", Na_ion["mg/L"], "Cl: ", Cl_ion["mg/L"], "SO4: ", SO4_ion["Molar0"]*SO4_ion["MW"]*1000, "K: ", K_ion["mg/L"])
    print("Acid (ppm): ", Acid_ppm, "Base (ppm): ", Base_ppm)

    bd_flow = evap_rate/ (opt_coc-1) - drift_rate - leak_rate
    mu_flow = bd_flow + drift_rate + leak_rate + evap_rate           #m3/hr
    res_time = sys_vol/mu_flow      #hr

    print("Acid (kg/hr): ", Acid_ppm*(mu_flow - evap_rate)/1000, "Base (kg/hr): ", Base_ppm*(mu_flow - evap_rate)/1000)
  
    coeff = 2.5
    x_Cae = ((coeff*float(Ca_ion["mg/L0"]) + float(M_alk_pH)) - ((coeff*float(Ca_ion["mg/L0"]) + float(M_alk_pH))**2 - 4*coeff*(1-1/float(SI_CaCO3))*float(Ca_ion["mg/L0"])*float(M_alk_pH))**0.5)/(2*coeff)
    Cae = Ca_ion["mg/L0"] - x_Cae
    Alke = float(M_alk_pH) - coeff*x_Cae

    x_Pae = fsolve(deltax,0.1,args = (Ca_ion["mg/L0"],PO4_ion["mg/L0"],SI_TCP))
    Pae = PO4_ion["mg/L0"]-2*x_Pae
    print(x_Pae)
    print("Peascort=",Pae)

    targetInh = 90

    if (inh_type != 'Single' and smInh_choice != 'None'):
        x1 = Corrosion_df[Corrosion_level]
        data=pd.DataFrame({'pH':[raw_data["pH"]]*5,'Ca':[raw_data["Ca"]]*5,'T_Alk':[raw_data["T_Alk"]]*5, 'Cl':[raw_data["Cl"]]*5, 'SO4':[raw_data["SO4"]]*5, 'Mg':[raw_data["Mg"]]*5, 'ProductDosage':CR_x,
                        'Avg_Conductivity':[raw_data["Avg_conductivity"]]*5})
        
        smDosg = []
        smDosCost = []

        sProdDos = 10

        ydata=CR_x
        log_y=np.log(ydata)
        if (smInh_choice == 'Select Best'):
            if (SI_CaCO3 > 1.2 and SI_TCP > 1000):
                smInhFiles = list(set(TCPMultiScaleProds).intersection(set(CaCO3MultiScaleProds)))
            elif SI_TCP > 1000:
                smInhFiles = TCPMultiScaleProds
            elif SI_CaCO3 > 1.2:
                smInhFiles = CaCO3MultiScaleProds
        if ((opt_pH > HEDP_pH_Level) and (Ca_ion["mg/L0"] > HEDP_Ca_Level)):
            smInhFiles = list(set(smInhFiles) - set(['544M']))

        for j in smInhFiles:
            # modelFile = 'C:\\Users\\arohanp\\OneDrive - Buckman\\ASCORT_final\\Model_Conda_'+j+'.sav'
            # Model = pickle.load(open(modelFile,'rb'))
            # xdata= Model.predict(data)
            # curve_fit = np.polyfit(xdata,log_y,1)
            # smcDosg = np.exp(curve_fit[1])*np.exp(curve_fit[0]*x1)
            smcDosg = 100
            indCost = list(costProducts['Impackt']).index("Impackt "+j)
            ind = list(formProducts['Impackt']).index("Impackt "+j)
            smsDosg_CaCO3 = newtonRaphsonDos(sProdDos,0.001,100,1)
            smsDosg_TCP = 0#newtonRaphsonDos(sProdDos,0.001,100,2)
            smsDosg = max(smsDosg_CaCO3, smsDosg_TCP)
            print("Corr Dos: ", smcDosg, " Scale Dos: ", smsDosg, smsDosg_CaCO3, smsDosg_TCP)
            smDosg.append(max(smsDosg, smcDosg))  
            smDosCost.append(max(smsDosg, smcDosg)*float(costProducts['Cur.Std.Cost/KG'][indCost]))

        min_val = min(smDosCost)
        min_idxs = [idx for idx, val in enumerate(smDosCost) if val == min_val]
        smInh_choice = smInhFiles[min_idxs[0]]
        smInh_conc = smDosg[min_idxs[0]]

    if (inh_type != 'Multi'):
        if (scInh_choice != 'None'):

            x1 = Corrosion_df[Corrosion_level]
            data=pd.DataFrame({'pH':[raw_data["pH"]]*5,'Ca':[raw_data["Ca"]]*5,'T_Alk':[raw_data["T_Alk"]]*5, 'Cl':[raw_data["Cl"]]*5, 'SO4':[raw_data["SO4"]]*5, 'Mg':[raw_data["Mg"]]*5, 'ProductDosage':CR_x,
                                'Avg_Conductivity':[raw_data["Avg_conductivity"]]*5})
            
            scDosg = []
            scDosCost = []
            scInhDos = 10

            ydata=CR_x
            log_y=np.log(ydata)
            for j in scInhFiles:
                # modelFile = 'C:\\Users\\arohanp\\OneDrive - Buckman\\ASCORT_final\Model_Conda_'+corrProdMap[j]+'.sav'
                # Model = pickle.load(open(modelFile,'rb'))
                # xdata= Model.predict(data)
                # curve_fit = np.polyfit(xdata,log_y,1)
                # sccDosg = np.exp(curve_fit[1])*np.exp(curve_fit[0]*x1)
                sccDosg = 0
                ind = list(formProducts['Impackt']).index("Impackt "+j)
                sind = list(formProducts['Impackt']).index("Impackt "+corrProdMap[j])
                if (j == '502C'):
                    scInhDos = sccDosg * formProducts['Phosphoric acid (75%)'][sind] / formProducts['Phosphoric acid (75%)'][ind]

                elif('501C'):
                    scInhDos = sccDosg * formProducts['BL9050 Zinc Chloride 49%'][sind] / formProducts['BL9050 Zinc Chloride 49%'][ind]

                else:
                    scInhDos = sccDosg * max (formProducts['Phosphoric acid (75%)'][sind] / formProducts['Phosphoric acid (75%)'][ind], formProducts['BL9050 Zinc Chloride 49%'][sind] / formProducts['BL9050 Zinc Chloride 49%'][ind])
                
                sindCost = list(costProducts['Impackt']).index("Impackt "+j)
                
                scDosg.append(scInhDos)
                scDosCost.append(scInhDos*float(costProducts['Cur.Std.Cost/KG'][sindCost]))

            scmin_val = min(scDosCost)
            scmin_idxs = [sidx for sidx, val in enumerate(scDosCost) if val == scmin_val]
            scInh_choice = scInhFiles[scmin_idxs[0]]
            scInh_conc = scDosg[scmin_idxs[0]]

        if (ssInh_choice != 'None'):

            sProdDos = 10
            ssDosg = []
            ssDosCost = []

            if (ssInh_choice == 'Select Best'):
                if (SI_CaCO3 > 1.2 and SI_TCP > 1000):
                    ssInhFiles = list(set(CaCO3SingleScaleProds).intersection(set(TCPSingleScaleProds)))
                elif SI_TCP > 1000:
                    ssInhFiles = TCPSingleScaleProds
                elif SI_CaCO3 > 1.2:
                    ssInhFiles = CaCO3SingleScaleProds

            if ((opt_pH > HEDP_pH_Level) and (Ca_ion["mg/L0"] > HEDP_Ca_Level)):
                ssInhFiles = list(set(ssInhFiles) - set(['521S', '524S']))

            for j in ssInhFiles:
                indCost = list(costProducts['Impackt']).index("Impackt "+j)
                ind = list(formProducts['Impackt']).index("Impackt "+j)
                ssDosg_CaCO3 = newtonRaphsonDos(sProdDos,0.001,100,1)
                ssDosg_TCP = 0#newtonRaphsonDos(sProdDos,0.001,100,2)
                ssInhDosg = max(ssDosg_CaCO3, ssDosg_TCP)
                ssDosg.append(ssInhDosg)
                ssDosCost.append(ssInhDosg*float(costProducts['Cur.Std.Cost/KG'][indCost]))

            ssmin_val = min(ssDosCost)
            ssmin_idxs = [sidx for sidx, val in enumerate(ssDosCost) if val == ssmin_val]
            ssInh_choice = ssInhFiles[ssmin_idxs[0]]
            ssInh_conc = ssDosg[ssmin_idxs[0]]

    if (inh_type == 'Best'):
        
        if (ssmin_val + scmin_val > min_val):
            ssInh_conc = 0
            ssInh_choice = 'None'
            scInh_conc = 0
            scInh_choice = 'None'
        
        else:
            smInh_choice = 'None'
            smInh_conc = 0

    print("BD Flow (m3/hr): ", bd_flow, "MU Flow (m3/hr): ", mu_flow, "Evaporation Rate (m3/hr): ", evap_rate, "Recirc Rate (m3/hr): ", recirc_rate)
    print("Drift rate (m3/hr): ", drift_rate, "Leak rate (m3/hr): ", leak_rate)
    print("HTI: ", HTI)

    op_days = a_op

    bd_flow1 = evap_rate/ (coc_base-1) - drift_rate - leak_rate
    mu_flow1 = bd_flow1 + drift_rate + leak_rate + evap_rate           #m3/hr
    res_time1 = sys_vol/mu_flow1      #hr
    HTI1 = 0.69314718056*sys_vol/(mu_flow1 - evap_rate)

    delta_mu_flow = mu_flow1 - mu_flow
    delta_bd_flow = bd_flow1 - bd_flow
    water_cost = 1  # $/m3
    bd_cost = 0.5
    water_saved_yr = delta_mu_flow * a_oph * op_days
    #water_savings_hr = delta_mu_flow * water_cost - delta_bd_flow * bd_cost
    water_savings_hr = delta_mu_flow * water_cost
    water_savings_yr = water_savings_hr * a_oph * op_days

    print("BD: ",bd_flow1,bd_flow, mu_flow1, mu_flow)
    print("Annual water savings: %s $/yr" % water_savings_yr)
    
    acid_kg_hr = Acid_ppm*(mu_flow - evap_rate)/1000
    acid_cost_hr = acid_kg_hr * pH_control_acid[pH_control_acid_option]['Cost']
    acid_cost_yr = acid_cost_hr * a_oph * op_days
    base_kg_hr = Base_ppm*(mu_flow - evap_rate)/1000
    base_cost_hr = base_kg_hr * pH_control_base[pH_control_base_option]['Cost']
    base_cost_yr = base_cost_hr * a_oph * op_days

    ssInh_g_hr = ssInh_conc *(mu_flow - evap_rate)
    ssInh_cost_hr = ssInh_g_hr * singleScale[ssInh_choice]['Cost']/1000
    ssInh_cost_yr = ssInh_cost_hr * a_oph * op_days
    print ("Scale Inh cost: %s $/yr" % ssInh_cost_yr)

    scInh_g_hr = scInh_conc *(mu_flow - evap_rate)
    scInh_cost_hr = scInh_g_hr * singleCorr[scInh_choice]['Cost']/1000
    scInh_cost_yr = scInh_cost_hr * a_oph * op_days
    print ("Corr Inh cost: %s $/yr" % scInh_cost_yr)

    smInh_g_hr = smInh_conc *(mu_flow - evap_rate)
    smInh_cost_hr = smInh_g_hr * multiFunc[smInh_choice]['Cost']/1000
    smInh_cost_yr = smInh_cost_hr * a_oph * op_days
    print ("Corr Inh cost: %s $/yr" % smInh_cost_yr)

    bsInh_g_hr = bsInh_conc *(mu_flow1 - evap_rate)
    bsInh_cost_hr = bsInh_g_hr * singleScale[bsInh_choice]['Cost']/1000
    bsInh_cost_yr = bsInh_cost_hr * a_oph * op_days
    print ("Base Case Corr Inh cost: %s $/yr" % bsInh_cost_yr)

    bcInh_g_hr = bCorrC *(mu_flow1 - evap_rate)
    bcInh_cost_hr = bcInh_g_hr * singleCorr[bcInh_choice]['Cost']/1000
    bcInh_cost_yr = bcInh_cost_hr * a_oph * op_days
    print ("Base Case Corr Inh cost: %s $/yr" % bcInh_cost_yr)

    bmInh_g_hr = bmInh_conc *(mu_flow - evap_rate)
    bmInh_cost_hr = bmInh_g_hr * multiFunc[bmInh_choice]['Cost']/1000
    bmInh_cost_yr = bmInh_cost_hr * a_oph * op_days
    print ("Corr Inh cost: %s $/yr" % bmInh_cost_yr)

    inh_cost_yr = ssInh_cost_yr + scInh_cost_yr + smInh_cost_yr - bsInh_cost_yr - bcInh_cost_yr - bmInh_cost_yr

    # with open('/ascort_trial/model_results.csv', 'w', newline='') as csvfile: 
    #     csvfile.truncate()
    #     csvwriter = csv.writer(csvfile) 
    #     csvwriter.writerow(['LSI', 'coc', 'pH', 'SI_CaCO3', 'SI_TCP', 'SI_Gyp', 'SI_SiO2','mu_flow','evap_rate','Operation_days','HTI']) 
    #     #['pH','Coc','LSI','Calcite Saturation','TCP Saturation','CaSO4 Saturation','Silica Saturation']
    #     csvwriter.writerow([LSI, opt_coc, opt_pH, SI_CaCO3, SI_TCP, SI_Gyp, SI_SiO2,mu_flow,evap_rate,a_op,HTI])
    
    data=pd.DataFrame({'LSI':[LSI],'coc':[opt_coc],'pH':[opt_pH],'SI_CaCO3':[SI_CaCO3],'SI_TCP':[SI_TCP],'SI_Gyp':[SI_Gyp],'SI_SiO2':[SI_SiO2],'mu_flow':[mu_flow],'evap_rate':[evap_rate],'Operation_days':[a_op],'HTI':[HTI],'Op_hours':[a_oph]})
    data.to_csv('model_results.csv')
    # data.to_csv('/app/model_results.csv')

    if 'key' not in st.session_state:
        st.session_state['key'] = 1
    
    placeholder.empty()
    
    #ASCORT Calculations End - Vinay

    ##------INTERNAL CALCULATIONS--------

if 'key' in st.session_state:
    df = pd.read_csv('model_results.csv')
    # df = pd.read_csv('/app/model_results.csv')
    # ['pH','Coc','LSI','Calcite Saturation','TCP Saturation','CaSO4 Saturation','Silica Saturation']
    LSI=df['LSI'][0]
    pH=df['pH'][0]
    coc=df['coc'][0]
    Phosphate_SI=df['SI_TCP'][0]
    CaSO4_SI=df['SI_Gyp'][0]
    Silica_SI=df['SI_SiO2'][0]
    CaCO3_SI=df['SI_CaCO3'][0]
    mu_flow=df['mu_flow'][0]
    evap_rate=df['evap_rate'][0]
    a_op=df['Operation_days'][0]
    a_oph = df['Op_hours'][0]
    HTI=df['HTI'][0]
    

    c_scale_corr_solid=0
    c_scale=0
    c_corr=0
    c_multi=0
    c_bio_solid=0
    c_ox=0
    c_Non=0
    c_disp_solid=0
    c_disp=0
    
    st.header('Corrosion and Scaling Products')
    scale_corrosion_sl=st.radio('What kind of product do you want',['Solid','Liquid'])
    if scale_corrosion_sl=='Solid':
        c_scale_corr_solid=1
        st.write('This product needs a multi drum dosage system.')
        if LSI>0.5:
            a_corr_solid='Impackt Inhibitor HW'
            Avg_Solid=str(40)
            a1= (40*(mu_flow - evap_rate)/1000)*a_oph*a_op
            Avg_doing_Solid=str(round(a1,1))
        else:
            a_corr_solid='Impackt Inhibitor SW'
            Avg_Solid=str(40)
            a1= (40*(mu_flow - evap_rate)/1000)*a_oph*a_op
            Avg_doing_Solid=str(round(a1,1))

    else:
        drum=st.radio('What kind of dosage system do you have?',['One drum','Multi drum'])
        if drum=='Multi drum':
            scale_required=st.radio('Is a scale inhibitor required?',['Yes','No'])
            if scale_required=='Yes':
                c_scale=1
                (a_scale,Avg_Scaling,Avg_doing_scale)=scale(LSI,Phosphate_SI,CaSO4_SI,Silica_SI,Total_Fe,mu_flow,evap_rate,a_op,a_oph)
                
            corrosion_required=st.radio('Is a corrosion inhibitor required?',['Yes','No'])    
            if corrosion_required=='Yes':
                c_corr=1
                (a_corr,Avg_corr,Avg_doing_corr)=corr(Sc_1,Sc_3,pH,mu_flow,evap_rate,a_op,a_oph)
                
        else:
            c_multi=1
            (a_multi,Avg_multi,Avg_doing_multi)=multi(Sc_3,LSI,mu_flow,evap_rate,a_op,a_oph)
    st.header('Biocide Products')
    ap=st.radio('Do you need a biocide product?',['Yes','No'])
    if ap=='Yes':
        biocide_sl=st.radio('What kind of Product do you want?',['Solid','Liquid'])
        if biocide_sl=='Solid':
            c_bio_solid=1
            if pH<8.5:
                a_bio_solid='Ox1 - PS-Duroklor56'
                Avg_oxbio_solid='_'
                Avg_doing_bioSolid="_"
            else:
                a_bio_solid='Ox2 - PS-Durobrom'
                Avg_oxbio_solid='_'
                Avg_doing_bioSolid="_"
        else:
            oxbiocide_needed=st.radio('Do you need an oxbiocide product?',['Yes','No'])
            if oxbiocide_needed == 'Yes':
                c_ox=1
                a_ox=ox(pH)
                Avg_ox='_'
                Avg_doing_Ox="_"
            else:
                
                nonoxbiocide_needed=st.radio('Do you need a nonoxbiocide product?',['Yes','No'])
                c_Non=1
                if nonoxbiocide_needed=='Yes':
                    a_Non=Non(Non_1,pH)
                    Avg_Non='_'
                    Avg_doing_Non="_"

    st.header('Dispersant Products')
    ap_1=st.radio('Do you need a dispersant product?',['Yes','No'])
    if ap_1=='Yes':
        
        Disp_sl=st.radio('What kind of product do you need?',['Solid','Liquid'])
        if Disp_sl=='Solid':
            c_disp_solid=1
            Ashrae=st.radio('Is Ashrae component add on was selected?',['Yes','No'])
            if Ashrae=='No':
                a_disp_solid='No Product Available'
                Avg_Disp_Solid='_'
                Avg_doing_DispSolid="_"
            else:
                Algae=st.radio('Is Algae an issue?',['Yes','No'])
                if Algae=='Yes':
                    a_disp_solid='Impackt Disp1'
                    Avg_Disp_Solid='_'
                    Avg_doing_DispSolid="_"
                else:
                    a_disp_solid='Impackt Disp2'
                    Avg_Disp_Solid='_'
                    Avg_doing_DispSolid="_"
        else:
            c_disp=1
            a_disp=disp(Disp_5)
            Avg_Disp='_'
            Avg_doing_Disp="_"
            
    if st.button(' Submit '):
    
        st.title('System Parameters')
        data=pd.DataFrame({'Parameters':['pH','Coc','Holding Time Index (hours)','LSI','Calcite Saturation','TCP Saturation','CaSO4 Saturation','Silica Saturation'],'Values':[str(round(pH,1)),str(round(coc,1)),str(round(HTI,1)),str(round(LSI,2)),str(round(CaCO3_SI,1)),str(round(Phosphate_SI,1)),str(round(CaSO4_SI,1)),str(round(Silica_SI,1))]},index=['1','2','3','4','5','6','7','8'])   
        st.table(data)
        #  c_scale_corr_solid=0
        #     c_scale=0
        #     c_corr=0
        #     c_multi=0
        #     c_bio_solid=0
        #     c_ox=0
        #     c_Non=0
        list_product_type=[]
        list_product=[]
        d=0
        Average_Conc=[]
        Average_dosing=[]
        if c_scale_corr_solid==1:
            d=d+1
            list_product_type.append('Scaling/Corrosion Solid')
            list_product.append(a_corr_solid)
            Average_Conc.append(Avg_Solid)
            Average_dosing.append(Avg_doing_Solid)

        if c_scale==1:
            d=d+1
            list_product_type.append('Scaling Product')
            list_product.append(a_scale)
            Average_Conc.append(Avg_Scaling)
            Average_dosing.append(Avg_doing_scale)
        if c_corr==1:
            d=d+1
            list_product_type.append('Corrosion Product')
            list_product.append(a_corr)
            Average_Conc.append(Avg_corr)
            Average_dosing.append(Avg_doing_corr)
        if c_multi==1:
            d=d+1
            list_product_type.append('Multifunctional Product')
            list_product.append(a_multi)
            Average_Conc.append(Avg_multi)
            Average_dosing.append(Avg_doing_multi)
        if c_bio_solid==1:
            d=d+1
            list_product_type.append('Solid Biocide Product')
            list_product.append(a_bio_solid)
            Average_Conc.append(Avg_oxbio_solid)
            Average_dosing.append(Avg_doing_bioSolid)
        if c_ox==1:
            d=d+1
            list_product_type.append('OxBiocide Product')
            list_product.append(a_ox)
            Average_Conc.append(Avg_ox)
            Average_dosing.append(Avg_ox)
        if c_Non==1:
            d=d+1
            list_product_type.append('NonOxBiocide Product')
            list_product.append(c_Non)
            Average_Conc.append(Avg_Non)
            Average_dosing.append(Avg_doing_Non)
        if c_disp_solid==1:
            d=d+1
            list_product_type.append('Dispersant Product')
            list_product.append(a_disp_solid)
            Average_Conc.append(Avg_Disp_Solid)
            Average_dosing.append(Avg_doing_DispSolid)
        if c_disp==1:
            d=d+1
            list_product_type.append('Dispersant Product')
            list_product.append(a_disp)
            Average_Conc.append(Avg_Disp)
            Average_dosing.append(Avg_doing_Disp)
        ind_list=[]
        for i in range(d):
            ind_list.append(i+1)
        data_1=pd.DataFrame({'Product Type':list_product_type,'Product List':list_product,'Average Concentration (ppm)':Average_Conc,'Average Dosing (Kg/year)':Average_dosing},index=ind_list)
        st.title('Treatment Program')
        st.table(data_1)
        print("The value of a_op is {}".format(a_op))
        #Avg. Dosing = (Avg. Concentration * (mu_flow - evap_rate)/1000)*a_oph*No. of days (kg/hr)
        #Slug Dosage = (System_Volume*Avg. Concentration)/1000


