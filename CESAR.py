import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

def ode_Cas12a(c, t, k): 
    """primitive kinetic model for Cas12a detection of target DNA
    concentrations in nM and time in minutes """

    Cas12a = c[0]
    crRNA = c[1]
    targetDNA = c[2]
    ssDNA = c[3]
    Cas12a_crRNA = c[4]
    Cas12a_crRNA_targetDNA = c[5]
    # rate constants   
    k_cr = k[0]
    k_combine = k[1]
    # enzymatic rate constants
    k_cat = k[2]
    K_M = k[3]

    Cas12adot =  - k_cr * Cas12a * crRNA
    crRNAdot = Cas12adot #- k_cr * Cas12a * crRNA
    Cas12a_crRNAdot = Cas12adot - k_combine * targetDNA * Cas12a_crRNA #k_cr * Cas12a * crRNA - k_combine * targetDNA * Cas12a_crRNA
    targetDNAdot = - k_combine * targetDNA * Cas12a_crRNA 
    Cas12a_crRNA_targetDNAdot = k_combine * targetDNA * Cas12a_crRNA
    ssDNAdot = k_cat * Cas12a_crRNA_targetDNA * ssDNA/(K_M + ssDNA) 
    return (Cas12adot,crRNAdot,targetDNAdot,ssDNAdot,Cas12a_crRNAdot,Cas12a_crRNA_targetDNAdot)
    

'''
Define constants here. The reaction was done at 20 nM Cas12a and 250 nM crRNA
'''
k_cr = 1.0
k_combine = 0.00183333333
k_cat = 10
K_M = 500.0
k0=(k_cr, k_combine, k_cat, K_M) 
plt.style.use('ggplot')

color_code = ['r','tomato','peachpuff',
              'darkgreen','limegreen','lawngreen',
              'indigo','b','lightskyblue']
line_style=['-','.']
fig,ax = plt.subplots(2,1)

raw = [418.34127712642345,555.56757658751,313.7268813598616,0.5]
for i in range(len(raw)):
    conc = raw[i]
    # initial condition (Cas12adot,crRNAdot,targetDNAdot,ssDNAdot,Cas12a_crRNAdot,Cas12a_crRNA_targetDNAdot)
    c0=(20,250,conc,200,0.0,0.0)
    time = np.linspace(0,120,10000)
    f1 = lambda y,t: ode_Cas12a(y, t, k0)
    sim = spi.odeint(f1,c0,time)
    # print('result of sim is ',sim)

    # Dictionary for plotting afterwards
    cdict={}
    cdict['Cas12a']=sim[:,0]
    cdict['crRNA']=sim[:,1]
    ssDNA = sim[:,3]
    cdict['ssDNA']=ssDNA
    cdict['Cas12a_crRNA']=sim[:,4]
    cdict['Cas12a_crRNA_targetDNA']=sim[:,5]

    datadict={'time': time, 'fluorescence': cdict}
    current_fluorescence = []
    for j in range(len(ssDNA)):
        current_fluorescence.append(100*(1-(ssDNA[j]-ssDNA[-1])/ssDNA[0]))  
    ax[1].plot(time, current_fluorescence,label='[targetDNA] = %.3f nM'%(conc),color=color_code[i],linewidth=2)
    ax[0].plot(time, ssDNA,label='[targetDNA] = %.3f nM'%(conc),color=color_code[i],linewidth=2)
    ax[0].legend()
    ax[1].legend()
    fluorescence = 100*(cdict['ssDNA'][0]-cdict['ssDNA'][-1])/cdict['ssDNA'][0]
        
ax[0].set_xlabel('time/min')
ax[0].set_ylabel('ssDNA remained/nM')
ax[0].set(title="ssDNA-time prediction given different target DNA concentration")
ax[1].set_xlabel('time/min')
ax[1].set_ylabel('fluorescence/%')
ax[1].set(title="fluorescence-time prediction given different target DNA concentration")
ax[0].grid(False)
ax[1].grid(False)

plt.tight_layout()
plt.savefig('./CESAR.png',dpi=500)
plt.show()
