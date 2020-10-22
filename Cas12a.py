import scipy.integrate as spi
import numpy as n
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

color_code = ['r','tomato','peachpuff',
              'darkgreen','limegreen','lawngreen',
              'indigo','b','lightskyblue']
line_style=['-','.']
#result=[]

a = 1
b = - 0.2
k_rpa_lib=[a * math.exp(b / t) for t in range(37,42)]

k_rpa_lib = [0.4, 0.6, 0.8]
def ode_Cas12a(c, t, k, k_rpa): 
    """primitive kinetic model for Cas12a detection of target DNA
    and collateral ssDNase activity.
    concentrations in nM
    time in minutes """
    if c[2]>1000:
        c[2]=1000
#	if t<10:
#		k_rpa=0
#	if 10<t<30:
#		k_rpa=k_rpa

    # single species
    Cas12a = c[0]; crRNA = c[1]; targetDNA = c[2]; collateralDNA = c[3];
    # complex species
    Cas12a_crRNA = c[4]; Cas12a_crRNA_targetDNA = c[5];
#    k_123=k_rpa

    # rate constants   
    k_c12cr = k[0]; k_c12crtarget = k[1]; 
    # enzymatic rate constants
    k_c = k[2]; K_M = k[3];
    
    Cas12adot =  - k_c12cr*Cas12a*crRNA
    crRNAdot = - k_c12cr*Cas12a*crRNA
    
    Cas12a_crRNAdot = k_c12cr*Cas12a*crRNA - k_c12crtarget * targetDNA * Cas12a_crRNA
    targetDNAdot = - k_c12crtarget * targetDNA * Cas12a_crRNA+k_rpa*targetDNA
    
    Cas12a_crRNA_targetDNAdot = k_c12crtarget * targetDNA * Cas12a_crRNA
    
    
    collateralDNAdot = - k_c * Cas12a_crRNA_targetDNA * collateralDNA/(K_M + collateralDNA) 
    
    return (Cas12adot,crRNAdot,targetDNAdot,collateralDNAdot,Cas12a_crRNAdot,
            Cas12a_crRNA_targetDNAdot)
    
    
    
fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax2 = fig.add_subplot(212)
ax1 = plt.subplot(111)
#ax2 = plt.subplot(212)


# Cas12a
'''

Define constants here. These come from plate reader experiments of Cas12a and 
afterwards fitting function to this. The reaction was done at 20 nM Cas12a and 
20 nM crRNA

'''
#k_c12cr = 1.0; k_c12crtarget = 0.1; k_c = 5; K_M = 20.0;
#k_c12cr = 0.000408423865267452; k_c12crtarget = 2782559.40220713; k_c = 1; K_M = 7.74263682681127
k_c12cr = 1.0; k_c12crtarget = 0.00183333333; k_c = 10; K_M = 1500.0
k0=(k_c12cr, k_c12crtarget, k_c, K_M) 

####
concentrations=[]
results=[]
for i in range(2,5):
#	k_rpa=k_rpa_lib[abc]
    for abc in range(len(k_rpa_lib)):
        
        conc=1000*(10**(-i*1))
        c0=(1,10,conc,200.0,0.0,0.0)
        #k_rpa=k_rpa_lib[abc]
        k_rpa=k_rpa_lib[abc]
        # print('rpatx is ,',k_rpa)
        time = np.linspace(0,120,10000)
        f1 = lambda y,t: ode_Cas12a(y, t, k0, k_rpa)
        sim = spi.odeint(f1,c0,time)
        # print('result of sim is ',sim)
    # Dictionary for plotting afterwards
        cdict={}
        cdict['Cas12a']=sim[:,0]
        cdict['crRNA']=sim[:,1]
        collateralDNA = sim[:,3]
        cdict['collateralDNA']=collateralDNA
        cdict['Cas12a_crRNA']=sim[:,4]
        cdict['Cas12a_crRNA_targetDNA']=sim[:,5]

        datadict={'time': time, 'fluorescence': cdict}


#### Labeling for target DNA concentration loop
        tmp = 3*i+abc-6
        print('%d+%d=%d'%(i,abc,tmp))
        ax1.plot(time, collateralDNA,label='[targetDNA] = '+str(conc)+' nM\nRPA temperature = %d °C' %(37+abc),color=color_code[(tmp)],linewidth=2)
        print(color_code[tmp])
####

        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
        # print(cdict['collateralDNA'])
        fluorescence=100*(cdict['collateralDNA'][0]-cdict['collateralDNA'][-1])/cdict['collateralDNA'][0]
    	#print(fluorescence)
        results.append(fluorescence)
        concentrations.append(conc)
        
ax1.set_xlabel('time (min)')
ax1.set_ylabel('ssDNA remained (nM)')
ax1.set(title="ssDNA-time prediction plot at different RPA temperature and target DNA concentration")
ax1.grid('on')

plt.show()
