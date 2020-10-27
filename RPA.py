import numpy as np
import matplotlib.pyplot as plt
import math

time = np.linspace(0,20,1000)


def amplification(template,p):
    global concentration
    concentration = template * (1+p) ** time
    print(concentration[-1])

end_concentration = []
amplification(0.5,0.4)
end_concentration.append(float(concentration[-1]))
plt.style.use('ggplot')
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(time,concentration,label='temperature at 37°C',linewidth=2,color='r')
ax1.set_xlabel('time/min')
ax1.set_ylabel('[RPA produc]/nM')
ax1.set(title="in sillico RPA plot")
ax1.grid(False)

amplification(0.5,0.42)
end_concentration.append(float(concentration[-1]))
ax1.plot(time,concentration,label='temperature at 39°C',linewidth=2,color='tomato')

amplification(0.5,0.38)
end_concentration.append(float(concentration[-1]))
ax1.plot(time,concentration,label='temperature at 42°C',linewidth=2,color='peachpuff')
plt.tight_layout()
ax1.legend()
ax1.set_xticks([0,5,10,15,20])
plt.show()

with open('./end_concentration.txt','w') as f:
    for i in end_concentration:
        f.write(str(i)+'\n')

