
# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Read data from file
Mjup = 318.
data = pd.read_csv('Jupiter_MRrel_CMF0_logFeH_0.dat', sep='\s+',header=0)
M_CMF0_logFeH0_Tint107 = data['M_tot[M_E]']
R_CMF0_logFeH0_Tint107 = data['R_tot[R_J]']
# Mass-radius plot
xmin = 0.04
xmax = 1.50
ymin = 0.78
ymax = 1.05
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', which='major', labelsize=14)
plt.title("GASTLI, solar envelope composition")
plt.plot(M_CMF0_logFeH0_Tint107/Mjup, R_CMF0_logFeH0_Tint107, color='deepskyblue',linestyle='solid',\
     linewidth=4, label=r'CMF = 0')
# Jupiter ref
plt.plot([1.], [1.], 'X', color='darkorange',label=r'Jupiter',markersize=14,\
     markeredgecolor='black')
plt.xlim((xmin,xmax))
plt.ylim((ymin,ymax))
plt.xlabel(r'Mass [$M_{Jup}$]',fontsize=14)
plt.ylabel(r'Radius [$R_{Jup}$]',fontsize=14)
plt.legend()
fig.savefig('Jupiter_MRrel.pdf',bbox_inches='tight',format='pdf', dpi=1000)
plt.close(fig)
