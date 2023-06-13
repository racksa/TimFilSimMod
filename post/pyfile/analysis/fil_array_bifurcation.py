import numpy as np
import matplotlib.pyplot as plt





fig = plt.figure(dpi=2000)
ax = fig.add_subplot(1,1,1)
ax.axvline(12, c='black', ls='dashed')
ax.axvline(20, c='black', ls='dashed')
ax.axvline(29.54, c='black', ls='dashed')
ax.annotate('Chaotic motion', (0.5,0.1), rotation=90)
ax.annotate('Beating with phase', (12.5,0.1), rotation=90)
ax.annotate('Synchronised beating', (20.5,0.1), rotation=90)
ax.annotate('Synchronised whirling', (30.04,0.1), rotation=90)

ax.set_xlim(0, 46)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$H/L$')
ax.axes.get_yaxis().set_visible(False)
# ax.set_ylabel(r'$1$')
plt.savefig('fig/fil_array_bifurcation.png', format='png')
plt.savefig('fig/fil_array_bifurcation.pdf', format='pdf')
plt.show()