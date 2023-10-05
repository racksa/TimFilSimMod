import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

endp = 46
fig = plt.figure(figsize=(12,3.5))
ax = fig.add_subplot(1,1,1)
fig2 = plt.figure(figsize=(9,0.7))
ax2 = fig2.add_subplot(1,1,1)
fig3 = plt.figure(figsize=(9,0.7))
ax3 = fig3.add_subplot(1,1,1)

# old
xticks  = [0, 12, 20, 29.54, endp]
x_ticks_label = ["H={:.2f}L".format(x) for x in xticks]
style = ['-.', ':', '--', '-']
markers = ['*', 'v', 's', 'o', '^']
markers_single = ['s', 'o', 's', '*', 'v']
width = [2, 2, 2, 3]
labels = ['Chaotic motion', 'Unsync. planar beating', 'Sync. planar beating', 'Whirling']

line_pos = [2.84, 14.2, 22.73, 31.25]

for i, pos in enumerate(line_pos):
	b=0
	if(i%2==0):
		b=-1
	ax.vlines(pos, b, b+1, color='black', ls=style[i], linewidth=1)
ax2.vlines(xticks[3], -1, 1, color='black', ls=style[2], linewidth=1)
ax3.vlines(xticks[1:-1], -1, 1, color='black', ls=style, linewidth=1)

	
for i, pos in enumerate(xticks[:-1]):
	ax.hlines(0, xticks[i], xticks[i+1], color='black', ls=style[i], linewidth=width[i], label = labels[i])
	ax.vlines(pos, -0.05, 0.05, color='black', linewidth=1)
#ax.annotate('Chaotic motion', (0.5,0.1), rotation=90)
#ax.annotate('Beating with phase', (12.5,0.1), rotation=90)
#ax.annotate('Synchronised beating', (20.5,0.1), rotation=90)
#ax.annotate('Synchronised whirling', (30.04,0.1), rotation=90)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_position('center')
ax.spines['bottom'].set_visible(False)
ax.set_xticks(xticks[:-1], x_ticks_label[:-1] )
ax.set_xlim(0, endp)
ax.set_ylim(-1, 1)
ax.axes.get_yaxis().set_visible(False)
fig.legend(loc='lower right')

# single rod
# ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
# ax2.spines['bottom'].set_position('center')
# ax2.spines['bottom'].set_visible(False)
# ax2.set_xticks(xticks[3], x_ticks_label[3] )
ax2.set_xlim(0, endp)
ax2.set_ylim(-1, 1)
ax2.axes.get_yaxis().set_visible(False)

H_list = np.array([125./44*x for x in range(1,17)])
print(H_list)
bifurcation_indices_single = [0, 10, -1]
labels_single = ['Planar beating', 'Whirling']
for i, ind0 in enumerate(bifurcation_indices_single[:-1]):
	ind1 = bifurcation_indices_single[i+1]
	ax2.plot(H_list[ind0:ind1], np.zeros(np.shape(H_list[ind0:ind1])), c='black', linestyle='', 
		  ms=10, marker=markers_single[i], label=labels_single[i])
ax2.set_xlabel(r'$H/L$')
# fig2.legend()


# 100 rods
# ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['left'].set_visible(False)
# ax3.spines['bottom'].set_position('center')
# ax3.spines['bottom'].set_visible(False)
# ax3.set_xticks(xticks[3], x_ticks_label[3] )
ax3.set_xlim(0, endp)
ax3.set_ylim(-1, 1)
ax3.axes.get_yaxis().set_visible(False)

H_list = np.array([125./44*x for x in range(1,17)])
bifurcation_indices_100 = [0, 4, 7, 10, -1]
for i, ind0 in enumerate(bifurcation_indices_100[:-1]):
	print(ind0)
	ind1 = bifurcation_indices_100[i+1]
	ax3.plot(H_list[ind0:ind1], np.zeros(np.shape(H_list[ind0:ind1])), c='black', linestyle='', 
		  ms=10, marker=markers[i], label=labels[i])
ax3.set_xlabel(r'$H/L$')
# fig3.legend()


fig.savefig('fig/fil_array_bifurcation.png', format='png', bbox_inches='tight')
fig.savefig('fig/fil_array_bifurcation.pdf', format='pdf', bbox_inches='tight')
fig2.savefig('fig/fil_array_bifurcation_single.png', format='png', bbox_inches='tight')
fig2.savefig('fig/fil_array_bifurcation_single.pdf', format='pdf', bbox_inches='tight')
fig3.savefig('fig/fil_array_bifurcation_100.png', format='png', bbox_inches='tight')
fig3.savefig('fig/fil_array_bifurcation_100.pdf', format='pdf', bbox_inches='tight')
plt.show()