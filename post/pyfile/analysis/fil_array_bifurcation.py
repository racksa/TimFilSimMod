import numpy as np
import matplotlib.pyplot as plt

endp = 55
fig = plt.figure(figsize=(12,3.5))
ax = fig.add_subplot(1,1,1)

xticks  = [0, 12, 20, 29.54, endp]
x_ticks_label = ["H={:.2f}L".format(x) for x in xticks]
style = ['-', ':', '--', '-.']
width = [2, 2, 2, 3]
labels = ['Chaotic motion', 'Planar beating', 'Synchronised planar beating', 'Whirling']

line_pos = [2.84, 14.2, 22.73, 31.25]

for i, pos in enumerate(line_pos):
	b=0
	if(i%2==0):
		print(i)
		b=-1
	ax.vlines(pos, b, b+1, color='black', ls=style[i], linewidth=1)
	
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
plt.legend(loc='lower right')
plt.savefig('fig/fil_array_bifurcation.png', format='png', bbox_inches='tight')
plt.savefig('fig/fil_array_bifurcation.pdf', format='pdf', bbox_inches='tight')
plt.show()