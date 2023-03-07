import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
# Pie chart, where the slices will be ordered and plotted c
# ounter-clockwise:
labels =  ['Updating solution', 'Applying force',  'Evaluating error', 'Updating Jacobian', 'Fast FCM']
sizes = [0.0233855, 0.000254995, 0.00735886, 0.0251631, 0.0240049]
explode = (0, 0, 0, 0, 0.1)  # only "explode" the 2nd slice (i.e. 'Hogs')

# ax.set_title("Time of Fast FCM")
ax.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90)
ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.savefig('fig/pie_chart_of_rod.eps', format='eps')
plt.show()