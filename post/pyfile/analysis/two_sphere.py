import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Nblob = 1400
# Fourier modes = 0.5 0
a = 20.533623615045787
v_isolate = 2.8145485528e-02
r_list = np.array([ 
                    42.58, 43, 43.4, 43.8,
                    44.2, 44.6, 45.0, 45.4, 45.8,
                    46.2, 46.6, 47.0, 47.4, 47.8,
                    49.0, 50.0, 55.0,\
                    60, 65, 70, 75, 80,\
                    85, 90, 95, 100, 105,\
                    120, 135, 150, 165, 180,\
                    200, 220, 240, 260, 280
                    ])
v_prl_list = np.array([
    4.3406154961e-02, 4.3330605504e-02, 4.3257478007e-02, 4.3184814740e-02,
    4.3112613376e-02, 4.3040872697e-02, 4.2969592026e-02, 4.2898770844e-02, 4.2828408515e-02,
    4.2758504110e-02, 4.2689056310e-02, 4.2620063376e-02, 4.2551523160e-02, 4.2483433147e-02,
    4.2281834790e-02, 4.2116847789e-02, 4.1330910102e-02,
    4.0604915553e-02, 3.9933188305e-02, 3.9311993202e-02, 3.8738259007e-02, 3.8209018617e-02,
    3.7706797263e-02, 3.7258035973e-02, 3.6844123261e-02, 3.6462003039e-02, 3.6108788790e-02,
    3.5196939310e-02, 3.4462238311e-02, 3.3860622583e-02, 3.3360424574e-02, 3.2938786281e-02,
    3.2470567208e-02, 3.2084481573e-02, 3.1760965101e-02, 3.1486128594e-02, 3.1249868122e-02])

v_ognl_list = np.array([
    3.9425866743e-02, 3.9304687476e-02, 3.9213277931e-02, 3.9118759013e-02,
    3.9023119346e-02, 3.8927971354e-02, 3.8832056698e-02, 3.8736686569e-02, 3.8637369907e-02,
    3.8543939067e-02, 3.8451308994e-02, 3.8359559074e-02, 3.8268756036e-02, 3.8178954995e-02,
    3.7915964803e-02, 3.7700578582e-02, 3.6749970332e-02,
    3.5958657052e-02, 3.5296898002e-02, 3.4737697762e-02, 3.4259881445e-02, 3.3847293250e-02,
    3.3487614337e-02, 3.3171372551e-02, 3.2891196575e-02, 3.2641270680e-02, 3.2416940368e-02,
    3.1863043569e-02, 3.1437938325e-02, 3.1100687974e-02, 3.0826456227e-02, 3.0599053635e-02,
    3.0350038381e-02, 3.0147127153e-02, 2.9978565764e-02, 2.9836288538e-02, 2.9714576130e-02])

v_prl_opposite_list = np.array([
    3.2234484902e-03, 3.6526457050e-03, 4.0739463096e-03, 4.4953573372e-03,
    4.9080292666e-03, 5.3167951974e-03, 5.6956523424e-03, 6.0554732715e-03, 6.3977933566e-03,
    6.7179726383e-03, 7.0319922228e-03, 7.3339519786e-03, 7.6233636424e-03, 7.9007005883e-03,
    8.6681911738e-03, 9.2421043933e-03, 1.1534329988e-02, 
    1.3202451195e-02, 1.4510513982e-02, 1.5574863446e-02, 1.6467676367e-02, 1.7231014115e-02,
    1.7893604482e-02, 1.8465839287e-02, 1.8982874927e-02, 1.9445157886e-02, 1.9861431328e-02,
    2.0896298251e-02, 2.1697776546e-02, 2.2338115520e-02, 2.2862018907e-02, 2.3298850198e-02,
    2.3779822733e-02, 2.4173827395e-02, 2.4502569812e-02, 2.4781069693e-02, 2.5020058832e-02])

v_ognl_opposite_list = np.array([
    1.4986735686e-02, 1.5292269218e-02, 1.5557519968e-02, 1.5795395553e-02,
    1.6012723441e-02, 1.6213859536e-02, 1.6401592798e-02, 1.6579429716e-02, 1.6749304213e-02,
    1.6907561611e-02, 1.7058641512e-02, 1.7203339474e-02, 1.7342302445e-02, 1.7476064369e-02,
    1.7850292778e-02, 1.8135831720e-02, 1.9300624423e-02,
    2.0196683655e-02, 2.0912402324e-02, 2.1488273591e-02, 2.1983688388e-02, 2.2407240979e-02,
    2.2774161088e-02, 2.3095499905e-02, 2.3379518310e-02, 2.3632538049e-02, 2.3859491358e-02,
    2.4419434665e-02, 2.4848306713e-02, 2.5187629461e-02, 2.5462971644e-02, 2.5690985128e-02,
    2.5940429351e-02, 2.6143556260e-02, 2.6312232538e-02, 2.6454574196e-02, 2.6576324411e-02])

########################################################################
# Batchelor (1982)
########################################################################
def A11(p, lam):
    return 1 - 60*lam**2/(1+lam)**4/p**4 - 192*lam**3*(5-22*lam**2+3*lam**4)/(1+lam)**8/p**8

def A12(p, lam):
    return 1.5/p - 2*(1+lam**2)/((1+lam)**2*p**3) + 1200*lam**3/(1+lam)**6/p**7

def B11(p, lam):
    return 1 - 68*lam**5/(1+lam)**6/p**6 - 32*lam**3*(10-9*lam**2+9*lam**4)/(1+lam)**8/p**8

def B12(p, lam):
    return 0.75/p + (1+lam**2)/(1+lam)**2/p**3

########################################################################
# Plot
########################################################################
vw_prl_list = v_prl_list/v_isolate
vw_ognl_list = v_ognl_list/v_isolate
vw_prl_opposite_list = v_prl_opposite_list/v_isolate
vw_ognl_opposite_list = v_ognl_opposite_list/v_isolate

ra_batchelor_list = np.linspace(2, 14, 100)
vw_batchelor_prl_list = A11(ra_batchelor_list, 1) + A12(ra_batchelor_list, 1)
vw_batchelor_ognl_list = B11(ra_batchelor_list, 1) + B12(ra_batchelor_list, 1)
vw_batchelor_prl_opposite_list = A11(ra_batchelor_list, 1) - A12(ra_batchelor_list, 1)
vw_batchelor_ognl_opposite_list = B11(ra_batchelor_list, 1) - B12(ra_batchelor_list, 1)

ax = plt.figure().add_subplot(1,1,1)
ax.scatter(r_list/a, vw_prl_list, facecolors='none', edgecolors='r', label='Parallel alignment')
ax.scatter(r_list/a, vw_ognl_list, facecolors='none', edgecolors='b', label='Othogonal alignment')
ax.scatter(r_list/a, vw_prl_opposite_list, facecolors='none', edgecolors='r', label='Parallel (opposite) alignment')
ax.scatter(r_list/a, vw_ognl_opposite_list, facecolors='none', edgecolors='b', label='Othogonal (opposite) alignment')
ax.plot(ra_batchelor_list, vw_batchelor_prl_list, c='r', label='Parallel alignment Batchelor(1982)')
ax.plot(ra_batchelor_list, vw_batchelor_ognl_list, c='b', label='Othogonal alignment Batchelor(1982)')
ax.plot(ra_batchelor_list, vw_batchelor_prl_opposite_list, c='r', label='Parallel alignment (opposite) Batchelor(1982)')
ax.plot(ra_batchelor_list, vw_batchelor_ognl_opposite_list, c='b', label='Othogonal alignment (opposite) Batchelor(1982)')
ax.axhline(y=1.0, color='grey', linestyle='-.')
ax.annotate('Force in same direction', (2.1, 1.02))
ax.annotate('Force in opposite direction', (2.1, 0.93))
ax.set_xlabel(r'$r/a$')
ax.set_ylabel(r'V/W')
ax.set_title('Two sphere settling speed')
ax.set_xlim(2, 3.2)
ax.set_ylim(0, 1.7)
red_patch = mpatches.Patch(color='r', label='Parallel alignment')
blue_patch = mpatches.Patch(color='b', label='Othogonal alignment')
circle_legend = ax.scatter([], [], facecolors='none', edgecolors='black', label='Rigid sphere data')
line_legend = mlines.Line2D([], [], c='black', label='Batchelor(1982)')
ax.legend(handles=[red_patch, blue_patch, circle_legend, line_legend])

plt.savefig('two_sphere.eps', format='eps')
plt.show()







