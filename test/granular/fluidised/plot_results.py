import numpy as np 
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(12.8, 9.6))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
plt.ion()
plt.show()

dragModels = ['DiFelice']
Uf_values = [1.85, 1.86, 1.87, 1.88, 1.89, 1.90]
for Uf in Uf_values:
    data = np.load('./results/fig_Uf_{}_DiFelice.npz'.format(Uf))
    t = data['t']
    xy = data['xy']
    xySol = data['xySol']

    ax1.plot(t, xy, linewidth=3.0, label='Uf = {:.2f}'.format(Uf))
    ax1.plot(t, np.ones_like(t)*xySol, 'k--')
    plt.pause(0.1)

dragModels = ['Ergun']
Uf_values = [1.60, 1.61, 1.62, 1.63, 1.64, 1.65]
for Uf in Uf_values:
    data = np.load('./results/fig_Uf_{}_Ergun.npz'.format(Uf))
    t = data['t']
    xy = data['xy']
    xySol = data['xySol']

    ax2.plot(t, xy, linewidth=3.0, label='Uf = {:.2f}'.format(Uf))
    ax2.plot(t, np.ones_like(t)*xySol, 'k--')
    plt.pause(0.1)

ax1.set_title('Di Felice')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position (cm)')
ax1.set_ylim((9.948, 9.952))
ax1.legend()
ax2.set_title('Ergun')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Position (cm)')
ax2.set_ylim((9.948, 9.952))
ax2.legend()
plt.tight_layout()
plt.savefig('fig_displacement.png')
plt.close()

