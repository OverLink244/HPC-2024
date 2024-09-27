import matplotlib.pyplot as plt

x = ['Origin', 'V1.0', 'V1.1.1', 'V1.1.2', 'V1.1.3', 'V1.2.1', 'V1.2.2','V1.2.3']
y1 = [7.47, 10.74, 11.27, 12.64, 3.27, 12.48,17.16, 24.13]
y2=[16.59,16.55,18.10,17.29,14.10,20.46,21.83,214.42]

plt.plot(x, y1, label='Small.conf', linestyle='-', color='green')
plt.plot(x, y2, label='Realworld.conf', linestyle='--', color='red')
plt.xlabel('version')
plt.ylabel('GFlops')
# plt.ylim([0, 200])
plt.legend(loc='best')
plt.savefig('PAC_result.png', dpi=500)
plt.show()
