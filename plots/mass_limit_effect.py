import matplotlib.pyplot as plt

# sorry for the spaghetti code, I know it can be written much more compactly but I was rushing.

gen = [1, 2, 3]

time_200 = [0.06462, 0.79871, 26.78112]
products_200 = [16, 110, 926]

time_250 = [0.06683, 1.07613, 117.04779]
products_250 = [16, 148, 4537]

time_300 = [0.06804, 1.02086, 571.16982]
products_300 = [16, 184, 13832]

fig, axes = plt.subplots(1, 2, sharex=True, figsize=(12, 6))

axes[0].plot(gen, products_200, color='magenta', marker='o', linestyle='-.', label='200')
axes[0].plot(gen, products_250, color='gray', marker='o', linestyle='-.', label='250')
axes[0].plot(gen, products_300, color='royalblue', marker='o', linestyle='-.', label='300')
axes[0].set_xlabel('Generation')
axes[0].set_ylabel('New Products')
axes[0].legend(loc='upper left')
axes[0].set_xticks([1, 2, 3])
axes[0].set_yscale('log')

axes[1].plot(gen, time_200, color='magenta', marker='o', linestyle='--', label='200')
axes[1].plot(gen, time_250, color='gray', marker='o', linestyle='--', label='250')
axes[1].plot(gen, time_300, color='royalblue', marker='o', linestyle='--', label='300')
axes[1].set_xlabel('Generation')
axes[1].set_ylabel('Time (seconds)')
axes[1].legend(loc='upper left')
axes[1].set_xticks([1,2,3])
axes[1].set_yscale('log')

plt.show()