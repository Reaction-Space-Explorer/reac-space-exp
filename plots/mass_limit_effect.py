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

axes[0].plot(gen, products_200, color='gray', marker='^', ms=9, mec='k', mew=0.2,
			 mfc='cornflowerblue', ls='--', label='200')
axes[0].plot(gen, products_250, color='gray', marker='o', ms=8, mec='k', mew=0.2,
			 mfc='deeppink', ls='-', label='250')
axes[0].plot(gen, products_300, color='gray', marker='*', ms=11, mec='k', mew=0.2,
			 mfc='cornflowerblue', ls='-.', label='300')
axes[0].set_xlabel('Generation', fontsize=13)
axes[0].set_ylabel('New Products', fontsize=13)
axes[0].legend(loc='upper left', title='Mass limit', fontsize=12, title_fontsize=12)
axes[0].set_xticks([1, 2, 3])
axes[0].set_yscale('log')
axes[1].plot(gen, time_200, color='gray', marker='^', ms=9, mec='k', mew=0.2,
			 mfc='cornflowerblue', ls='--', label='200')
axes[1].plot(gen, time_250, color='gray', marker='o', ms=8, mec='k', mew=0.2,
			 mfc='deeppink', ls='-', label='250')
axes[1].plot(gen, time_300, color='gray', marker='*', ms=11, mec='k', mew=0.2,
			 mfc='cornflowerblue', ls='-.', label='300')
axes[1].set_xlabel('Generation', fontsize=13)
axes[1].set_ylabel('Time (seconds)', fontsize=13)
axes[1].legend(loc='upper left', title='Mass limit', fontsize=12, title_fontsize=12)
axes[1].set_xticks([1,2,3])
axes[1].set_yscale('log')



# make ticks larger:
axes[0].tick_params(axis='both', labelsize=13)
axes[1].tick_params(axis='both', labelsize=13)

plt.tight_layout()
plt.savefig('mass_limit_effect.jpg', dpi=300)
plt.show()
