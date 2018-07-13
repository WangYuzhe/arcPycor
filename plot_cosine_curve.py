# 2018-7-11
# Wang Yuzhe

import os
import numpy as np
import matplotlib.pyplot as plt

Dir = r"E:\cloud\AcademicWriting_WYZ\paper_PyCoreg\results\run_sample_1arc"

x_bin1 = np.loadtxt(os.path.join(Dir, "x_bin1.txt"))
y_bin1 = np.loadtxt(os.path.join(Dir, "y_bin1.txt"))
sigma_bin1 = np.loadtxt(os.path.join(Dir, "sigma_bin1.txt"))

fig = plt.figure(figsize=(6.5,5))
plt.errorbar(x_bin1, y_bin1, sigma_bin1)

plt.xlabel('Aspect (degrees)')
plt.ylabel(r'dh/tan($\alpha$) (m)')

ax = plt.gca()
ax.set_xlabel('Aspect (degrees)', fontsize=11)
ax.set_ylabel(r'dh/tan($\alpha$) (m)', fontsize=11)
ax.set_xlim([0, 370])
#ax.set_ylim([-5, 45])
ax.set_xticks(np.arange(0, 370, 40))
ax.set_xticklabels(np.arange(0, 370, 40))
ax.tick_params(axis='both', labelsize=11)
plt.show()

plt.tight_layout()