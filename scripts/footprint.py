#
#
#  Copyright 2016 Tom Deakin, University of Bristol
#
#  This file is part of mega-stream.
#
#  mega-stream is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  mega-stream is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with mega-stream.  If not, see <http://www.gnu.org/licenses/>.
#
#
#  This aims to investigate the limiting factor for a simple kernel, in particular
#  where bandwidth limits not to be reached, and latency becomes a dominating factor.
#
#

import numpy
import matplotlib.pyplot as plt

Nm = 64
Ni = [8, 16, 32, 64, 128, 256, 512, 1024]
Nj = [4, 8, 16, 32, 64]

data = numpy.zeros((8,5))

for i in range(len(Ni)):
  for j in range(len(Nj)):
    data[i][j] = 8.0 * 1.0E-6 * (
      2.0 * Ni[i]*(Nj[j]**3)*Nm +
      3.0 * Ni[i]*(Nj[j]**2)*Nm +
      3.0 * Ni[i] +
      (Nj[j]**3)*Nm
    )

# zero too big results
data[3][4] = 0.0
data[4][4] = 0.0
data[5][4] = 0.0
data[6][3:] = 0.0
data[7][3:] = 0.0
        
fig, ax = plt.subplots()
plt.pcolor(data, cmap='GnBu')
ax.set_xticks(numpy.arange(data.shape[1]) + 0.5)
ax.set_yticks(numpy.arange(data.shape[0]) + 0.5)
ax.set_xticklabels(Nj)
ax.set_yticklabels(Ni)
ax.set_xlabel('Middle size')
ax.set_ylabel('Inner size')
plt.title('Outer size=64')
cbr = plt.colorbar()
cbr.ax.set_ylabel('Footprint MB')

# Add data labels
for i in range(data.shape[1]):
  for j in range(data.shape[0]):
    if (data[j][i] != 0.0):
      plt.text(i + 0.5, j + 0.5, '%.1f' % (data[j][i]),
        ha='center', va='center',
        size='small', color='black', weight='bold')
    else:
      plt.text(i + 0.5, j + 0.5, '-',
        ha='center', va='center',
        size='small', color='black', weight='bold')

#fig.set_tight_layout(True)
plt.savefig('footprint.pdf')

