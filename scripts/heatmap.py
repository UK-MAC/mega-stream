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

data = numpy.zeros((8,5))
data[0] = [71103.0, 114238.5, 94292.4, 92105.7, 52930.6]
data[1] = [147649.4, 223801.5, 251318.1, 227114.9, 196121.0]
data[2] = [252762.3, 311192.7, 294210.3, 227833.1, 185339.1]
data[3] = [310676.5, 395393.0, 302705.0, 195018.7, 0.0]
data[4] = [351479.6, 332399.7, 241249.2, 183720.3, 0.0]
data[5] = [439309.4, 294268.8, 191220.3, 168287.6, 0.0]
data[6] = [411714.6, 212903.5, 167718.5, 0.0, 0.0]
data[7] = [270262.7, 181380.7, 145228.9, 0.0, 0.0]

data *= 1.0E-3

fig, ax = plt.subplots()
plt.pcolor(data, cmap='GnBu')
ax.set_xticks(numpy.arange(data.shape[1]) + 0.5)
ax.set_yticks(numpy.arange(data.shape[0]) + 0.5)
ax.set_xticklabels([4, 8, 16, 32, 64])
ax.set_yticklabels([8, 16, 32, 64, 128, 256, 512, 1024])
ax.set_xlabel('Middle size')
ax.set_ylabel('Inner size')
plt.title('Outer size=64')
cbr = plt.colorbar()
cbr.ax.set_ylabel('Bandwidth GB/s')

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
plt.savefig('heatmap.pdf')

