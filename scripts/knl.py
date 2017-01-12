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
from matplotlib.ticker import ScalarFormatter

inner = numpy.array([96.4, 248.3, 296.9, 302.7, 237.6, 190.5, 165.5])
inner_labels = [8, 16, 32, 64, 128, 256, 512]
outer = numpy.array([240.2, 236.8, 239.2, 240.1, 239.2])
outer_labels = [64, 128, 256, 512, 1024]
middle = numpy.array([351.5, 345.1, 236.7, 184.2])
middle_labels = [4, 8, 16, 32]

fig, ax = plt.subplots()
ax.set_xscale('log', basex=2)
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.ylabel('Bandwidth GB/s')
plt.xlabel('Size')
plt.title('Varying one parameter\nDefaults: Inner=128, Middle=16, Outer=64')

plt.plot(inner_labels, inner, '-x', label='Inner')
plt.plot(outer_labels, outer, '-+', label='Outer')
plt.plot(middle_labels, middle, '-o', label='Middle')

plt.legend()

plt.savefig('knl.pdf')

