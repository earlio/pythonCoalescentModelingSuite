## EZ 2018 this file required no modifications to work with Python 3. earl.io
'''
Copyright (C) 2011, 2014 by Xuebing Wu and Robert C. Berwick

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''

from array import *
from random import *
from math import *
import pylab

def poisson(mean):
    sum = -log(random())
    i = 0
    while sum < mean:
        sum += -log(random())
        i = i + 1
    return i
def muller(pop, U, fitness,mutmax):
    minimum = [0]; wMax = 1.0
    popn = []
    for i in range(pop): popn.append(0)
    generation = [0]
    cur_generation = 0
    while minimum[-1] < mutmax:
        for i in range(pop): popn[i] = popn[i] + poisson(U)
        i = 0
        tmpPopn = []
        while i < pop:
            n = choice(popn)
            if random() < fitness**n / wMax:
                tmpPopn.append(n)
                i += 1
        popn = tmpPopn
        cur_generation += 1
        if min(popn) > minimum[-1]:
            minimum.append(min(popn))
            generation.append(cur_generation)
            wMax = fitness**minimum[-1]
    #pylab.plot(generation,minimum)
    #pylab.show()
    return generation,minimum

