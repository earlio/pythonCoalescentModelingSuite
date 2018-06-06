## EZ 2018 This file required only print statement updates to work with Python 3.
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
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
python script to simulate drift and selection

author: Xuebing Wu (wuxbl@mit.edu)
date: 10-10-2011

Usage:
    python diploid_simulation.py

    Check the help message to change parameters:

    python diploid_simulation.py -h

'''


import pylab,optparse,scipy.stats

def drift(freq,N):
    '''
    Draw N samples from a multinomial distribution with parameter freq, 
    return the genotype numbers
    
    freq = frequency of [AA,Aa,aa]
    '''

    # first, draw AA genotypes
    nAA = scipy.stats.binom.rvs(N,freq[0])
    
    # second, draw Aa genotypes
    nAa = 0
    if freq[2] == 0:    # if aa frequency is 0, can't do binomial(n,1)
        nAa = N-nAA
    elif nAA < N:   # if AA doesn't take the whole population, can't do binomial(0,p)
        nAa = scipy.stats.binom.rvs(N-nAA,freq[1]/(1-freq[0]))
    naa = N-nAA-nAa

    return nAA,nAa,naa

def mutation(gN,u,N):
    '''
    simulate mutations
    N = population size
    gN = genotype numbers of [AA,Aa,aa]
    u = mutation rate of [A2a,a2A]
    '''

    # number of A allele
    nA = gN[0]*2+gN[1]  
      
    # draw A2a mutations if nA > 0. can't do binomial(0,p)
    A2a = 0
    if nA > 0:
        A2a = scipy.stats.binom.rvs(nA,u[0])
    
    # draw a2A mutations if na > 0. can't do binomial(0,p)
    a2A = 0
    if nA < 2*N:
        a2A = scipy.stats.binom.rvs(2*N-nA,u[1])

    # total number of A 
    nA = nA - A2a + a2A

    
    # return the frequency of A
    return float(nA)/N/2    
    
def migration(p,f_migrant,p_migrant):
    '''
    f_migrant: fraction of the new population that will be migrants
    p_migrant: allele A frequency of the migrants     
    '''
    return p * (1 - f_migrant) + f_migrant * p_migrant

def selection(p,fitness,freq_dependent):
    '''
    compute selection-adjusted frequency based on allele A frequency and and fitness
          
    '''
    
    # compute Hardy-Weinberg genotype frequency before selection
    HDfreq = pylab.array([p*p,2*p*(1-p),(1-p)*(1-p)])

    # frequency dependent selection
    if freq_dependent >= 0:
        # equation 7.19, on page 218 of Hamilton's book
        # fitness = 1-HDfreq*(1-fitness)
        # Rice's book, p15
        fitness[0] = 1 - 3*HDfreq[1]+3*HDfreq[2]
        fitness[1] = 1 - freq_dependent*HDfreq[1]
        fitness[2] = 1 - 3*HDfreq[1]+3*HDfreq[0] 
        
    # compute population mean fitness
    population_mean_fitness = sum(HDfreq*fitness)

    # adjust frequency
    freq_new = HDfreq*fitness/population_mean_fitness

    return freq_new    

def diploid_simulation(N,p,u,fitness,f_migrant,p_migrant,t,freq_dependent):
    '''
    simulation for t generation
    starting population size N
    starting allele A frequency p
    relative fitness of [AA, Aa, aa]: fitness
    '''

    # record actual frequency of [AA,Aa,aa] for each generation
    freq_matrix = pylab.zeros((t+1,3))    
    
    # the initial frequency
    freq_matrix[0,] = pylab.array([p*p,2*p*(1-p),(1-p)*(1-p)])
    
    # for each generation
    for i in range(t):
        # migration
        p = migration(p,f_migrant,p_migrant)
        # selection
        freq_selec = selection(p,fitness,freq_dependent)
        # drift
        nAA,nAa,naa = drift(freq_selec,N)
        # mutation
        p = mutation([nAA,nAa,naa],u,N)
        
        freq_matrix[i+1,] = pylab.array([nAA,nAa,naa])/float(N)

    return freq_matrix

def plotfreq(axes,freq_matrix,title):
    for i in range(freq_matrix.shape[0]):
        axes.plot(freq_matrix[i,:])
    axes.set_ylim(0,1)
    axes.set_title(title)
    axes.set_xlabel("generation")
    axes.set_ylabel("frequency")

def simulation(N,k,t,p,u0,fitness,f_migrant,p_migrant,freq_dependent):                    
    # compute relative fitness
    fitness = pylab.fromstring(fitness,sep=',')
    print('relative fitness [AA,Aa,aa]=\t',fitness)

    # get mutation rate
    u = pylab.fromstring(u0,sep=',')
    if len(u) == 1:
        u = pylab.array([u[0],u[0]])
    print('mutation rate [A2a,a2A] =',u)

    # run simulation
    freq_matrix = pylab.zeros((k,t+1,3))

    for i in range(k):
        freq_matrix[i,] = diploid_simulation(N,p,u,fitness,f_migrant,p_migrant,t,freq_dependent)

    return freq_matrix

def draw(fig,freq_matrix):
    
    # compute allele A frequency
    freq_A = freq_matrix[:,:,0] + freq_matrix[:,:,1]/2    
    #print freq_A

    # compute allele A fixation probability
    #print freq_A[:,options.t-1]
    pAfix = sum(freq_A[:,freq_matrix.shape[1]-1] == 1)/float(freq_matrix.shape[0])
    print('allele A fixation probability = ',pAfix)
    # 
    
    
    #pylab.hist(freq_A[:,freq_matrix.shape[1]-1],bins=20)
    #pylab.title('final allele A frequency')
    
    # plot allele A frequency
    plotfreq(fig.add_subplot(211),freq_A,"frequency of allele A1 at each generation, p(fix)="+str(pAfix))

    # plot histogram of final allele A frequency
    axes = fig.add_subplot(212)
    axes.hist(freq_A[:,freq_matrix.shape[1]-1],bins=20,range=(0,1))
    axes.set_title("Distribution of final frequency of allele A1")
    axes.set_xlabel("frequency")
    axes.set_ylabel("number of populations")

    fig.subplots_adjust(hspace=.5,wspace=.5)
    pylab.show()
            
def main():

    # Parse command line    
    parser = optparse.OptionParser( usage="%prog [options] " )
    parser.add_option("-N","--population-size", type='int',dest="N",default=10,
                      help="Number of individuals in the population (default=10) " )
    parser.add_option("-p","--allele-A-freq", dest="p",type='float',default=0.5,
                      help="Starting allele A frequency (default=0.5)" ) 
    parser.add_option("-f","--relative-fitness",dest="f",default='1.,1.,1.',
                      help="Relative fitness (default=1.0,1.0,1.0). The relative fitness of AA, Aa, and aa" )     
    parser.add_option("-u","--mutation-rate",dest="u",default='0.0,0.0',
                      help="Rate of A-to-a and a-to-A mutation, seperate by ',' (default=0.0,0.0) ") 
    parser.add_option("-t","--num-generation", type='int',dest="t",default=100,
                      help="Number of generations to simulate (default=100) ")    
    parser.add_option("-k","--num-simulation", type='int',dest="k",default=20,
                      help="Number of simulations to run (default=20) ")
    parser.add_option("-m","--f-migrant", type='float',dest="m",default=0.0,
                      help="fraction of migrants each generation (default=0.0) ")    
    parser.add_option("-q","--p-migrant", type='float',dest="q",default=1.0,
                      help="frequency of A in source population (default=1.0) ")
    parser.add_option("-d","--freq-dependent",type='float',dest="d",default=-1.0,
                      help="frequency-dependent selection (s)" )  
    options, args = parser.parse_args()

    # get starting frequency of A (p)
    p =  0.5
    if 0 <= options.p <= 1:
        p = options.p
        print('starting frequency of A: p =',p)
    else:
        print("p should be in the range [0,1]!")
        exit()

    f_migrant = options.m
    p_migrant = options.q
    print('fraction migrants =', f_migrant)
    print('allele A frequency of migrants =', p_migrant)
        
    
    k = options.k
    t = options.t
    N = options.N
    
    print('population size =',N)
    print('simulate',t,'generations')
    print('simulate',k,'populations')   

    freq_matrix = simulation(N,k,t,p,options.u,options.f,f_migrant,p_migrant,options.d)
    draw(pylab.figure(),freq_matrix)
    pylab.show()


if __name__ == "__main__":
    main()
