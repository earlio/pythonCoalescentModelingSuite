## EZ 2018 Changes to this file were required to make it Python 3 compatible.
## search for 'EZ' in this file to find those changes (other than print statement updates). earl.io
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

python script to simulate runaway sexual selection

date: 10-10-2011

usage:
    python sexsel.py

check help message:
    python sexsel.py -h


'''

import pylab,random,sys,optparse

def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack 
    result = [None] * k
    for i in xrange(k): # xrange: faster, less memory
        j = _int(_random() * n)
        result[i] = population[j]
    return result


def compute_mate_selection_probability(genotype_numbers,mate_selection_preference):
    '''
    compute mate selection probability based on mate preference and genotype numbers
    '''
    #print 'genotype_numbers',genotype_numbers
    #print 'mate_selection_preference',mate_selection_preference
    mate_selection_probability = pylab.zeros(mate_selection_preference.shape)
    for i in range(4):
        mate_selection_probability[i,] = mate_selection_preference[i,]*genotype_numbers
        mate_selection_probability[i,] = mate_selection_probability[i,]/sum(mate_selection_probability[i,])
    #print 'mate_selection_probability',mate_selection_probability
    return mate_selection_probability

def initial_population(genotype_numbers):
    '''
    generate individuals based on frequency specified genotype numbers
    '''

    # total number of individuals to generate
    ## EZ 
    # apparently we need to explicitly cast sum as int now (Python 2 to 3)
    n_individual = int(sum(genotype_numbers))

    # compute tally
    tally = [0]*5
    for i in range(4):
        tally[i+1] =tally[i] + genotype_numbers[i]

    # initialize population
    population = pylab.zeros(n_individual,dtype=int)
    for i in range(4):
        print("EZ tally'i' = ", tally)
        print("EZ population size = ", len(population))
        print("EZ pop[0] = ", population[0])
        ## EZ here we assign population members o to x the value 0, x to x2 the value 1, etc.
        # apparently this assignment no longer works in Python 3:
        #population[tally[i]:tally[i+1]] = i  
        for j in range(int(tally[i]), int(tally[i+1])):
            population[j] = i
 
    return population

def mate_selection(population,mate_selection_probability):
    '''
    select mate based on selection probability
    '''

    # number of individuals need to select mate
    n_individual = len(population)

    # compute tally of mate_selection_probability
    mate_selection_probability_tally = pylab.zeros((4,5))
    for i in range(4):
        mate_selection_probability_tally[:,i+1] = mate_selection_probability_tally[:,i] + mate_selection_probability[:,i] 

    # to store the selected mates for each individual
    selected_mate = pylab.zeros(n_individual,dtype=int)

    # generate a list of random number
    r = pylab.rand(n_individual)

    # select mate for each individual in the population
    for i in range(n_individual):
        for j in range(4):
            if mate_selection_probability_tally[population[i],j] < r[i] and mate_selection_probability_tally[population[i],j+1] >= r[i]:
                selected_mate[i] = j
                break

    return selected_mate

def natural_selection(selected_mate,survival_rate):
    '''
    randomly kill individuals with low survival rate
    killed individuals will have genotype -1
    '''
    r = pylab.rand(len(selected_mate))
    selected_mate[survival_rate[selected_mate] < r] = -1
    #print 'killed',sum(selected_mate == -1),'individuals'
    return selected_mate  

def breeding(female_population,survived_mate,fix_population_size):
    '''
    generate offspring
    '''

    male_offspring = []
    female_offspring = []

    idx = pylab.array(range(len(survived_mate)))
    survived = idx[survived_mate > -1]
    
    if fix_population_size:
        ## sample from survived males
        survived = sample_wr(survived, len(survived_mate))

    for i in survived:
        #make diploid genotype
        diploid_string = genotype[female_population[i]]+genotype[survived_mate[i]]
        #make a male haploid genotype
        male_haploid = diploid2haploid(diploid_string)
        #make a female haploid genotype
        female_haploid = diploid2haploid(diploid_string)
        #add to output
        male_offspring.append(genotype2number[male_haploid])
        female_offspring.append(genotype2number[female_haploid])

    return pylab.array(male_offspring),pylab.array(female_offspring)

def diploid2haploid(diploid):
    '''
    generate haploid from diploid in string form
    '''
    t = random.sample([0,4],1)[0] # T1 or T2
    p = random.sample([2,6],1)[0] # P1 or P2
    #print diploid,t# diploid[t:t+1]
    haploid = diploid[t:t+2]+diploid[p:p+2]
    return haploid

def genotype_count_freq(population):
    '''
    count genotype number and frequency
    '''
    genotype_count = pylab.zeros(4)
    for i in range(4):
        genotype_count[i] = sum(population == i)
    genotype_freq = genotype_count/float(sum(genotype_count))
    
    return genotype_count,genotype_freq

def simulation(N,n_generation,s,v,gnm,gnf,fix):    
    # starting genotype frequency
    male_genotype_freq = pylab.fromstring(gnm,sep=',')
    male_genotype_freq = male_genotype_freq/sum(male_genotype_freq)
    female_genotype_freq = pylab.fromstring(gnf,sep=',')
    female_genotype_freq = female_genotype_freq/sum(female_genotype_freq)
    
    # mate selection preference
    mate_selection_preference = pylab.fromstring(s.replace('\n','\t'),dtype='float',sep='\t') 
    mate_selection_preference.shape = [4,4]

    print('mate_selection_preference\n',mate_selection_preference)
    
    # male survival rate
    male_survival_rate = pylab.array([1,1,v,v])
    print('male_survival_rate\n',male_survival_rate)

    '''
    initialization
    '''

    # starting genotype numbers
    male_genotype_number = N*male_genotype_freq
    female_genotype_number = N*female_genotype_freq

    # compute mate selection probability
    mate_selection_probability = compute_mate_selection_probability(male_genotype_number,mate_selection_preference)

    # initialize male and female populations based on genotype numbers
    male_population = initial_population(male_genotype_number)
    female_population = initial_population(male_genotype_number)
    
    # compute initial frequency
    male_genotype_freqs = pylab.zeros((n_generation+1,4))
    male_genotype_freqs[0,:] = male_genotype_freq
    female_genotype_freqs = pylab.zeros((n_generation+1,4))
    female_genotype_freqs[0,:] = female_genotype_freq


    '''
    simulation
    '''

    population_size = [N]
    for i in range(n_generation):

        # sex selection
        selected_mate = mate_selection(female_population,mate_selection_probability)
        #print 'selected mate:',genotype_count_freq(selected_mate )

        # natural selection
        survived_mate = natural_selection(selected_mate,male_survival_rate)

        # breeding
        male_population,female_population = breeding(female_population,survived_mate,fix)

        # update genotype number and frequency
        male_genotype_number, male_genotype_freqs[i+1,:] = genotype_count_freq(male_population)
        female_genotype_number, female_genotype_freqs[i+1,:] = genotype_count_freq(female_population)

        # update population size
        population_size.append(sum(male_genotype_number))        

        # update mate selection probability
        mate_selection_probability = compute_mate_selection_probability(male_genotype_number,mate_selection_preference)


    return male_genotype_freqs,female_genotype_freqs,population_size

def draw(fig,male_genotype_freqs,female_genotype_freqs):    
    '''
    visualization
    '''
    axes = fig.add_subplot(211)
    # plot male frequency
    axes.plot(male_genotype_freqs)
    axes.set_ylim(0,1)
    axes.legend(genotype,bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=len(genotype), mode="expand", borderaxespad=0.)
    #axes.set_title('Change in Genotype Frequencies: Males')
    axes.set_xlabel('Male Generation')
    axes.set_ylabel('Genotype frequencies')

    # plot female frequency
    axes2 = fig.add_subplot(212)
    axes2.plot(female_genotype_freqs)
    axes2.set_ylim(0,1)
    axes2.legend(genotype,bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=len(genotype), mode="expand", borderaxespad=0.)
    #axes2.set_title('Change in Genotype Frequencies: Females')
    axes2.set_xlabel('Female Generation')
    axes2.set_ylabel('Genotype frequencies')
    fig.subplots_adjust(hspace=.5)

# for the simplicity of programming, genotypes are coded as numbers: 0,1,2,3
genotype = ['T1P1', 'T1P2', 'T2P1', 'T2P2']
genotype2number = {'T1P1':0, 'T1P2':1, 'T2P1':2, 'T2P2':3}
    
def main():

    # Parse command line    
    parser = optparse.OptionParser( usage="%prog [options] " )
    parser.add_option("-N","--population-size", type='int',dest="N",default=1000,
                      help="Population size (N male and N female)(default N=1000) " )
    parser.add_option("-k","--num-generation", dest="k",type='int',default=20,
                      help="Number of generations to simulate (default=20)" ) 
    parser.add_option("-v","--survival probability",dest="v",type='float',default='0.0',
                      help="Survival rate of male with long tail (T2) (0<=v<=1,default v=0.0, no natural selection against long tail (T2))" )       
    parser.add_option("-s","--mate-selection-preference",dest="s",type='string',default='0.25\t0.25\t0.25\t0.25\n0\t0\t0.5\t0.5\n0.25\t0.25\t0.25\t0.25\n0\t0\t0.5\t0.5',
                      help="Mate selection preference matrix (T1P1,T1P2,T2P1,T2P2)" ) 
    parser.add_option("-f","--fix-population-size",action='store_true',dest="f",
                      help="Fix population size" )  
    parser.add_option("--genotype-number-male",dest="gnm",default='250,250,250,250',
                      help="genotype number in male: T1P1,T1P2,T2P1,T2P2 (default:250,250,250,250)" )                                             
    parser.add_option("--genotype-number-female",dest="gnf",default='250,250,250,250',
                      help="genotype number in female: T1P1,T1P2,T2P1,T2P2 (default:250,250,250,250)" )   
                      
    options, args = parser.parse_args()

    '''
    parameters 
    '''

    # number of generations to simulate
    n_generation = options.k

    # population size (N male and N female)
    N = options.N

    male_genotype_freqs,female_genotype_freqs, population_size = simulation(N,n_generation,options.s,options.v,options.gnm,options.gnf,options.f)
    
    draw(pylab.figure(),male_genotype_freqs,female_genotype_freqs)
    pylab.show()



if __name__ == "__main__":
    main()
