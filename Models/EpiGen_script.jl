
#Credits: this code, and its embedded documentation, was ported to Julia from the original Matlab implementation of [EpiGen](https://github.com/LevineLab/EpiGen).
#Please refer to [their paper](https://www.biorxiv.org/content/10.1101/637272v1) if you use this model. And if in doubt then please use the original code.

using StatsBase, Distributions, DataFrames, Dates, Random

using CairoMakie

"""
    one_plot(results)

Plot some of the output from `EpiGen`.
"""
function one_plot(results)
    f = Figure()
    ax = Axis(f[1,1], xlabel="time",ylabel="various units")
    lines!(ax,results[:,2],label="fitness")
    lines!(ax,results[:,5],label="genetic distance")
    lines!(ax,results[:,6],label="epigenetic distance")
    f[1, 2] = Legend(f, ax, framevisible = false)
    f
end

"""
    EpiGen(N=500,time=200)

Run the EpiGen model from [their paper](https://www.biorxiv.org/content/10.1101/637272v1) 

```
results=EpiGen(500,200)
one_plot(results)
```
"""
function EpiGen(N=500,time=200)

#N=1000 #the number of individuals in the population
#time=1000 #the number of generations to run the simulation
t_split = 40 #if patch='TRUE', t_split is the number of generations (40 in this case) the population stays in one environment or the other
radius=1 #distance from the hypersphere origin
nugen = 10 #number of genetic (HT) mutations for each generation
nuepi = 90 #number of epigenetic (LT) mutations for each generation
epiback = 0.0205 #epigenetic reversion rate
mlimitgen=2 #maximum effect a genetic mutation can have
mlimitepi = 0.3 #maximum effect an epigenetic mutation can have
epigenetics = true #switch to toggle effects of epigenetics on or off
patch = true #switch to determine if population switches between one environment and another. If TRUE, population goes in and out of selection environment at a frequency of t_split generations. If FALSE, only stays in the selection environment.
mode = "negative" #When patch = 'TRUE', this will run the model in two modes that determines the selection rules for environment 2: random or negative. 'negative' turns on randomly sampling individuals weighted by the reciprocal number of genetic mutations when in environment 2 (stabilizing selection). 'random' turns on randomly sampling indviduals in environment 2 with no weighting. Random sampling weighted by fitness always occurrs in environment 1 
runnumber = 1 #replicate number of model run

patch_intervals = fill(0.0, t_split) #Create cell array of matrices to store interval runs
patch_removed_intervals = fill(0.0, t_split) #Create cell array of matrices to store interval runs of removed individuals. These are individuals that were not sampled to go into the next generation.
results_mat = fill(0.0, t_split) #Create cell array of matrices to store results matrix runs
results_removed_mat = fill(0.0, t_split) #Create cell array of matrices to store results matrix runs of removed individuals
epi_mat = fill(0.0, t_split) #Create cell array of matrices to store epigenetic distribution run
gen_mat = fill(0.0, t_split) #Create cell array of matrices to store genetic distribution run
epi_removed_mat = fill(0.0, t_split) #Create cell array of matrices to store epigenetic distribution run
gen_removed_mat = fill(0.0, t_split) #Create cell array of matrices to store genetic distribution run

##

"""
Driver is a vector of 0's and 1's and has length = time in generations. The loop below will 
step through every element of Driver. When patch = 'TRUE' and Driver == 1, then select
= 1, and the population is randomly sampled to the next generation weighted
by fitness, i.e. is in environment 1. When patch = 'TRUE', mode = 'negative', and
Driver == 0, then the population is randomly sampled to the next generation weighted 
by the reciprocal of the number of genetic mutations (stabilizing selection),
i.e. is in environment 2.
"""

nt=time+t_split-1
Driver=fill(1,nt)
for q in t_split:t_split*2:time 
    Driver[q:q+t_split-1].=0
end

##

#Population information
population_co=fill(0.0,(N,25)) #create populations
population_co[:,1] .= 1:N #Index for individuals
population_co[:,2] .=1 #initial distance from the optimum is radius
population_co[:,3] .=exp.((-population_co[:,2].^2)./2) #Fitness - Using a Gaussian fitness function
#Column 4 is for distance traveled by genetic mutation
#Column 5 is for number of genetic mutations
#Column 6 is last time point of genetic mutations
#Column 7 is for distance traveled by epigenetic mutations
#Column 8 is for number of epigenetic mutations
#Column 9 is last time point of epigenetic mutations

#Information of individuals in the population that were removed from the main reproducing population
population_removed=zeros(N,25) ##create populations
#Column 1 #Index for individuals
#Column 2 #initial distance from the optimum is radius
#Column 3 #fitness function
#Column 4 is for distance traveled by genetic mutation
#Column 5 is for number of genetic mutations
#Column 6 is last time point of genetic mutations
#Column 7 is for distance traveled by epigenetic mutations
#Column 8 is for number of epigenetic mutations
#Column 9 is last time point of epigenetic mutations

#matrix to populate in loop
population_disgen=fill(0.0,(N,nt)) #distribution of genetic mutations of kept population (i.e. sampled to go onto next generation)
population_disepi=fill(0.0,(N,nt)) #distribution of epigenetic mutations of kept population
population_removed_disgen = fill(0.0,(N,nt)) #distribution of genetic mutations of removed population (i.e. not sampled to go onto next generation)
population_removed_disepi = fill(0.0,(N,nt)) #distribution of epigenetic mutations of removed population

##

#Initialize geometric model
#dimensions of spherical phenotypic space 
#Z-values have been solved manually, dimensions solved for are 5, 10, 15,
#20, 25, 30, 35, 40
#If some other value of n is needed, need to solve the integral for Z and add the value to Z.mat
zmat=[5 0.75; 10 1.1641; 15 1.4663; 20 1.7162; 25 1.9342; 30 2.1299; 35 2.3092; 40 2.4755];

#Store Z-value for 5 dimension. We use 5 dimensions for all model runs
Z5=zmat[1,2]

#the angle between the mutational vector and the vector running from the current phenotype to the origin
#Angles can take on values between -pi/2 to pi/2
phi_values=collect(range(-pi/2,pi/2,50000))

#Calculate the n = 5 probability density function for the random angle
prob_density5=Z5*cos.(phi_values) .^(zmat[1,1] - 2)

#Initialize Results of kept population
#Column 1 generations
#Column 2 mean population fitness
#Column 3 mean number of genetic mutations
#Column 4 mean number of epigenetic mutations
#Column 5 mean distance traveled by genetic mutations
#Column 6 mean distance traveled by epigenetic mutations
results_co=zeros(nt,28)
results_co[:,1] .= 1:nt
#colname={'time','meanfitness','nogenmut','noepimut','distgen','distepi','net_fi_carbon_gen','net_fi_p_gen','net_fi_i_gen','net_fi_carbon_epi','net_fi_p_epi','net_fi_i_epi'};
#colindex={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

#Initialize Results of removed population (i.e. individuals not sampled to
#next generation
results_removed=zeros(nt,28)
results_removed[:,1] .=1:nt

##


#select = 1 #Switch to randomly resample population either not weighted (0) or weighted (1) by fitness 

####
#Looping over time
for i in 1:length(Driver)

    #Check if patch is TRUE or FALSE. If patch_tf == 1, run model with only
    #1 type of selection (i.e. 1 environment) where sampling is weighted by
    #fitness. If patch_ft == 0, run the model for two types of selection (i.e. two
    #environments). When Driver = 1, sample weighted by fitness. When
    #Driver = 0, sample weighted by the reciprocal number of genetic
    #mutations
    patch_tf = !patch

    if patch_tf
        select = 1 #only sample with weighting as fitness
    elseif (!patch_tf) && (Driver[i] == 1)
        select = 1 #resample population with weight by fitness
    elseif (!patch_tf) && (Driver[i] == 0)
        select = 0 #resample population with other weighting parameters (e.g., by the reciprocal number of genetic mutations)
    end

    ####
    #Make genetic mutations, mutational supply is set to nugen/1 unit time
    #Distances from optimum are calculated individually
    
    seed=1000*Dates.second(now())+Dates.millisecond(now())
    Random.seed!(seed)
    
    #Vector of genetic mutations with each value being the magnitude of
    #their effect
    mut_effect=0 .+(mlimitgen*radius).*rand(1,nugen) 
    
    seed=1000*Dates.second(now())+Dates.millisecond(now())
    Random.seed!(2*seed)
    
    #Randomly sample 10 angles
    phi5=sample(phi_values, ProbabilityWeights(prob_density5), nugen; replace=true)    

    ###Assign mutations to individuals###
    
    index_mut=1:N
    
    seed=1000*Dates.second(now())+Dates.millisecond(now())
    Random.seed!(3*seed)
    
    #Randomly sample 10 individuals to genetically mutate without replacement (multiple hits are not allowed)
    index_mutated=sample(index_mut,nugen; replace=false) 

    #retrieve mutants from population 
    mutated_orgs =population_co[index_mutated,:]
    #mutated_orgs =dropdims(mutated_orgs,dims=1)
   
    ###calculate distance traveled after genetic mutation
    dist_travelled5 = mutated_orgs[:,2] - sqrt.(mutated_orgs[:,2].^2 + mut_effect[:].^2 .+ 2.0 *mutated_orgs[:,2].*mut_effect[:].*sin.(phi5[:]))

    #Calculate distance from the optimum phenotype
    dist5 = sqrt.(mutated_orgs[:,2].^2 + mut_effect[:].^2 + 2.0*sin.(phi5[:]).*mutated_orgs[:,2].*mut_effect[:])
    
    ####
    #Distance from the optimum phenotype based on oligo conditions;
    #oligo == 1 ; n = 5 dimensions of the hypersphere
        
    #store distance from optimum for mutated individuals
    mutated_orgs[:,2] .= dist5
    
    #rename distance_travelled5 vector to dist_travelled_gen
    dist_travelled_gen = dist_travelled5
        
    #NOT USED IN THIS VERSION OF MODEL CALCULATIONS - IGNORE
    #Using law of sines, calculate alternative angle.
    fi5 = asin.((sin.(phi5[:] .+ (pi/2)).*(mut_effect[:]))./dist5)

    #Insert fi5
    mutated_orgs[:,10] .= mutated_orgs[:,10] .+ fi5

    #Save mean distance traveled toward optimum by genetic mutation
    results_co[i,25] = mean(mutated_orgs[:,2])
    #Save std of distance traveled toward optimum by genetic mutation
    results_co[i,26] = std(mutated_orgs[:,2])
    #Calculates fitness for the mutants after genetic mutation
    mutated_orgs[:,3] .= exp.((-mutated_orgs[:,2].^2)/2)
    ###Cumulative distance travelled by genetic mutations
    mutated_orgs[:,4] .= mutated_orgs[:,4] + dist_travelled_gen
    #Increase the number of genetic mutations by one
    mutated_orgs[:,5] .= mutated_orgs[:,5] .+ 1
    #save distance travelled after genetic mutations
    mutated_orgs[:,6] .= dist_travelled_gen
    #Storing the mutated individuals
    population_co[index_mutated,:].=mutated_orgs[:,:]
    #Storing distribution of genetic mutations
    for g=1:length(dist_travelled_gen)
        population_disgen[index_mutated[g], i] = dist_travelled_gen[g]
    end
   
    #If the model run includes epigenetic effects, run the following
    tf = epigenetics

    if tf == 1
        
        #Make epigenetic mutations, mutational supply is set to Nu.epi / 1 unit of time
        
        seed=1000*Dates.second(now())+Dates.millisecond(now())
        Random.seed!(4*seed)
            
        #Mutational effects of epigenetics
        mut_effect_epi=(mlimitepi*radius).*rand(1, nuepi)
        
        #Random angle for epigenetic mutations
        
        seed=1000*Dates.second(now())+Dates.millisecond(now())
        Random.seed!(5*seed)
            
        #Randomly draw random angles for epigenetic mutations
        phi_epi5=sample(phi_values, ProbabilityWeights(prob_density5), nuepi; replace=true)
        
        seed=1000*Dates.second(now())+Dates.millisecond(now())
        Random.seed!(6*seed)
            
        #retrieve epimutants from population w/out replacement
        index_mutated_epi = sample(index_mut, nuepi, replace=false)
        mutated_orgs_epi = population_co[index_mutated_epi, :]
        
        ###distance traveled after epimutations
        dist_travelled_epi5 = mutated_orgs_epi[:,2] - sqrt.(mutated_orgs_epi[:,2].^2 + mut_effect_epi[:].^2 + 2.0*mutated_orgs_epi[:,2].*mut_effect_epi[:].*sin.(phi_epi5[:]))

        #Distance from the optimum phenotype epimutations
        dist_epi5 = sqrt.(mutated_orgs_epi[:,2].^2 + mut_effect_epi[:].^2 + 2.0*sin.(phi_epi5[:]).*mutated_orgs_epi[:,2].*mut_effect_epi[:])

        #Distance from the optimum phenotype based on oligo conditions;

        #Store distance from the optimum phenotype 
        mutated_orgs_epi[:,2] .= dist_epi5 

        #rename dist_travelled_epi5 to dist_travelled_e
        dist_travelled_e = dist_travelled_epi5;
        
        #Save mean distance travelled from optimum by epigenetic mutation
        results_co[i,27] = mean(mutated_orgs_epi[:,2])
    
        #Save std of distance travelled from optimum by epigenetic mutation
        results_co[i,28] = std(mutated_orgs_epi[:,2])
        
        #Calculate fitness for mutants
        mutated_orgs_epi[:,3] .= exp.((-mutated_orgs_epi[:,2].^2)/2)
        
        #Cumulative distance travelled by epigenetic mutations
        mutated_orgs_epi[:,7] .= mutated_orgs_epi[:,7] + dist_travelled_e
        
        #Increase the number of epigenetic mutations by one
        mutated_orgs_epi[:,8] .= mutated_orgs_epi[:,8] .+ 1
        
        #Store distance travelled
        mutated_orgs_epi[:,9] .= dist_travelled_e

        #CALCULATE ALTERNATIVE ANGLE - IGNORE THIS METRIC - NOT USED IN THIS MODEL VERSION
        #Using law of sines
        #fi_epi5 = asin.((sin.(phi_epi5[:] /+ pi/2).*(mut_effect_epi[:]))./dist_epi5[:])

        #Insert fi5 for epi
        #mutated_orgs_epi[:,18] = mutated_orgs_epi[:,18] + fi_epi5
        
        #Storing mutated individuals
        population_co[index_mutated_epi,:] .= mutated_orgs_epi
        
        #Distribution of epigenetic mutations        
        for e=1:length(dist_travelled_e)
            population_disepi[index_mutated_epi[e], i] = dist_travelled_e[e]
        end
    
    end
    
    
    #As of now, multiple hits of either genetic and epigenetic mutations are not possible. 
    #However, a single individual can have two hits if one is genetic and one epigenetic mutation. 

    ######
    #Population regulation and reproduction

    #Soft selection is acting -> population size stays the same
    #Sampling individuals to the next time point with replacement
    
    ofsp_index = 1:N;
    
    if select == 0 #if true, population is in environment 2
        
        mode_type = (mode=="random") #If true, sample without weighting
        #Sampling the next generation, no weighting
        if mode_type == 1
        
        seed=1000*Dates.second(now())+Dates.millisecond(now())
        Random.seed!(7*seed)
                
        sampled_index = sample(ofsp_index, N, replace=true)
        else
        #Sampling the next generation, negative weighting (i.e. the
        #reciprocal number of genetic mutations)
        reciprocal_weights = 1 ./ (0.1 .+ population_co[:,5])
        
        seed=1000*Dates.second(now())+Dates.millisecond(now())
        Random.seed!(8*seed)
            
        sampled_index = sample(ofsp_index, ProbabilityWeights(reciprocal_weights), N; replace=true)

        end
        
    elseif select == 1 #if true, population is in environment 1
        
        #Sampling the next generation, weighting with fitness
        seed=1000*Dates.second(now())+Dates.millisecond(now())
        Random.seed!(9*seed)

        sampled_index = sample(ofsp_index, ProbabilityWeights(population_co[:,3]), N; replace=true)
       
    end

    #Population after reproduction and regulation
    population_co .= population_co[sampled_index,:]
    population_disgen .= population_disgen[sampled_index,:]
    population_disepi .= population_disepi[sampled_index,:]
    
    #Removed individuals from population
    pop_removed_index = setdiff(ofsp_index, sampled_index)
    population_removed = population_co[pop_removed_index,:]
    
    population_removed_disgen = population_disgen[pop_removed_index,:]
    population_removed_disepi = population_disepi[pop_removed_index,:]
    
    ######Store results removed population
   
    #Mean fitness of the population
    results_removed[i,2]= mean(population_removed[:,3]) 
    #Store the mean number of genetic mutations present in an individual
    results_removed[i,3]= mean(population_removed[:,5]) 
    #Store the standard deviation of mean genetic mutations present in an individual
    results_removed[i,23]= std(population_removed[:,5])
    #Store the mean number of epigenetic mutations present in an individual
    results_removed[i,4]= mean(population_removed[:,8])
    #Store the standard deviation of mean epigenetic mutations present in an individual
    results_removed[i,24]= std(population_removed[:,8])
    #Store the cumulative distance travelled towards the optimum by genetic mutations
    results_removed[i,5]= mean(population_removed[:,4])
    #Store the cumulative distance travelled towards the optimum by epigenetic mutations
    results_removed[i,6]= mean(population_removed[:,7])
    
    
    ######Store results sampled (non-removed) population

    #Mean fitness of the population
    results_co[i,2]= mean(population_co[:,3])
    #Store the mean number of genetic mutations present in an individual
    results_co[i,3]= mean(population_co[:,5]) 
    #Store the standard deviation of mean genetic mutations present in an individual
    results_co[i,23]= std(population_co[:,5])
    #Store the mean number of epigenetic mutations present in an individual
    results_co[i,4]= mean(population_co[:,8])
    #Store the standard deviation of mean epigenetic mutations present in an individual
    results_co[i,24]= std(population_co[:,8])
    #Store the cumulative distance travelled towards the optimum by genetic mutations
    results_co[i,5]= mean(population_co[:,4])
    #Store the cumulative distance travelled towards the optimum by epigenetic mutations
    results_co[i,6]= mean(population_co[:,7])
    #Store mean genetic mutation fi5
    results_co[i,7] = mean(population_co[:,10])
    #Store mean epigenetic mutation fi5
    results_co[i,15] = mean(population_co[:,18])
              
    ######
    #Loss of epigenetic effects
    #Next epigenetic changes are lost at a given rate 
    #(This needs to be slower than the forward mutation rate of epigenetic effects)
    #Epigenetic changes revert independently (always at a probability of: epi.backmutation in each generation)
    
    #Number of epigenetic mutations present
    no_epi_mut = sum(population_co[:,8])

    if (no_epi_mut>0)
        #Check which individuals have epigenetic mutation
        epi_p_index = findall(population_co[:,8] .> 0)
        epi_present = population_co[epi_p_index,:]
        epi_present_gen = population_disgen[epi_p_index,:];
        epi_present_epi = population_disepi[epi_p_index,:];
        
        #Generate backmutations probabilistically
        n=Int(sum(epi_present[:,8]))
        if n > 0
            p=Binomial(1,epiback)
            epi_backmuts = rand(p,n)
        else
            epi_backmuts = fill(0.0,n)
        end

        #Check for each backmutation whether reversion occurs or not
        ###Mapping from individual mutations to individuals
        
        #Construct a helping vector with as many factors as there are inds with the equal 
        #length of mutation vector
        
        help_vec = vcat(fill.( 1:length(epi_present[:,8]) ,  Int.(epi_present[:,8]) )...)
                
        ####
        #Sum by each factor
        ##Same function as "tapply" in R
        
        df = DataFrame(a=epi_backmuts, b=help_vec)
        gdf = groupby(df, :b)
        tmp = combine(gdf, :a => sum)
        list_of_reversions = tmp.a_sum
                
        #Logical vector if backmutation occurs in a given individual (TRUE) or not
        epi_back_mut_index = findall(list_of_reversions .> 0)
        
        #Take only those individuals that had backmutations
        epi_back_mut_ind = epi_present[epi_back_mut_index,:]
        epi_back_mut_ind_gen = epi_present_gen[epi_back_mut_index,:]
        epi_back_mut_ind_epi = epi_present_epi[epi_back_mut_index,:]
        
        #How many backmutations occurred in each individual that had backmutations
        number_of_reversions = list_of_reversions[epi_back_mut_index]
        
        #if any backmutations occurred
        if sum(list_of_reversions.>0) > 0
            
            seed=1000*Dates.second(now())+Dates.millisecond(now())
            Random.seed!(10*seed)
                    
            for b=1:length(number_of_reversions)
                
                lost_mutation=zeros(1,number_of_reversions[b])
                
                current_org = [epi_back_mut_ind_epi[b]]
                current_org = current_org[findall((!).(current_org .== 0))]
                #randomly sample epi effects from the distribution
                if length(current_org) >= number_of_reversions[b]
                    lost_mutation[1,:]=sample(current_org, number_of_reversions[b],replace=false)
                end
                       
                #Now individuals need to lose mutation
                if !isempty(lost_mutation)
                    
                    #loop through all epimutations
                    for count=1:length(lost_mutation) 
                        #ensure current org contains epimutations to be
                        #lost
                        logic = sum(current_org .== lost_mutation[1,count]) 
                        
                        if logic == 1
                            #gather epimutation values to be lost
                            r = epi_back_mut_ind_epi[b,:]

                            #find indicies of lost mutations
                            idx2 = findall(r.==lost_mutation[1,count])
                                 
                            #if only 1 mutation is found,
                            #set to zero
                            if length(idx2) == 1 
                                epi_back_mut_ind_epi[b,idx2] .= 0
                            
                            #if more than 1 epimutation, loop through the
                            #number of epimutations to lose
                            elseif length(idx2) > 1
                                [epi_back_mut_ind_epi[b,idx2[l]] = 0 for l=1:length(idx2)]
                            end
                        end
                    end
                end
                             
               #Subtract the effect of the mutations from distance travelled by epigenetics
               epi_back_mut_ind[b,7]= epi_back_mut_ind[b,7]-sum(lost_mutation[1,:])
               #Add the effect of the mutation from distance to optimum
               epi_back_mut_ind[b,2]= epi_back_mut_ind[b,2]+sum(lost_mutation[1,:])
               #Subtract the number of epigenetic mutations
               epi_back_mut_ind[b,8]=epi_back_mut_ind[b,8]-number_of_reversions[b]               
               #Recalculate fitness for backmutated individuals
               epi_back_mut_ind[b,3] = exp((-epi_back_mut_ind[b,2].^2)/2)
            end
                          
            #Replace backmutated individuals in the population_co dataframe
            epi_present[epi_back_mut_index,:]=epi_back_mut_ind
            population_co[epi_p_index,:].=epi_present
            
            #Replace the epigenetic mutational distributions of the 
            #backmutated individuals in the population_disepi dataframe
            epi_present_epi[epi_back_mut_index,:]=epi_back_mut_ind_epi
            population_disepi[epi_p_index,:]=epi_present_epi
        end
    end    
end

return results_co
end #function EpiGen()

