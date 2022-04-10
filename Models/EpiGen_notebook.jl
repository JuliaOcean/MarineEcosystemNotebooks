### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5aa7f3b3-efd6-47f3-adcf-a69553bef200
begin
	using StatsBase, Distributions, DataFrames, Dates, Random, UnPack
	
	using CairoMakie, PlutoUI

	"Done with packages."
end

# ╔═╡ d5e545be-9e40-48ce-9ccc-c694ae621d4e
PlutoUI.TableOfContents()

# ╔═╡ 1caab1b2-b815-11ec-3f7b-2945e8c38199
md"""

!!! credits
    This code, and its embedded documentation, was ported to Julia from the original Matlab implementation of [EpiGen](https://github.com/LevineLab/EpiGen). Please refer to [their paper](https://www.biorxiv.org/content/10.1101/637272v1) if you use this model. And if in doubt then please use the original code.
"""

# ╔═╡ fb6f8673-4ca9-4417-af77-977ec31947fd
md"""## Model Run Example"""

# ╔═╡ 3d66d246-057c-4564-9d34-b18028ad29e6
begin
	parameters=(
	N=[500], #the number of individuals in the population
	time=[500], #the number of generations to run the simulation
	t_split = [40], #if patch='TRUE', t_split is the number of generations (40 in this case) the population stays in one environment or the other
	radius=[1], #distance from the hypersphere origin
	nugen = [10], #number of genetic (HT) mutations for each generation
	nuepi = [90], #number of epigenetic (LT) mutations for each generation
	epiback = [0.0205], #epigenetic reversion rate
	mlimitgen=[2], #maximum effect a genetic mutation can have
	mlimitepi = [0.3], #maximum effect an epigenetic mutation can have
	epigenetics = [true], #switch to toggle effects of epigenetics on or off
	patch = [true], #switch to determine if population switches between one environment and another. If TRUE, population goes in and out of selection environment at a frequency of t_split generations. If FALSE, only stays in the selection environment.
	mode = ["negative"], #When patch = 'TRUE', this will run the model in two modes that determines the selection rules for environment 2: random or negative. 'negative' turns on randomly sampling individuals weighted by the reciprocal number of genetic mutations when in environment 2 (stabilizing selection). 'random' turns on randomly sampling indviduals in environment 2 with no weighting. Random sampling weighted by fitness always occurrs in environment 1 
	runnumber = [1] #replicate number of model run
	)
	
	"Default parameters set up"
end

# ╔═╡ b8ab5a19-3f5f-45a3-a70d-d0f607a6555d
begin
	𝑷=DataFrame(group=String[],name=String[],default=Float64[],factors=Array[],
		long_name=String[],unit=String[])
	push!(𝑷,("main","N",1000,[0.2, 0.5, 1.0, 2.0, 3.0, 5.0],"number of individuals","unitless"))
	push!(𝑷,("main","time",200,[0.2, 0.5, 1.0, 2.0, 3.0, 5.0],"number of generations","unitless"))
	push!(𝑷,("main","t_split",40,[0.2, 0.5, 1.0, 2.0, 3.0, 5.0],"environment duration","unitless"))
	push!(𝑷,("main","epiback",1.0,[0.0205, 0.2, 0.5],"epigenetic reversion rate","unitless"))
	push!(𝑷,("main","epigenetics",true,[true, false],"epigenetics on or off","unitless"))
	push!(𝑷,("main","nugen",10,[0.5, 1.0, 2.0, 5.0],"genetic mutations / generation","unitless"))
	push!(𝑷,("main","nuepi",90,[0.5, 1.0, 2.0, 5.0],"epigenetic mutations / generation","unitless"))
	
	𝑉=[𝑷.default[i]*𝑷.factors[i] for i in 1:length(𝑷.default)]
	𝑷
	
	md"""### Modify Parameters

	Parameter name | Value | unit
	----|----|----
	$(𝑷.long_name[1]) | $(@bind 𝑄_1 Select(𝑉[1]; default=𝑷.default[1]))  |  $(𝑷.unit[1])
	$(𝑷.long_name[2]) | $(@bind 𝑄_2 Select(𝑉[2]; default=𝑷.default[2]))  |  $(𝑷.unit[2])
	$(𝑷.long_name[3]) | $(@bind 𝑄_3 Select(𝑉[3]; default=𝑷.default[3]))  |  $(𝑷.unit[3])
	----|----|----
	$(𝑷.long_name[6]) | $(@bind 𝑄_6 Select(𝑉[6]; default=𝑷.default[6]))  |  $(𝑷.unit[6])
	$(𝑷.long_name[7]) | $(@bind 𝑄_7 Select(𝑉[7]; default=𝑷.default[7]))  |  $(𝑷.unit[7])
	$(𝑷.long_name[4]) | $(@bind 𝑄_4 Select(𝑉[4]; default=𝑷.default[4]))  |  $(𝑷.unit[4])
	$(𝑷.long_name[5]) | $(@bind 𝑄_5 Select(𝑉[5]; default=𝑷.default[5]))  |  $(𝑷.unit[5])
	----|----|----

	### Update Model Run
	
	$(@bind update_param PlutoUI.Button("Run Model"))
	"""
end

# ╔═╡ ea6a560b-e466-4727-bf2f-bdc5d220d3fc
md"""## Packages and Functions"""

# ╔═╡ 2977ac5a-8408-47b3-bd5b-03910486130f
"""
    EpiGen(parameters,storage)

Main function that runs the model over `parameters.Driver`.

```
storage=setup_storage(parameters)
results=EpiGen(parameters,storage)
one_plot(results)
```
"""
function EpiGen(parameters,storage)	

	N=parameters.N[1]
	time=parameters.time[1]
	t_split=parameters.t_split[1]

	radius=parameters.radius[1]
	nugen=parameters.nugen[1]
	nuepi=parameters.nuepi[1]
	epiback=parameters.epiback[1]
	mlimitgen=parameters.mlimitgen[1]
	mlimitepi=parameters.mlimitepi[1]	

	epigenetics=parameters.epigenetics[1]
	patch=parameters.patch[1]
	mode=parameters.mode[1]
	runnumber=parameters.runnumber[1]

	@unpack patch_intervals, patch_removed_intervals = storage
	@unpack results_mat, results_removed_mat = storage
	@unpack epi_mat, gen_mat, epi_removed_mat, gen_removed_mat = storage
	@unpack patch_intervals, patch_removed_intervals = storage

	@unpack Driver, population_co, population_removed = storage
	@unpack population_disgen, population_disepi = storage
	@unpack population_removed_disgen, population_removed_disepi = storage
	@unpack zmat, prob_density5, phi_values, results_co, results_removed = storage
	
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

# ╔═╡ 944a2695-642b-4a81-930f-08351f5a5290
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
	#xlims!(ax, [0, 1000])
	f
end

# ╔═╡ 811043d0-7b74-4f2c-b34d-fefe2156a657
"""
    setup_storage(parameters)

Set up and initialize arrays etc for storing state variables and results.
"""
function setup_storage(parameters)
	@unpack N, time, t_split = parameters # equivalent to: a,b = pa.a,pa.b
	
	ni=N[1]
	nt_0=time[1]
	nt_s=t_split[1]
	
	patch_intervals = fill(0.0, nt_s) #to store interval runs
	patch_removed_intervals = fill(0.0, nt_s) #to store interval runs of removed individuals. These are individuals that were not sampled to go into the next generation.
	results_mat = fill(0.0, nt_s) #to store results matrix runs
	results_removed_mat = fill(0.0, nt_s) #to store results matrix runs of removed individuals
	epi_mat = fill(0.0, nt_s) #to store epigenetic distribution run
	gen_mat = fill(0.0, nt_s) #to store genetic distribution run
	epi_removed_mat = fill(0.0, nt_s) #to store epigenetic distribution run
	gen_removed_mat = fill(0.0, nt_s) #to store genetic distribution run
	
	"""
	Driver is a vector of 0's and 1's and has length = time in generations. The loop below will 
	step through every element of Driver. When patch = 'TRUE' and Driver == 1, then select
	= 1, and the population is randomly sampled to the next generation weighted
	by fitness, i.e. is in environment 1. When patch = 'TRUE', mode = 'negative', and
	Driver == 0, then the population is randomly sampled to the next generation weighted 
	by the reciprocal of the number of genetic mutations (stabilizing selection),
	i.e. is in environment 2.
	"""
	nt=nt_0+nt_s-1
	Driver=fill(1,nt)	 
	[Driver[q:q+nt_s-1].=0 for q in nt_s:nt_s*2:nt_0]
	
	#Population information
	population_co=fill(0.0,(ni,25)) #create populations
	population_co[:,1] .= 1:ni #Index for individuals
	population_co[:,2] .=1 #initial distance from the optimum is radius
	population_co[:,3] .=exp.((-population_co[:,2].^2)./2) #Fitness - Using a Gaussian fitness function
	#Column 4 is for distance traveled by genetic mutation
	#Column 5 is for number of genetic mutations
	#Column 6 is last time point of genetic mutations
	#Column 7 is for distance traveled by epigenetic mutations
	#Column 8 is for number of epigenetic mutations
	#Column 9 is last time point of epigenetic mutations
	
	#Information of individuals in the population that were removed from the main reproducing population
	population_removed=zeros(ni,25) ##create populations
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
	population_disgen=fill(0.0,(ni,nt)) #distribution of genetic mutations of kept population (i.e. sampled to go onto next generation)
	population_disepi=fill(0.0,(ni,nt)) #distribution of epigenetic mutations of kept population
	population_removed_disgen = fill(0.0,(ni,nt)) #distribution of genetic mutations of removed population (i.e. not sampled to go onto next generation)
	population_removed_disepi = fill(0.0,(ni,nt)) #distribution of epigenetic mutations of removed population

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
	
	return (patch_intervals=patch_intervals,patch_removed_intervals=patch_removed_intervals,
	results_mat=results_mat,results_removed_mat=results_removed_mat,
	epi_mat=epi_mat,gen_mat=gen_mat,epi_removed_mat=epi_removed_mat,gen_removed_mat=gen_removed_mat,
	Driver=Driver, population_co=population_co, population_removed=population_removed,
	population_disgen=population_disgen, population_disepi=population_disepi,
	population_removed_disgen=population_removed_disgen, population_removed_disepi=population_removed_disepi,
	zmat=zmat, prob_density5=prob_density5, phi_values=phi_values,
	results_co=results_co, results_removed=results_removed
	)	
end

# ╔═╡ 59d02328-34b6-4304-822f-04fa2fde457e
begin
	update_param

	#modify parameter values within nml
	parameters.N[1]=𝑄_1
	parameters.time[1]=𝑄_2
	parameters.t_split[1]=𝑄_3
	parameters.epiback[1]=𝑄_4
	parameters.epigenetics[1]=𝑄_5
	parameters.nugen[1]=𝑄_6
	parameters.nuepi[1]=𝑄_7
		
	#initialize storage space
	storage=setup_storage(parameters)

	#run model
	results=EpiGen(parameters,storage)
		
end

# ╔═╡ 18b6392d-2559-4e24-9620-0881536717fb
one_plot(results)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
UnPack = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"

[compat]
CairoMakie = "~0.7.5"
DataFrames = "~1.3.2"
Distributions = "~0.25.53"
PlutoUI = "~0.7.38"
StatsBase = "~0.33.16"
UnPack = "~1.0.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "c933ce606f6535a7c7b98e1d86d5d1014f730596"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "5.0.7"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "StaticArrays"]
git-tree-sha1 = "4a0de4f5aa2d5d27a1efa293aeabb1a081e46b2b"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.7.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "5a4168170ede913a2cd679e53c2123cb4b889795"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.53"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d064b0340db45d48893e7604ec95e7a2dc9da904"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.5.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "169c3dc5acae08835a573a8a3e25c62f689f8b5c"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.6.5"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageIO]]
deps = ["FileIO", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "464bdef044df52e6436f8c018bea2d48c40bb27b"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.1"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "bcf640979ee55b652f3b01650444eb7bbe3ea837"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "63de3b8a5c1f764e4e3a036c7752a632b4f0b8d1"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.16.6"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "c5fb1bfac781db766f9e4aef96adc19a729bc9b2"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.2.1"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "70e733037bbf02d691e78f95171a1fa08cdc6332"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.2.1"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e8185b83b9fc56eb6456200e873ce598ebc7f262"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.7"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "eb4dbb8139f6125471aa3da98fb70f02dc58e49c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.14"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a121dfbba67c94a5bec9dde613c3d0cbcf3a12b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.3+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "28ef6c7ce353f0b35d0df0d5930e0d072c1f5b9b"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMD]]
git-tree-sha1 = "7dbc15af7ed5f751a82bf3ed37757adf76c32402"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.1"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "9cc2955f2a254b18be655a4ee70bc4031b2b189e"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "87e9954dfa33fd145694e42337bdd3d5b07021a6"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "aaa19086bc282630d82f818456bc40b4d314307d"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.4"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─d5e545be-9e40-48ce-9ccc-c694ae621d4e
# ╟─1caab1b2-b815-11ec-3f7b-2945e8c38199
# ╟─fb6f8673-4ca9-4417-af77-977ec31947fd
# ╟─3d66d246-057c-4564-9d34-b18028ad29e6
# ╟─b8ab5a19-3f5f-45a3-a70d-d0f607a6555d
# ╟─18b6392d-2559-4e24-9620-0881536717fb
# ╟─59d02328-34b6-4304-822f-04fa2fde457e
# ╟─ea6a560b-e466-4727-bf2f-bdc5d220d3fc
# ╟─5aa7f3b3-efd6-47f3-adcf-a69553bef200
# ╟─2977ac5a-8408-47b3-bd5b-03910486130f
# ╟─944a2695-642b-4a81-930f-08351f5a5290
# ╟─811043d0-7b74-4f2c-b34d-fefe2156a657
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
