#### Bioturbation in Julia ####
include("bioturbationFunctions.jl")

## Objects and parameters 
bd = 1500               # Bulk density [kg m-3]
nlayers = 200           # Number of soil layers
layer_thickness = 0.01  # Initial layer thickness [m]
bleaching_depth = 0.005 # Thickness of surface layer where bleaching occurs [m]
soildepth = 2           # Thickness of the entire soil profile [m]
grains_per_layer = 150  # Initial number of layers per soil layer
BT_pot = 10             # Potential bioturbation rate [kg m-2]
dd = 1                  # Depth parameter [m-1]
dd_exch = 10            # Depth parameter for subsurface mixing [m-1]
ntime = 10000           # Amount of simulation years
soil = zeros(Float64, nlayers, 4)   # Array of the soil profile. Index 1: thickness [m]. Index 2: mass [kg]. Index 3: midpoint depth [m]. Index 4: cumulative depth [m].
ages = fill(Int[],nlayers)          # Array of arrays of OSL age particles


# Simulations

# A run with bioturbation by mounding
initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)

for t in ProgressBar(1:ntime)
    
    BT_mounding(BT_pot, "grd")

    update_layers()
    update_OSL()
end

# Update soil layer properties
depth = 0
for l in 1:nlayers
    depth += soil[l,1]/2
    soil[l,3] = depth
    depth += soil[l,1]/2
end
soil[:,4]=cumsum(soil[:,1])

plot_age_depth_profile(ages, soil,5,3)

# Save output as .jld2 files
save(string(wd,"output_soil.jld2"), "soil", soil)
save(string(wd,"output_ages.jld2"), "ages", ages)


# Exp. 2: Different ratios of mounding and mixing
depth_function = "grd"
ntime = 10000
dd_exch = 10
BT_pot = 10
ratios = [0,.25,.5,.75,.9,.95,.99,1]

# Create final output arrays
output_soil = zeros(length(ratios), nlayers,4)
output_ages = fill(Int[],length(ratios), nlayers)
runtimes = zeros(length(ratios))

# Run simulations
for i1 in 1:length(ratios)
    t1 = time()
    println("Run $i1 of $(length(ratios))")
    println(string("Ratio: $(ratios[i1])"))
    initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)

    for t in 1:ntime
        if t%100 == 0
            println("t ",t)
        end
        BT_mixing(BT_pot * ratios[i1], depth_function)
        BT_mounding(BT_pot * (1-ratios[i1]), depth_function)
        update_layers()
        update_OSL()
        update_layers()
    end

    depth = 0
    for l in 1:nlayers
        depth += soil[l,1]/2
        soil[l,3] = depth
        depth += soil[l,1]/2
    end
    soil[:,4]=cumsum(soil[:,1])

    # Write output to final output arrays
    output_soil[i1,:,:] = soil
    output_ages[i1,:] = ages

    runtimes[i1] = time() - t1
end

runtimes_min = runtimes / 60
runtimes_min[6]
wd = "C:\\Simulations\\Bioturbation\\Julia\\"
save(string(wd,"exp2_output_soil.jld2"), "output_soil", output_soil)
save(string(wd,"exp2_output_ages.jld2"), "output_ages", output_ages)


# Exp 3: Upheaval
mixing_depth = [.1,.25,.5,.75]
frequency = vcat([1:9,(1:9)*10,(1:9)*100,(1:5)*1000]...)
ntime = 10000

# Create final output arrays
output_soil = zeros(length(mixing_depth), length(frequency), nlayers,4)
output_ages = fill(Int[],length(mixing_depth), length(frequency), nlayers)
runtimes = zeros(length(mixing_depth) * length(frequency))

ind=0
for md in 1:length(mixing_depth)
    for freq in ProgressBar(1:length(frequency))
        t1=time()
        ind+=1

        initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)

        for t in 1:ntime
            BT_upheaval(md, freq, t)

            update_OSL()
            update_layers
        end

        depth = 0
        for l in 1:nlayers
            depth += soil[l,1]/2
            soil[l,3] = depth
            depth += soil[l,1]/2
        end
        soil[:,4]=cumsum(soil[:,1])

        # Write output to final output arrays
        output_soil[md,freq,:,:] = soil
        output_ages[md,freq,:] = ages

        runtimes[ind] = time() - t1
    end
end
runtimes_min = runtimes / 60

wd = "C:\\Simulations\\Bioturbation\\Julia\\"
save(string(wd,"exp3_output_soil.jld2"), "output_soil", output_soil)
save(string(wd,"exp3_output_ages.jld2"), "output_ages", output_ages)


# Exp 4. All processes together
upheaval_frequencies = [100,1000,2000,5000]
ntime = 10000
dd_exch = 10
output_soil = zeros(length(upheaval_frequencies), nlayers,4)
output_ages = fill(Int[],length(upheaval_frequencies),  nlayers)
runtimes = zeros(length(upheaval_frequencies))

for freq in 1:length(upheaval_frequencies)
    t1=time()

    initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)


    for t in ProgressBar(1:ntime)
        BT_upheaval(0.75, freq, t)
        BT_mounding(10 * 0.1, "grd")
        BT_mixing(10 * 0.9, "grd")
        update_layers()
        
        update_OSL()
        update_layers()
    end

    depth = 0
    for l in 1:nlayers
        depth += soil[l,1]/2
        soil[l,3] = depth
        depth += soil[l,1]/2
    end
    soil[:,4]=cumsum(soil[:,1])

    # Write output to final output arrays
    output_soil[freq,:,:] = soil
    output_ages[freq,:] = ages

    runtimes[freq] = time() - t1
end

wd = "C:\\Simulations\\Bioturbation\\Julia\\"
save(string(wd,"exp4_output_soil.jld2"), "output_soil", output_soil)
save(string(wd,"exp4_output_ages.jld2"), "output_ages", output_ages)




# Exp. 5: Repeatability simulations
depth_function = "grd"
ntime = 10000
dd_exch = 10
BT_pot = 10
ratios = [1,1,1,1,1,0,0,0,0,0] # 5 repetitions for each process

# Create final output arrays
output_soil = zeros(length(ratios), nlayers,4)
output_ages = fill(Int[],length(ratios), nlayers)
runtimes = zeros(length(ratios))

# Run simulations
for i1 in 1:length(ratios)
    t1 = time()
    println("Run $i1 of $(length(ratios))")
    println(string("Ratio: $(ratios[i1])"))
    initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)

    for t in ProgressBar(1:ntime)
        if t%100 == 0
            println("t ",t)
        end
        BT_mixing(BT_pot * ratios[i1], depth_function)
        BT_mounding(BT_pot * (1-ratios[i1]), depth_function)
        update_layers()
        update_OSL()
        update_layers()
    end

    depth = 0
    for l in 1:nlayers
        depth += soil[l,1]/2
        soil[l,3] = depth
        depth += soil[l,1]/2
    end
    soil[:,4]=cumsum(soil[:,1])

    # Write output to final output arrays
    output_soil[i1,:,:] = soil
    output_ages[i1,:] = ages

    runtimes[i1] = time() - t1
end

runtimes_min = runtimes / 60
runtimes_min[6]
wd = "C:\\Simulations\\Bioturbation\\Julia\\"
save(string(wd,"exp5_output_soil.jld2"), "output_soil", output_soil)
save(string(wd,"exp5_output_ages.jld2"), "output_ages", output_ages)



# Exp 6. Calibrate termite profile
calib_results = hcat(0, 0, 0, -1)
termites = CSV.read("C:\\Users\\Marijn\\sciebo\\Research\\Bioturbation\\Data\\Termites_Ghana\\De_termites_Ghana_complete.csv", DataFrame)
depths = unique(termites[:,"depth"])
ntime = 10
BT_pot = 10
dd_ref = 1

# Fractions to multiply input to
cal_depthfunction = ["grd"]
cal_BT_pot = [.25,.5,1,2,4]
cal_dd_pot = [.25,.5,1,2,4]




for i1 in ProgressBar(1:length(cal_depthfunction))
    for i2 in ProgressBar(1:length(cal_BT_pot))
        for i3 in ProgressBar(1:length(cal_dd))
            initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)

            # Simulations
            for t in ProgressBar(1:ntime)
                BT_mounding(BT_pot * cal_BT_pot[i2], cal_depthfunction[i1], cal_dd[i3])
                update_layers()
                update_OSL()
                update_layers()
            end

            # Calibration
            total_agreement = 0
            for d in 1:length(depths)

                it_layer = 1
                depth_ref = 0
                depth_dist = 9999

                for l in 1:nlayers
                    depth_ref += soil[l,1]/2
                    if abs(depth_ref-depths[d])<depth_dist
                        depth_dist = abs(depth_ref-depths[d])
                        it_layer +=1
                    else
                        break
                    end
                    depth_ref += soil[l,1]/2
                end

                termites_sub = termites[termites[:,"depth"].==depths[d],:]
                termites_ages = termites_sub[:,"De"]./termites_sub[:,"DR"]*1000  
                model_ages = copy(ages[it_layer])
                # filter!(e->e≠ntime,model_ages)
                total_agreement += agreement_distributions(termites_ages, model_ages)
            end

            calib_results = vcat(i1, BT_pot* cal_BT_pot[i2], dd_ref, total_agreement)
        end
    end
end




ind=findall(x->x==maximum(cal_results), cal_results)
println("Best result: df: $(cal_depthfunction[ind[1]]), BTpot_corr: $(cal_BT_pot[ind[2]]), dd_corr: $(cal_dd[ind[3]])")

plot()




# Exp XX. delta dose rate
De = fill(Float64[],nlayers)
DR = zeros(Float64,nlayers,1)
initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)
for t in ProgressBar(1:ntime)
    BT_mounding(10,"grd")
    update_layers()
    update_OSL()
    update_layers()
end


ages_De = copy(De)
for l in 1:nlayers
    ages_De[l] = ages_De[l] ./DR[l]
end





using ColorSchemes
cg = cgrad(:viridis,range(0,stop=1,length=nlayers))



plot(0,xlims = [0, 10000], ylims = [0, 10000])
cols = cgrad(:autumn1, nlayers)

xlabel!("Real age")
ylabel!("De age")
for l in 1:nlayers
    plot!(ages[l], ages_De[l],c=cols[l])
end
Plots.abline!(1,0, lwd = 3, color="black")




plot(ages_De, ages)

ages[12]
ages[13]
De[13]
Dims(ages)