# Bioturbation in Julia ####

using StatsBase
using BenchmarkTools
using CSV
using DataFrames
using HDF5, JLD2, FileIO
using ProgressBars
using Random
using Plots
using CSV
using KernelDensity, KernelDensitySJ

# Soil parameters
bd = 1500
nlayers = 200
layer_thickness = 0.01
bleaching_depth = 0.005
dz_tol = 0.55
soildepth = 2
grains_per_cm = 150

# index 1: thickness
# index 2: mass
# index 3: midpoint thickness
# index 4: cumulative thickness
soil = zeros(Float64, nlayers, 4)
soil2 = soil
ages = fill(Int[],nlayers)
De = fill(Float64[],nlayers)
DR = zeros(Float64,nlayers,1)
# Functions
initialise_soil = function(soildepth, layer_thickness, bleaching_depth, grains_per_cm)
    remainingsoil = soildepth
    for l in 1:nlayers
        if l==1 
            soil[l,1] = bleaching_depth
        else
            if remainingsoil > layer_thickness
                soil[l,1] = layer_thickness
            elseif remainingsoil <= 0
                soil[l,1] = 0
            else 
                soil[l,1] = remainingsoil
            end
        end
        remainingsoil -= soil[l,1]
    
        if l == nlayers && remainingsoil > 0
            soil[l,1] += remainingsoil
        end

        soil[l,2] = bd*soil[l,1] 
        ages[l] = fill(0, ceil(Int32, grains_per_cm * soil[l,1] * 100))
        De[l] = fill(0, ceil(Int32, grains_per_cm * soil[l,1] * 100))
        DR[l] = 10 - l*0.01
    end
end

transfer_OSL = function(layer, otherlayer, P_layer, P_otherlayer)
    ages_from = ages[layer]
    ages_to = ages[otherlayer]

    De_from = De[layer]
    De_to = De[otherlayer]

    ind_l = sample(0:1,ProbabilityWeights([1-P_layer,P_layer]),length(ages_from)).==1
    ind_ol = sample(0:1,ProbabilityWeights([1-P_otherlayer,P_otherlayer]),length(ages_to)).==1
    
    ages[layer] = vcat(ages_from[map(!,ind_l)], ages_to[ind_ol])
    ages[otherlayer] = vcat(ages_from[ind_l], ages_to[map(!,ind_ol)]) 

    De[layer] = vcat(De_from[map(!,ind_l)], De_to[ind_ol])
    De[otherlayer] = vcat(De_from[ind_l], De_to[map(!,ind_ol)]) 
end

recalculate_layer_thicknesses = function()
    for ll in 1:nlayers
        soil[ll,1] = soil[ll,2] / bd
    end
end

update_layers = function()
    for l in 1:(nlayers-1)

        if soil[l,2]/bd > (layer_thickness*(1+dz_tol))
            # Split
            if l < nlayers-2
                # merge bottom two layers and move layers down to create space
                soil[nlayers,2] += soil[nlayers-1,2]
                soil[nlayers-1,2] = 0
                transfer_OSL(nlayers-1, nlayers, 1,0)

                for ll in (nlayers-1):-1:(l+2)
                    soil[ll,2] = soil[ll-1,2]
                    soil[ll-1,2] = 0
                    transfer_OSL(ll-1, ll, 1,0)
                end
            end

            # Even splitting of the layers, except when resulting layer thickness exceeds initial layer thickness
            # In that case, the splitted layer becomes layer_thickness and the remainder is moved down for further splitting
            split_ratio = 0.5
            if (soil[l,2]/bd)/2 > layer_thickness
                split_ratio =  layer_thickness / (soil[l,2]/bd)
            end
            soil[l+1,2] = soil[l,2] * (1-split_ratio)
            soil[l,2] *= split_ratio
            transfer_OSL(l, l+1, (1-split_ratio), 0)
        end
        

        # Combine
        if l != 1 && soil[l,2]/bd < (layer_thickness*(1-dz_tol)) 
            # Determine thinnest neighbouring layer to merge with
            l2 = l+1
            if (l-1)>1 && soil[l-1,2] < soil[l+1,2]
                l2 = l-1
            end

            soil[l2,2] += soil[l,2]
            soil[l,2] = 0
            transfer_OSL(l, l2, 1,0)

            # Move layers up
            for ll in l:(nlayers-1)
                soil[ll,2] = soil[ll+1,2]
                soil[ll+1,2] = 0
                transfer_OSL(ll+1, ll, 1,0)
            end
        end
        recalculate_layer_thicknesses()
    end
end

update_OSL = function()
    # update thickness top layer to bleaching depth
    if soil[1,2]/bd > bleaching_depth # Layer too thick, give to layer below
        sep_fraction = (soil[1,2]/bd-bleaching_depth)/(soil[1,2]/bd)
        transport = sep_fraction * soil[1,2]
        soil[1,2] -= transport
        soil[2,2] += transport
        transfer_OSL(1,2,sep_fraction,0)
    end

    if soil[1,2]/bd < bleaching_depth # Layer too thin, take from layer below
        sep_fraction = (bleaching_depth-soil[1,2]/bd)/(soil[2,2]/bd)
        transport = sep_fraction * soil[2,2]
        soil[2,2] -= transport
        soil[1,2] += transport
        transfer_OSL(2,1,sep_fraction,0)
    end
    recalculate_layer_thicknesses()

    ages[1] = ages[1] .*0 # Bleach particles in the top layer
    De[1] = De[1] .* 0
    for l in 2:nlayers # Add a year to particles in all other layers
        ages[l] = ages[l] .+ 1 
        De[l] = De[l] .+ DR[l] 
    end
    
end

BT_layer_index = function(depth, depth_low, _depth_function, _corr_dd)
    tot_thickness = sum(soil[1:nlayers,1])

    if _depth_function == "grd" # Gradational / linear decline in rate
        dd = 1 * _corr_dd
        dd_ref = dd
        if depth < (1/dd)
            if depth_low > (1/dd)
                depth_low = 1/dd
            end
            layer_index = -dd/2 * (depth_low^2 - depth^2) + (depth_low - depth)
            total_index = -dd/2 * (1/dd)^2 + 1/dd
        else
            layer_index = 0
            total_index = 1
        end
    end
    
    if _depth_function == "exp" # Exponential decline in rate
        dd = 6 * _corr_dd
        dd_ref = dd
        layer_index = exp(-dd*depth) - exp(-dd * (depth_low))
        total_index = 1 - exp(-dd*tot_thickness)
    end

    if _depth_function == "abr" # Constant rate, or abrupt decline in rate
        dd = 1 * _corr_dd
        dd_ref = dd
        if depth < dd
            if depth_low > dd
                depth_low = dd
            end
            layer_index = depth_low - depth
            total_index = dd
        else
            layer_index = 0
            total_index = 1
        end
    end

    return layer_index / total_index
end

BT_mixing = function(_BT_pot, _depth_function)
    depth = 0

    for l in 1:nlayers
        # Calculate BT occurring in this layer
        BT_layer = _BT_pot * BT_layer_index(depth, depth + soil[l,1], _depth_function)
        if BT_layer > soil[l,2]
            BT_layer = soil[l,2]
        end

        depth += soil[l,1]/2

        if BT_layer > 0
            # Calculate exchange with other layers

            #exchange index
            tot_exch_index = -1/dd_exch * (exp(-dd_exch*(depth))-1)+ -1/dd_exch*(exp(-dd_exch*(soildepth-depth))-1)
            otherdepth = 0

            for ol in 1:nlayers
                z_upp = otherdepth
                z_low = otherdepth + soil[ol,1]
                if(ol<l) # above
                    ol_exch_index = -1/dd_exch*(exp(-dd_exch*(depth-z_upp)) - exp(-dd_exch * (depth - z_low)));
                end

                if(ol>l) # below
                    ol_exch_index = -1/dd_exch*(exp(-dd_exch*(z_low - depth)) - exp(-dd_exch * (z_upp - depth)));
                end

                otherdepth += soil[ol,1]

                if(l!=ol) # calculate exchange
                    interlayer_exchange = BT_layer * ol_exch_index / tot_exch_index
                    dmass_l = min(soil[l,2],interlayer_exchange / 2)
                    dmass_ol = min(soil[ol,2],interlayer_exchange / 2)

                    P_transfer_l = dmass_l / soil[l,2]
                    P_transfer_ol = dmass_ol / soil[ol,2]
                    transfer_OSL(l, ol, P_transfer_l, P_transfer_ol)
                end
            end
        end
        depth += soil[l,1]/2
    end
end

BT_mounding = function(_BT_pot, _depth_function, _corr_dd = 1)
    depth = 0
    for l in 2:nlayers
        # Calculate BT occurring in this layer
        BT_layer = _BT_pot * BT_layer_index(depth, depth + soil[l,1], _depth_function, _corr_dd)
        if BT_layer > soil[l,2]
            BT_layer = soil[l,2]
        end

        if BT_layer > 0
            # Transport to upper layer
            P_transfer_l = BT_layer / soil[l,2]
            soil[l,2] -= BT_layer
            soil[1,2] += BT_layer
            transfer_OSL(l, 1, P_transfer_l, 0)
        end
        depth += soil[l,1]
    end
end

BT_upheaval = function(_md, _freq, t)
    if t % _freq==0 & t != ntime
        particles_contributed = zeros(Int64, nlayers)
        mixed_particles = Array{Int64}(undef, 0)
        tot_particles = 0
        depth_mixed = 0

        # Take up particles
        for l in 1:nlayers
            if depth_mixed < _md
                if depth_mixed + soil[l,1] < _md
                    particles_contributed[l] = length(ages[l])
                    mixed_particles = vcat(mixed_particles, ages[l])
                    ages[l] = []
                else
                    particles_contributed[l] = round(Int32,length(ages[1]) * (_md - depth_mixed) / soil[l,1])
                    ages[l] = shuffle(ages[l])
                    mixed_particles = vcat(mixed_particles, ages[l][1:particles_contributed[l]])
                    ages[l] = ages[l][(particles_contributed[l]+1):length(ages[l])]
                end
                depth_mixed += soil[l,1]
            end
            if depth_mixed >= _md
                break
            end
        end
        
        # Mix
        mixed_particles = shuffle(mixed_particles)
        tot_particles = sum(particles_contributed)

        # And give back to layers
        count = 1
        for l in 1:nlayers
            ages[l] = vcat(ages[l],mixed_particles[count:(count+particles_contributed[l]-1)])
            count += particles_contributed[l]
            if count >= tot_particles
                break
            end
        end
    end
end

KDE = function(data)
    avg = mean(data)
    stdev = std(data)
    IQR = iqr(data)
    sigma = 0.9*min(stdev,IQR/1.34)*length(data)^(-0.2)
    
    MAX = -99999999999
    MIN = 9999999999
    
    N = length(data)
    
    for i in 1:N
        if MAX < data[i]
            MAX = data[i]
        end
        if MIN > data[i]
            MIN = data[i]
        end
    end
    
    MIN = floor(Int32,MIN)
    MAX = ceil(Int32, MAX)
    
    #=nsteps = MAX - MIN
    if nsteps<500
        nsteps = 500
    end
    while nsteps > 512
        nsteps /= 2
    end
    nsteps = round(Int32, nsteps)
    =#
    nsteps = 2000 # every 5 years
    result = zeros(Float64, nsteps,2)
    x = zeros(nsteps,1)
    y = zeros(nsteps,1)
    
    
    # x[1] = MIN
    x[1] = 0
    for i in 2:nsteps
        # x[i] = x[i - 1] + ((MAX - MIN) / nsteps)
        x[i] = x[i-1]+ 5
    end
    
    # kernel density estimation
    c = 1.0 / (sqrt(2 * pi * sigma * sigma))
    
    for i in 1:N
        for j in 1:nsteps
            y[j] = y[j] + 1.0 / N * c * exp(-(data[i] - x[j]) * (data[i] - x[j]) / (2 * sigma * sigma))
        end
    end
    
    for i in 1:nsteps
        result[i, 1] = x[i]
        result[i, 2] = y[i]
    end
    return result
end

agreement_distributions = function(ages_field, ages_model) 
    # Silverman's approach, own function
    #kde1=KDE(ages_field)
    #kde2=KDE(ages_model)
    #kde1=hcat(kde1,kde2[:,2])
    #kde1 = hcat(kde1,zeros(2000,1))

    # Sheather and Jones bandwidth estimation
    bw1=bwsj(ages_field)
    bw2=bwsj(ages_model)

    kde1=kde(ages_field,bandwidth=bw1,boundary=(0,10000),npoints=2001)
    kde2=kde(ages_model,bandwidth=bw1,boundary=(0,10000),npoints=2001)

    kde1=hcat(kde1.x,kde1.density,kde2.density,zeros(length(kde1.x)))

    for i in 1:size(kde1)[1]
        kde1[i,4] = minimum(kde1[i,2:3])
    end

    agreement = sum(sum(kde1[:,4])/sum(kde1[:,2]) + sum(kde1[:,4])/sum(kde1[:,3]))/2
    #plot(kde1[:,1], kde1[:,2])
    if isnan(agreement) | isinf(agreement)
    agreement = 0
    end
    return(agreement)
end


calib_shift_and_zoom = function(para_number, zoom_factor, orig_value)
    mid_ratio = 1

end


                double mid_ratio = 0;
                if (calib_ratios.GetLength(1) % 2 == 0) { mid_ratio = (calib_ratios[para_number, Convert.ToInt32(calib_ratios.GetLength(1) / 2) - 1] + calib_ratios[para_number, Convert.ToInt32((calib_ratios.GetLength(1) / 2))]) / 2; }
                else { mid_ratio = calib_ratios[para_number, Convert.ToInt32(calib_ratios.GetLength(1) / 2 - 0.5)]; }
                Debug.WriteLine("mid ratio is " + mid_ratio);
                Double best_ratio = best_parameters[para_number] / orig_par_value;
                Debug.WriteLine("best ratio is " + best_ratio);
                if (best_parameters[para_number] == calib_ratios[para_number, 0] * orig_par_value | best_parameters[para_number] == calib_ratios[para_number, calib_ratios.GetLength(1) - 1] * orig_par_value)
                {
                    //the best parameter ratio (and thus value) was on the edge of the range. We must shift our range sideways (we keep the same ratio between upper and lower ratio - are you still with me?)
                    Debug.WriteLine(" currentpara value was on edge of range");

                    for (int ratio = 0; ratio < calib_ratios.GetLength(1); ratio++)
                    {
                        Debug.WriteLine(" setting ratio " + calib_ratios[para_number, ratio] + " to " + calib_ratios[para_number, ratio] * (best_ratio / mid_ratio));
                        calib_ratios[para_number, ratio] = calib_ratios[para_number, ratio] * (best_ratio / mid_ratio);
                    }
                }
                else
                {
                    //the best parameter ratio (and thus value) NOT on the edge of the range. We must shift to the best observed value and then zoom IN
                    Debug.WriteLine(" currentpara value was NOT on edge of range");
                    for (int ratio = 0; ratio < calib_ratios.GetLength(1); ratio++)
                    {
                        Debug.Write(" setting ratio " + calib_ratios[para_number, ratio] + " to " );
                        calib_ratios[para_number, ratio] = best_ratio + (((calib_ratios[para_number, ratio]/mid_ratio)*best_ratio)-best_ratio)/ zoom_factor;
                        Debug.WriteLine(calib_ratios[para_number, ratio]);
                    }
                }
            }
            catch { Debug.WriteLine(" problem adapting parameters and ratios "); }
        }

write_CSV_output = function(directory, soils, ages)
    _soil = DataFrame(soil, ["Thickness", "Mass", "cumThickness"])
    _ages = zeros(Int32,sum(length, ages), 2)
    index = 1
    for l in 1:nlayers
        n = length(ages[l])
        _ages[index:index+(n-1),1] .= l
        _ages[index:index+(n-1),2] = ages[l]
        index += n
    end
    _ages = DataFrame(_ages,["Layer", "Age"])
    CSV.write(string(directory,"soils.csv"), _soil)
    CSV.write(string(directory,"ages.csv"), _ages)
end




# Simulations

# Exp 1. Testing different depth functions for different bioturbation types
type = ["mounding", "mixing"]
depth_function = ["grd", "exp", "abr"]
ntime = 10000
dd_exch = 10
BT_pot = 10

# Create final output arrays
output_soil = zeros(length(type), length(depth_function), nlayers,4)
output_ages = fill(Int[],length(type), length(depth_function), nlayers)
runtimes = zeros(length(type)*length(depth_function))

# Run simulations
ind = 0
for i1 in 1:length(type)
    for i2 in 1:length(depth_function)
        t1 = time()
        ind += 1
        println("Run $ind of $(length(type)*length(depth_function))")
        println(string("Process: $(type[i1]). Depth function: $(depth_function[i2])"))
        initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_cm)

        for t in 1:ntime
            if t%100 == 0
                println("t ",t)
            end
            if type[i1] == "mixing"
                BT_mixing(BT_pot,depth_function[i2])
            end
            if type[i1] == "mounding"
                BT_mounding(BT_pot, depth_function[i2])
            end
        
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
        output_soil[i1,i2,:,:] = soil
        output_ages[i1,i2,:] = ages

        runtimes[ind] = time() - t1
    end
end

runtimes_min = runtimes / 60
runtimes_min[6]
wd = "C:\\Simulations\\Bioturbation\\Julia\\"
save(string(wd,"exp1_output_soil.jld2"), "output_soil", output_soil)
save(string(wd,"exp1_output_ages.jld2"), "output_ages", output_ages)


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
    initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_cm)

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

        initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_cm)

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

    initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_cm)


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
    initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_cm)

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
termites = CSV.read("C:\\Users\\Marijn\\sciebo\\Research\\Bioturbation\\Data\\Termites_Ghana\\De_termites_Ghana_complete.csv", DataFrame)
depths = unique(termites[:,"depth"])
ntime = 4000
BT_pot = 10
dd_ref = 1
cal_depthfunction = ["grd", "exp"]

# Fractions to multiply input to
calib_ratios=hcat([.25,.5,1,2,4],[.25,.5,1,2,4])
calib_results = hcat(0, 0, -1)



for(it in 1:3)
    depth_function = "grd"
    for i2 in ProgressBar(1:length(cal_BT_pot))
        for i3 in ProgressBar(1:length(cal_dd))
            initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_cm)

            # Simulations
            for t in ProgressBar(1:ntime)
                BT_mounding(BT_pot * cal_BT_pot[i2], depth_function, cal_dd[i3])
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

            calib_results = vcat(BT_pot_* cal_BT_pot[i2], dd_ref, total_agreement)
            # cal_results[i1,i2,i3] = total_agreement
        end
    end
    zoom_and_shift_ratios
end




ind=findall(x->x==maximum(cal_results), cal_results)
println("Best result: df: $(cal_depthfunction[ind[1]]), BTpot_corr: $(cal_BT_pot[ind[2]]), dd_corr: $(cal_dd[ind[3]])")

plot()




# Exp XX. delta dose rate
initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_cm)
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


cols = :gray()


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