# Packages
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

# Functions
initialise_soil = function(soildepth, layer_thickness, bleaching_depth, grains_per_layer) # Fill the soil and age arrays based on the parameters
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
        ages[l] = fill(0, ceil(Int32, grains_per_layer))
    end
end

transfer_OSL = function(layer, otherlayer, P_layer, P_otherlayer) # Transfer OSL particles in between layers, based on transfer probability
    ages_from = ages[layer]
    ages_to = ages[otherlayer]

    ind_l = sample(0:1,ProbabilityWeights([1-P_layer,P_layer]),length(ages_from)).==1
    ind_ol = sample(0:1,ProbabilityWeights([1-P_otherlayer,P_otherlayer]),length(ages_to)).==1
    
    ages[layer] = vcat(ages_from[map(!,ind_l)], ages_to[ind_ol])
    ages[otherlayer] = vcat(ages_from[ind_l], ages_to[map(!,ind_ol)]) 
end

recalculate_layer_thicknesses = function() # Update layer thicknesses after transport of material
    for ll in 1:nlayers
        soil[ll,1] = soil[ll,2] / bd
    end
end

update_layers = function() # Split or combine layers based on their thickness and the thickness tolerance (55%)
    for l in 1:(nlayers-1)

        if soil[l,2]/bd > (layer_thickness*(1+0.55))
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
        if l != 1 && soil[l,2]/bd < (layer_thickness*(1-0.55)) 
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

update_OSL = function() # Update OSL properties: thickness of bleaching layer and OSL ages
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
    for l in 2:nlayers # Add a year to particles in all other layers
        ages[l] = ages[l] .+ 1 
    end
end

BT_layer_index = function(depth, depth_low, _depth_function, _dd) # Determine fraction of bioturbation occuring in a layer, depending on the depth function
    tot_thickness = sum(soil[1:nlayers,1])

    if _depth_function == "grd" # Gradational / linear decline in rate
        if depth < (1/_dd)
            if depth_low > (1/_dd)
                depth_low = 1/_dd
            end
            layer_index = -_dd/2 * (depth_low^2 - depth^2) + (depth_low - depth)
            total_index = -_dd/2 * (1/_dd)^2 + 1/_dd
        else
            layer_index = 0
            total_index = 1
        end
    end
    
    if _depth_function == "exp" # Exponential decline in rate
        layer_index = exp(-_dd*depth) - exp(-_dd * (depth_low))
        total_index = 1 - exp(-_dd*tot_thickness)
    end

    if _depth_function == "abr" # Constant rate, or abrupt decline in rate
        if depth < _dd
            if depth_low > _dd
                depth_low = _dd
            end
            layer_index = depth_low - depth
            total_index = _dd
        else
            layer_index = 0
            total_index = 1
        end
    end

    return layer_index / total_index
end

BT_mixing = function(_BT_pot, _depth_function) # Bioturbation by subsurface mixing
    depth = 0

    for l in 1:nlayers
        # Calculate BT occurring in this layer
        BT_layer = _BT_pot * BT_layer_index(depth, depth + soil[l,1], _depth_function, _dd)
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

BT_mounding = function(_BT_pot, _depth_function) # Bioturbation by mounding
    depth = 0
    for l in 2:nlayers
        # Calculate BT occurring in this layer
        BT_layer = _BT_pot * BT_layer_index(depth, depth + soil[l,1], _depth_function, _dd)
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

BT_upheaval = function(md, freq, t) # Bioturbation by upheaval
    if t % freq==0 & t != ntime # determine is there is an upheaval event in this timestep

        # Pepare for mixing
        particles_contributed = zeros(Int64, nlayers)
        mixed_particles = Array{Int64}(undef, 0)
        tot_particles = 0
        depth_mixed = 0

        # Take up particles from source layers
        for l in 1:nlayers
            if depth_mixed < md
                if depth_mixed + soil[l,1] < md
                    particles_contributed[l] = length(ages[l])
                    mixed_particles = vcat(mixed_particles, ages[l])
                    ages[l] = []
                else
                    particles_contributed[l] = round(Int32,length(ages[1]) * (md - depth_mixed) / soil[l,1])
                    ages[l] = shuffle(ages[l])
                    mixed_particles = vcat(mixed_particles, ages[l][1:particles_contributed[l]])
                    ages[l] = ages[l][(particles_contributed[l]+1):length(ages[l])]
                end
                depth_mixed += soil[l,1]
            end
            if depth_mixed >= md
                break
            end
        end
        
        # Mix OSL particles
        mixed_particles = shuffle(mixed_particles)

        # And give back to layers
        tot_particles = sum(particles_contributed)

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

plot_age_depth_profile = function(ages, soil, aggr_layers, every_nth_layer, xlim = [0,ntime], ylim = [0,1.5])
    plt = plot(0,  xlims = xlim, ylims = ylim, xlabel="Age [a]", ylabel="Depth [m]",plot_title="Bioturbated fraction [-]", plot_titlefontsize = 12); yflip!(true)
    modes = Array{Float64}(undef, 0)
    depths = Array{Float64}(undef, 0)
    fbios = Array{Float64}(undef, 0)
    for l in 1:round(Int32,nlayers/aggr_layers)
        lay_sel = range((l-1)*aggr_layers+1, (l)*aggr_layers)
        data = vcat(copy(ages[lay_sel])...)
        l1=length(data)
        filter!(e->e≠ntime,data)
        fbio = (length(data)/l1)*10000
        depth = mean(soil[lay_sel,3])

        bw = bwsj(data)
        mode = NaN
        if(!isnan(bw))
            kde_ages = KernelDensity.kde(data,bandwidth=bw,boundary=(0,ntime*1.1),npoints=round(Int64,ntime*1.1/5+1))
            mode = kde_ages.x[kde_ages.density.==maximum(kde_ages.density)]
            kde_ages.density = kde_ages.density./maximum(kde_ages.density)
            kde_ages.density[kde_ages.density.<0.01] .= NaN # remove entries with P<0.05. Can also be set to other values
            kde_ages.density = kde_ages.density.*-0.1 .+ depth

            lab_plot = "Age distributions"
            if l>1
                lab_plot = false
            end

            if (l-1) %every_nth_layer == 0
                plot!(plt, kde_ages.x,kde_ages.density, c=:black, label = lab_plot)    
            end
        end
        modes=vcat(modes,mode)
        depths = vcat(depths,depth)
        fbios = vcat(fbios,fbio)
    end
    plot!(plt, modes, depths,label = "Modes",lw=2,c=:black)
    plot!(plt, fbios, depths, label = "Bioturbated fraction", c=:red, legend=:bottomleft)
    copy_ticks(current(),([0,10000],["0","1"]))

    return(plt)
end

    

# Code to add extra axis on top, sourced from https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
function twiny(sp::Plots.Subplot)
    sp[:top_margin] = max(sp[:top_margin], 30Plots.px)
    plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
    twinsp = sp.plt.subplots[end]
    twinsp[:xaxis][:mirror] = true
    twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
    Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
    twinsp
end

# # twiny(plt::Plots.Plot = current()) = twiny(plt[1])

function copy_ticks(sp::Plots.Subplot,xticks)
    ptx = twinx(sp)
    plot!(ptx,xlims=xlims(plt),ylims=ylims(plt),yformatter=_->"",xformatter=_->"", grid=false)
    pty = twiny(sp)
    plot!(pty,xlims=xlims(plt),ylims=ylims(plt),yformatter=_->"",xticks=xticks, grid=false)
end
# copy_ticks(plt::Plots.Plot,xticks) = copy_ticks(plt[1],xticks)



#=KDE = function(data) # Function for calculating kernel density functions
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
=#

#=
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
    kde2=kde(ages_model,bandwidth=bw2,boundary=(0,10000),npoints=2001)

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
=#