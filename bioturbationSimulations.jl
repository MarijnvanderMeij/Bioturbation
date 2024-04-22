#### Bioturbation in Julia ####
include("bioturbationFunctions.jl")


## Objects and parameters 
bd = 1500               # Bulk density [kg m-3]
nlayers = 70           # Number of soil layers
layer_thickness = 0.02  # Initial layer thickness [m]
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
ntime = 10000
initialise_soil(soildepth, layer_thickness, bleaching_depth, grains_per_layer)
t1 = time()
for t in ProgressBar(1:ntime)
    
    BT_mounding(BT_pot, "grd", 1)

    update_layers()
    update_OSL()
end
t2 = time()
t2-t1

# A run with bioturbation by subsurface mixing
ntime = 5000
dd = 4.5
dd_exch = 5
soil = zeros(Float64, nlayers, 4)   # Array of the soil profile. Index 1: thickness [m]. Index 2: mass [kg]. Index 3: midpoint depth [m]. Index 4: cumulative depth [m].
ages = fill(Int[],nlayers)          # Array of arrays of OSL age particles
soil, ages = initialise_soil(soil, ages, soildepth, layer_thickness, bleaching_depth, grains_per_layer)

t1 = time()
for t in ProgressBar(1:ntime)
    
    soil, ages = BT_mixing(soil, ages, BT_pot, "exp", 6)

    soil, ages = update_layers(soil, ages)
    soil, ages = update_OSL(soil, ages)
end
depth = 0
for l in 1:nlayers
    depth += soil[l,1]/2
    soil[l,3] = depth
    depth += soil[l,1]/2
end
soil[:,4]=cumsum(soil[:,1])
t2 = time()
t2-t1
plot_age_depth_profile(ages, soil, 5, 2)

# 1. dd - 6
# 2. dd = 4.5, dd_exch = 10
# 3. dd = 4.5, dd_exch = 5
s_dd45 = soil
a_dd45 = ages



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
termites = CSV.read("C:\\Users\\wmeij\\Simulations\\bioturbation\\De_termites_complete.csv", DataFrame)
termites[:,"Ages"] = termites[:,"De"]./termites[:,"DR"]*1000  

depths = unique(termites[:,"depth"])


ntime = 4000
BT_pot = 1
dd_ref = 1.02
dd_exch = 10
grains_per_layer = 300

filter!(e->e <= dd_ref,depths)

# Fractions to multiply input to
cal_depthfunction = ["grd", "abr", "exp"]
cal_BT_pot = 3:.25:6.0 # [0.5, 1, 1.5]#[.1,.25,.5,1,2,4]
cal_dd = 1 # cal_dd = .5:.1:1# [.1,.25,.5,1,2,4,10]
cal_rel_process = 0:.1:1

calib_results_termites = zeros(Float64, length(cal_depthfunction)*length(cal_BT_pot)*length(cal_dd)*length(cal_rel_process), 7) 
ct = 0
for i1 in 1:length(cal_depthfunction)
    for i2 in 1:length(cal_BT_pot)
        for i3 in 1:length(cal_dd)
            for i4 in 1:length(cal_rel_process)

                ct += 1
                calib_results_termites[ct,1] = i1
                calib_results_termites[ct,2] = cal_BT_pot[i2] * BT_pot

                if cal_depthfunction[i1]=="exp"
                    calib_results_termites[ct,3] = -log(0.0025) / (dd_ref * cal_dd[i3])
                end

                if cal_depthfunction[i1]=="grd"
                    calib_results_termites[ct,3] = cal_dd[i3] * 1 / dd_ref
                end

                if cal_depthfunction[i1]=="abr"
                    calib_results_termites[ct,3] = cal_dd[i3] * dd_ref
                end
                calib_results_termites[ct,4] = cal_rel_process[i4]
            end
        end
    end
end
size(calib_results_termites)

bd = 1500               # Bulk density [kg m-3]
nlayers = 75           # Number of soil layers
layer_thickness = 0.02  # Initial layer thickness [m]
bleaching_depth = 0.005 # Thickness of surface layer where bleaching occurs [m]
soildepth = 1.5           # Thickness of the entire soil profile [m]
grains_per_layer = 300  # Initial number of layers per soil layer



ntime = 4000
@time Threads.@threads for i in 1:size(calib_results_termites)[1]   
    s, a = BT_run(calib_results_termites, i)
    calib_results_termites[i,5:7] .= BT_calibration(s, a, termites, true, true, true)
    if i%100==0
        GC.gc()
    end
end


CSV.write("C:\\Users\\wmeij\\Simulations\\bioturbation\\calibration_termites_squaredError.csv", Tables.table(calib_results_termites), writeheader=false)

calib_results_termites = CSV.read("C:\\Users\\wmeij\\Simulations\\bioturbation\\calibration_termites_squaredError.csv", DataFrame, header=false)
calib_results_termites = Matrix(calib_results_termites)
x = calib_results_termites[:,2]
y = vec(sum(calib_results_termites[:,5:7], dims = 2))

# plot(x, y, xlabel = "BT_pot", ylabel = "calibration error", legend= false)
# scatter!(x, y)
# calib_results_termites = CSV.read("C:\\Users\\wmeij\\Simulations\\bioturbation\\calibration_termites_squaredError.csv", DataFrame, header = false)
calib_results_termites_df = DataFrame(; depthfunction = calib_results_termites[:,1],
BT_pot = calib_results_termites[:,2],
cal_dd = calib_results_termites[:,3],
cal_rel_process = calib_results_termites[:,4],
err_mode = calib_results_termites[:,5],
err_fbio = calib_results_termites[:,6],
err_iqr = calib_results_termites[:,7],
# err_total = vec(sum(calib_results_termites[:,5:7], dims = 2) ./ 3)
err_total = vec(sum.(eachrow(calib_results_termites[:, 5:7])) ./ 3)
)

using RCall
R"""
library(ggplot2)
library(dplyr)
"""

R"""
rdf <- $calib_results_termites_df
cal_termites_gg = ggplot(rdf, aes(x = cal_rel_process, y = err_total, color = as.factor(BT_pot), linetype = as.factor(depthfunction))) + 
geom_line() + 
labs(x = "Percentage mounding", y = "Calibration error", color = "BTpot [kg/m2/a]", linetype = "depth function") + 
scale_linetype_manual(labels = c("Gradational", "Abrupt", "Exponential"), values = c(1,2, 3))
cal_termites_gg
# ggsave("Figures\\calibration_curve_termites.png", cal_termites_gg)
"""

# run with best parameter
y=vec(sum.(eachrow(calib_results_termites[:, 5:7])) ./ 3)
ind_termites = findall(y .== minimum(y))[1]
calib_results_termites[ind_termites,:]
s, a = BT_run(calib_results_termites, ind_termites)

fig = plot_age_depth_profile(a,s,1, 1, true)
fig = add_calibration_PDFs(termites, true)
savefig(fig,"Figures\\termites_calibrated.pdf")


# calibrate worms dataset
worms = CSV.read("C:\\Users\\wmeij\\Simulations\\bioturbation\\De_worms_complete.csv", DataFrame)
worms[:,"Ages"] = worms[:,"De"]./worms[:,"DR"]*1000  .- 3800

depths = unique(worms[:,"depth"])

ntime = 13200-3800
BT_pot = 1 # kg m-2
dd_ref = 0.6
dd_exch = 10

filter!(e->e <= dd_ref,depths)

# Fractions to multiply input to
cal_depthfunction = ["grd", "exp", "abr"]
cal_BT_pot = .25:.25:3 # [0.5, 1, 1.5]#[.1,.25,.5,1,2,4]
cal_dd = 1 # cal_dd = .5:.1:1# [.1,.25,.5,1,2,4,10]
cal_rel_process = 0:.1:1

calib_results_worms = zeros(Float64, length(cal_depthfunction)*length(cal_BT_pot)*length(cal_dd)*length(cal_rel_process), 7) 
ct = 0
for i1 in 1:length(cal_depthfunction)
    for i2 in 1:length(cal_BT_pot)
        for i3 in 1:length(cal_dd)
            for i4 in 1:length(cal_rel_process)

                ct += 1
                calib_results_worms[ct,1] = i1
                calib_results_worms[ct,2] = cal_BT_pot[i2] * BT_pot

                if cal_depthfunction[i1]=="exp"
                    calib_results_worms[ct,3] = -log(0.0025) / (dd_ref * cal_dd[i3])
                end

                if cal_depthfunction[i1]=="grd"
                    calib_results_worms[ct,3] = cal_dd[i3] * 1 / dd_ref
                end

                if cal_depthfunction[i1]=="abr"
                    calib_results_worms[ct,3] = cal_dd[i3] * dd_ref
                end
                calib_results_worms[ct,4] = cal_rel_process[i4]
            end
        end
    end
end

bd = 1500               # Bulk density [kg m-3]
nlayers = 50            # Number of soil layers
layer_thickness = 0.02  # Initial layer thickness [m]
bleaching_depth = 0.005 # Thickness of surface layer where bleaching occurs [m]
soildepth = 1           # Thickness of the entire soil profile [m]
grains_per_layer = 300  # Initial number of layers per soil layer
# BT_pot = 10             # Potential bioturbation rate [kg m-2]
# dd = 1                  # Depth parameter [m-1]
# dd_exch = 10            # Depth parameter for subsurface mixing [m-1]

soil_error = 0
age_error = 0

@time Threads.@threads for i in 1:size(calib_results_worms)[1]    
    s, a = BT_run(calib_results_worms, i)
    try
        calib_results_worms[i,5:7] .= BT_calibration(s, a, worms, true, true, true)   
    catch
        global soil_error
        global age_error
        soil_error = s
        age_error = a
        println("error in iteration ", i)
        calib_results_worms[i,5:7] .= NaN
    end
    if(i%100==0)
        GC.gc()
    end
end


CSV.write("C:\\Users\\wmeij\\Simulations\\bioturbation\\calibration_worms_squaredError.csv", Tables.table(calib_results_worms), writeheader=false)

# calib_results_worms = CSV.read("C:\\Users\\wmeij\\Simulations\\bioturbation\\calibration_worms_squaredError.csv", DataFrame, header=false)
# calib_results_worms = Matrix(calib_results_worms)
x = calib_results_worms[:,4]
y = vec(sum(calib_results_worms[:,5:7], dims = 2))
# lot(x, y, xlabel = "rel mixing", ylabel = "calibration error", legend= false)
# scatter!(x, y)

calib_results_worms_df = DataFrame(; depthfunction = calib_results_worms[:,1],
BT_pot = calib_results_worms[:,2],
cal_dd = calib_results_worms[:,3],
cal_rel_process = calib_results_worms[:,4],
err_mode = calib_results_worms[:,5],
err_fbio = calib_results_worms[:,6],
err_iqr = calib_results_worms[:,7],
err_total = vec(sum(calib_results_worms[:,5:7], dims = 2) ./ 3)
)
calib_results_worms_df_mat = Matrix{Float64}, (calib_results_worms_df)
sum.(eachrow(calib_results_worms[:,5:7]))

using RCall
R"""
library(ggplot2)
library(dplyr)
"""

R"""
rdf <- $calib_results_worms_df
cal_termites_gg = ggplot(rdf, aes(x = cal_rel_process, y = err_total, color = as.factor(BT_pot), linetype = as.factor(depthfunction))) + 
geom_line() + 
labs(x = "Percentage mounding", y = "Calibration error", color = "BTpot [kg/m2/a]", linetype = "depth function") + 
scale_linetype_manual(labels = c("Gradational", "Abrupt", "Exponential"), values = c(1,2, 3))
cal_termites_gg
# ggsave("Figures\\calibration_curve_worms.png", cal_termites_gg)
"""

# calib_results_worms[156,:]
# run with best parameter
ind_worms = findall(y .== minimum(y))[1]
# BT_pot_cal = x[]
calib_results_worms[ind_worms,:]

s, a = BT_run(calib_results_worms, ind_worms)

fig = plot_age_depth_profile(a,s,1, 1, true)
fig = add_calibration_PDFs(worms, true)
savefig(fig,"Figures\\worms_calibrated.pdf")


## Calibrate ants dataset
ants = CSV.read("C:\\Users\\wmeij\\Simulations\\bioturbation\\De_ants_complete.csv", DataFrame)
ants[:,"Ages"] = ants[:,"De"]./ants[:,"DR"]*1000 

depths = unique(ants[:,"depth"])

ntime = 50000
BT_pot = 1 # kg m-2
dd_ref = 0.5
dd_exch = 10

filter!(e->e <= dd_ref,depths)

# Fractions to multiply input to
cal_depthfunction = ["grd", "exp", "abr"]
cal_BT_pot = .25:.25:4 # [0.5, 1, 1.5]#[.1,.25,.5,1,2,4]
cal_dd = 1 # cal_dd = .5:.1:1# [.1,.25,.5,1,2,4,10]
cal_rel_process = 0:.1:1

calib_results_ants = zeros(Float64, length(cal_depthfunction)*length(cal_BT_pot)*length(cal_dd)*length(cal_rel_process), 7) 
ct = 0
for i1 in 1:length(cal_depthfunction)
    for i2 in 1:length(cal_BT_pot)
        for i3 in 1:length(cal_dd)
            for i4 in 1:length(cal_rel_process)

                ct += 1
                calib_results_ants[ct,1] = i1
                calib_results_ants[ct,2] = cal_BT_pot[i2] * BT_pot

                if cal_depthfunction[i1]=="exp"
                    calib_results_ants[ct,3] = -log(0.0025) / (dd_ref * cal_dd[i3])
                end

                if cal_depthfunction[i1]=="grd"
                    calib_results_ants[ct,3] = cal_dd[i3] * 1 / dd_ref
                end

                if cal_depthfunction[i1]=="abr"
                    calib_results_ants[ct,3] = cal_dd[i3] * dd_ref
                end
                calib_results_ants[ct,4] = cal_rel_process[i4]
            end
        end
    end
end
size(calib_results_ants)


bd = 1500               # Bulk density [kg m-3]
nlayers = 11            # Number of soil layers
layer_thickness = 0.05  # Initial layer thickness [m]
bleaching_depth = 0.005 # Thickness of surface layer where bleaching occurs [m]
soildepth =   0.55       # Thickness of the entire soil profile [m]
grains_per_layer = 300  # Initial number of layers per soil layer


soil_error = 0
age_error = 0
ntime = 50000
@time Threads.@threads for i in 1:size(calib_results_ants)[1]    
    s, a = BT_run(calib_results_ants, i)
    try
        calib_results_ants[i,5:7] .= BT_calibration(s, a, ants, true, true, true)   
    catch
        global soil_error
        global age_error
        soil_error = s
        age_error = a
        println("error in iteration ", i)
        calib_results_ants[i,5:7] .= NaN
    end
    if(i%100==0)
        GC.gc()
    end
end

CSV.write("C:\\Users\\wmeij\\Simulations\\bioturbation\\calibration_ants_squaredError.csv", Tables.table(calib_results_ants), writeheader=false)

# calib_results_ants = CSV.read("C:\\Users\\wmeij\\Simulations\\bioturbation\\calibration_ants_squaredError.csv", DataFrame, header=false)
# calib_results_ants = Matrix(calib_results_ants)
x = calib_results_ants[:,4]
y = vec(sum(calib_results_ants[:,5:7], dims = 2))
# plot(x, y, xlabel = "rel mixing", ylabel = "calibration error", legend= false)
# scatter!(x, y)

calib_results_ants_df = DataFrame(; depthfunction = calib_results_ants[:,1],
BT_pot = calib_results_ants[:,2],
cal_dd = calib_results_ants[:,3],
cal_rel_process = calib_results_ants[:,4],
err_mode = calib_results_ants[:,5],
err_fbio = calib_results_ants[:,6],
err_iqr = calib_results_ants[:,7],
err_total = vec(sum(calib_results_ants[:,5:7], dims = 2) ./ 3)
)
# calib_results_ants_df_mat = Matrix{Float64}, (calib_results_ants)
# sum.(eachrow(calib_results_ants[:,5:7]))


using RCall
R"""
library(ggplot2)
library(dplyr)
"""

R"""
rdf <- $calib_results_ants_df
cal_ants_gg = ggplot(rdf, aes(x = cal_rel_process, y = err_total, color = as.factor(BT_pot), linetype=as.factor(depthfunction))) + 
geom_line() + 
labs(x = "Percentage subsurface mixing", y = "Calibration error", color = "BTpot [kg/m2/a]")
cal_ants_gg
# ggsave("Figures\\calibration_curve_ants.png", cal_ants_gg)
"""


# run with best parameter
ind_ants = findall(y .== minimum(y))[1]
calib_results_ants[ind_ants,:]
# BT_pot_cal = x[]

s, a = BT_run(calib_results_ants, ind_ants)

fig = plot_age_depth_profile(a,s,2, 1, true)
fig = add_calibration_PDFs(ants, true)
savefig(fig,"Figures\\ants_calibrated.pdf")
















N = 15
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]
cl = hcat(calib_results, vcat(repeat([NaN], outer = Int32(floor(N/2))), moving_average(calib_results[:,4], N), repeat([NaN], outer = Int32(floor(N/2)))))
plot(cl[:,2],cl[:,4:5])

cols=["red","black","blue"]

for df in [1.0,2.0,3.0]
    for dd in cal_dd
        temp = calib_results[(calib_results[:,1].==df) .& (calib_results[:,3].==dd),:]
        x = temp[:,2]
        y = temp[:,4]
        p=plot!(x,y,title = dd, label = string("df: ", df),color=cols[Int32(df)])
    end
end
current()

subset(calib_results(), :       )

countmap(calib_results[:,1].==df .* calib_results[:,2].==dd)
ind=findall(x->x==maximum(filter(!isnan,cl[:,4])), cl[:,4])
cl[ind,:]
ind = findall(x->x==0,calib_results[:,4])
ind=findall(calib_results[:,4]->calib_results[:,4]==0)
cl[ind,:]


soil, ages = BT_run(ind[1])
BT_calibration(soil, ages, termites, "Cumulative Distribution")

plot_age_depth_profile(ages, soil, 5, 3)
add_calibration_PDFs(termites)
current()

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