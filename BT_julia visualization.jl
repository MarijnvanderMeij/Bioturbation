# Plotting
using Plots, Colors, ColorSchemes
using HDF5, JLD2, FileIO
using StatsBase, DataStructures
using KernelDensity, KernelDensitySJ
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


get_modes = function(data_ages, data_soil, aggr_layers)
    modes = zeros(nlayers,3)
    fill!(modes, NaN)

    for l in 1:round(Int32,nlayers/aggr_layers)
        lay_sel = range((l-1)*aggr_layers+1, (l)*aggr_layers)

        data = vcat(copy(data_ages[lay_sel])...)
        l1=length(data)
        filter!(e->e≠10000,data)
        nsf = (1- length(data)/l1)*10000
        depth = mean(data_soil[lay_sel,3])
        if length(data) > 2
            bw = bwsj(data)
            dens = kde(data, bwsj)
            mode = dens[argmax(dens[:,2]),1]
        else
            mode = NaN
        end
        modes[l,:] = [mode, depth, nsf]       
    end
    return(modes)
end



# Exp 1. Comparison processes and depth functions
type = ["mounding", "mixing"]
depth_function = ["grd", "exp", "abr"]
wd = "C:\\Simulations\\Bioturbation\\Julia\\"
output_soil = load(string(wd,"exp1_output_soil.jld2"))["output_soil"]
output_ages = load(string(wd,"exp1_output_ages.jld2"))["output_ages"]


for i in 1:2
    for j in 1:3
        write_CSV_output(string("C:\\Simulations\\Bioturbation\\Julia\\exp1_",type[i],"_",depth_function[j],"_"),
        output_soil[i,j,:,:], output_ages[i,j,:,:])
    end
end



# write_CSV_output("C:\\Simulations\\Bioturbation\\Julia\\exp1_",output_soil, output_ages)
nlayers = size(output_soil)[3]

vars = 1:10
cols = distinguishable_colors(length(depth_function)+1, [RGB(1,1,1)])[2:end]
pcols = map(col -> (red(col), green(col), blue(col)), cols)

aggr_layers = 5
#using PlotlyJS

modes_plotting = zeros(6,200,3)
ind = 0
for i1 in 1:length(type)
    plot(0, xlims=[0,10000], ylims = [0,1.5], legend=:bottomleft, show=true); yflip!(true)

    for i2 in 1:length(depth_function)
        ind += 1
        mds = get_modes(output_ages[i1,i2,:], output_soil[i1,i2,:,:],aggr_layers)
        modes_plotting[ind,:,:] = mds
        plot!(mds[:,1], mds[:,2], legend=:topright, color = i2, label = string("Age ",depth_function[i2]))
        plot!(mds[:,3], mds[:,2], linestyle=:dash, color = i2, label = string("NSF ",depth_function[i2]))
    end
    
end
current()
plot()




# Compare PDFs mixing with depth
indices = 2:5:100
lgth = length(indices)
P=plot()
cols = ["red", "blue", "green"]

pdf_all=zeros(3,lgth,2000,2)

for df in 1:3
    pdf_lines = zeros(2000,length(indices)+1)

    for l in 1:length(indices)
        data = copy(output_ages[2,df,indices[l]])
        filter!(e->e≠10000,data)
        pdf = KDE(data)
        # pdf[:,2] /= maximum(pdf[:,2])
        pdf_all[df,l,:,:] = pdf

        if(l==1)
            pdf_lines[:,1] = pdf[:,1]
        end
        pdf_lines[:,l+1] = pdf[:,2]
    end
    plot!(P, pdf_lines[:,1],pdf_lines[:,2:lgth+1], 
    lc = reshape( range(parse(Colorant, cols[df]), stop=colorant"white",length=lgth), 1, lgth ))
end


p1 = plot(pdf_all[1,1,:,1], pdf_all[1,:,:,2]',lc = reshape( range(parse(Colorant, cols[1]), stop=colorant"white",length=lgth), 1, lgth ))
p2 = plot(pdf_all[2,1,:,1], pdf_all[2,:,:,2]',lc = reshape( range(parse(Colorant, cols[2]), stop=colorant"white",length=lgth), 1, lgth ))
p3 = plot(pdf_all[3,1,:,1], pdf_all[3,:,:,2]',lc = reshape( range(parse(Colorant, cols[3]), stop=colorant"white",length=lgth), 1, lgth ))
plot(p1,p2,p3, layout = (1,3))

display(P)


p1=plot(pdf_lines[:,1],pdf_lines[:,2:lgth+1], 
lc = reshape( range(colorant"red", stop=colorant"white",length=lgth), 1, lgth ))




plot(ylims = [0,1], xlims=[0,10000],cbar=true)

plot(pdf_lines[:,1],pdf_lines[:,2:lgth+1], line_z=1:lgth,
lc = reshape( range(colorant"red", stop=colorant"white",length=lgth), 1, lgth ))




p=plot(pdf_lines[:,2:lgth+1], line_z = (0:lgth)',
lc = reshape( range(colorant"red", stop=colorant"white",length=lgth), 1, lgth ))





p
f(x; b) = x.^2 .+ b
x = collect(0:10)
Y = [f(x, b = b) for b in -10:2:10]
p


plot!plot(legend=false, ratio=1, xlims=extrema(x),cbar=true)


plot!(p, pdf, col = CList[Int32((l-1)/5+1)])



using DataFrames, StatsPlots, Distributions 
x = 1:10
linear(b0, b1, x) = b0 + b1*x
B0 = [0,10]
B1 = range(.5, stop=1, length=5)
y = [linear.(b0, b1, x) for b1 in B1 for b0 in B0]
y = vcat(y...)
df = DataFrame(group = repeat(0:1,inner=50), slope = repeat(repeat(B1, inner=10), outer =2),
    x = repeat(x,outer=10),  y = y)
@df df plot(:x, :y, group = (:group,:slope), legend=false, palette = cgrad(:viridis,:slope))


XX, YY, SL = zeros(10,10), zeros(10,10), zeros(10)
i = 1
for d in groupby(df,[:group,:slope])
    XX[:,i], YY[:,i], SL[i] = d.x, d.y, d.slope[1]
    i += 1
end
I = sortperm(SL)
plot(legend=false, ratio=1, xlims=extrema(x),cbar=true)
plot!(XX[I,:], YY[I,:], line_z = SL[I]')




p



# Exp. 4
upheaval_frequencies = [100,1000,2000,5000]
depth_function = "grd"
aggr_layers = 5
wd = "C:\\Simulations\\Bioturbation\\Julia\\"
output_soil = load(string(wd,"exp4_output_soil.jld2"))["output_soil"]
output_ages = load(string(wd,"exp4_output_ages.jld2"))["output_ages"]
nlayers = size(output_soil)[2]


for i in 1:4
    
    write_CSV_output(string("C:\\Simulations\\Bioturbation\\Julia\\exp4_f",upheaval_frequencies[i],"_"),
    output_soil[i,:,:], output_ages[i,:,:])
end




p=plot(0, xlims=[0,10000], ylims = [0,1.5], legend=:bottomleft); yflip!(true)
for freq in 1:length(upheaval_frequencies)
    modes = get_modes(output_ages[freq,:], output_soil[freq,:,:],aggr_layers)
    plot!(p, modes[:,1], modes[:,2], legend=:topright, color = freq, label = string("Age ",upheaval_frequencies[freq]))
    plot!(p, modes[:,3], modes[:,2], linestyle=:dash, color = freq, label = string("NSF ",upheaval_frequencies[freq]))
end
p


# Exp 5
ratios = [1,1,1,1,1,0,0,0,0,0]
depth_function = "grd"
wd = "C:\\Simulations\\Bioturbation\\Julia\\"
output_soil = load(string(wd,"exp5_output_soil.jld2"))["output_soil"]
output_ages = load(string(wd,"exp5_output_ages.jld2"))["output_ages"]
nlayers = size(output_soil)[2]
aggr_layers=5

p=plot(0, xlims=[0,10000], ylims = [0,1.5], legend=:bottomleft); yflip!(true)
for i in 1:length(ratios)
    modes = get_modes(output_ages[i,:], output_soil[i,:,:],aggr_layers)
    plot!(p, modes[:,1], modes[:,2], legend=:topright, color = i, label = string("Age ",ratios[i]))
    plot!(p, modes[:,3], modes[:,2], linestyle=:dash, color = i, label = string("NSF ",ratios[i]))
end
p

#=for i in length(ratios)
    if i==1
        for l in 1:200
            
    else
    end
end
=#


# Compare distributions
kde1 = KDE(output_ages[5,2])
kde2 = KDE(output_ages[5,85])


kde1=hcat(kde1,kde2[:,2])
kde1 = hcat(kde1,zeros(2000,1))

for i in 1:2000
    kde1[i,4] = minimum(kde1[i,2:3])
end

plot(kde1[:,1], kde1[:,2:4], lw  = [2 2 1])


sum(sum(kde1[:,4])/sum(kde1[:,2]) + sum(kde1[:,4])/sum(kde1[:,3]))/2
