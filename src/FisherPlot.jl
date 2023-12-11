module FisherPlot

using LaTeXStrings
using Makie
using Base: @kwdef

function ellipseparameterization(a::Float64, b::Float64, θ::Float64)
    t = LinRange(0,2π, 200)
    x = Array(a .* sin.(θ) .* cos.(t) + b .* cos.(θ) .* sin.(t))
    y = Array(a .* cos.(θ) .* cos.(t) - b .* sin.(θ) .* sin.(t))
    return x, y
end

function gaussian(μ::Float64, σ::Float64, x)
    return 1/(sqrt(2π*σ^2))*exp(-0.5*(x-μ)^2/σ^2)
end

function ellipseparameters(covmatrix::Matrix{Float64}, i::Int64, j::Int64)
    σi = sqrt(covmatrix[i,i])
    σj = sqrt(covmatrix[j,j])
    σij = covmatrix[i,j]
    θ = (atan(2σij,(σi^2-σj^2)))/2
    a = sqrt((σi^2+σj^2)/2+sqrt(((σi^2-σj^2)^2)/4+σij^2))
    if i == j
        b = 0.
    else
        b = sqrt((σi^2+σj^2)/2-sqrt(((σi^2-σj^2)^2)/4+σij^2))
    end
    return σi, σj, a, b, θ
end

function preparecanvas(LaTeX_array, limits, ticks, probes, colors, PlotPars::Dict)
    matrix_dim = length(LaTeX_array)
    #TODO: add the textsize to PlotPars
    figure = Makie.Figure(textsize = 40, font = PlotPars["font"],
        backgroundcolor=PlotPars["backgroundcolor"])
    ga = figure[1, 1] = GridLayout()

    titlelayout = GridLayout(figure[0, 1], halign = :center, tellwidth = false)

    if haskey(PlotPars, "plottitle")
        Label(titlelayout[1, 1], PlotPars["plottitle"], fontsize = PlotPars["plottitlesize"], halign = :center, font = PlotPars["plottitlefont"])
    end
    for i in 1:matrix_dim
        for j in 1:i
            if i == j
                ax = Axis(ga[i,i],
                    width = PlotPars["sidesquare"], height = PlotPars["sidesquare"],
                    xticklabelsize = PlotPars["dimticklabel"],
                    yticklabelsize = PlotPars["dimticklabel"], yaxisposition = (:right),
                    xlabel = L"%$((LaTeX_array)[i])", xlabelsize = PlotPars["parslabelsize"],
                    ylabelsize = PlotPars["PPmaxlabelsize"], spinewidth=3., xtickwidth = 3.0, ytickwidth = 3.0, xticksize=10., yticksize=10.,  yticks = [0,1], ylabel = L"P/P_\mathrm{max}",
                    xticklabelrotation = PlotPars["xticklabelrotation"],
                    xticks = ([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]],
                        [string(myi) for myi in round.([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]], sigdigits = 3)]),
                        backgroundcolor=PlotPars["backgroundcolor"], alignmode = Inside())
                Makie.ylims!(ax, (-0.0,1.05))
                Makie.xlims!(ax, (limits[i,1],limits[i,2]))
                Makie.hideydecorations!(ax, ticks = false, ticklabels = false, label = false)

                ax.alignmode = Mixed(right = -102)
                if i != matrix_dim
                    hidexdecorations!(ax, ticks = true, ticklabels = true,  label = true)
                else
                    hidexdecorations!(ax, ticks = false, ticklabels = false,  label = false)
                end
            else
                ax = Axis(ga[i,j], width = PlotPars["sidesquare"], height = PlotPars["sidesquare"],spinewidth=3., xtickwidth = 3.0, ytickwidth = 3.0, xticksize=10., yticksize=10.,
                    xticklabelsize = PlotPars["dimticklabel"], yticklabelsize = PlotPars["dimticklabel"],
                    ylabel = L"%$(LaTeX_array[i])", xlabel = L"%$((LaTeX_array)[j])",
                    ylabelsize = PlotPars["parslabelsize"], xlabelsize = PlotPars["parslabelsize"], xticklabelrotation = PlotPars["xticklabelrotation"],
                    yticks = ([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]],
                        [string(myi) for myi in round.([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]], sigdigits = 3)]),
                    xticks = ([ticks[j,1], 0.5*(ticks[j,1]+ticks[j,2]), ticks[j,2]],
                        [string(myi) for myi in round.([ticks[j,1], 0.5*(ticks[j,1]+ticks[j,2]), ticks[j,2]], sigdigits = 3)]),
                        yticklabelpad=8, backgroundcolor=PlotPars["backgroundcolor"], alignmode = Inside())
                Makie.ylims!(ax, (limits[i,1],limits[i,2]))
                Makie.xlims!(ax, (limits[j,1],limits[j,2]))
                if i == matrix_dim
                    hidexdecorations!(ax, ticks = false, ticklabels = false,  label = false)
                else
                    hidexdecorations!(ax, ticks = true, ticklabels = true,  label = true)
                end
                if j == 1
                    hideydecorations!(ax, ticks = false, ticklabels = false,  label = false)
                    Legend(ga[1,matrix_dim],
                    [PolyElement(color = color, strokecolor = color, strokewidth = 1) for color in colors],
                    probes, padding = PlotPars["legendpadding"], patchlabelgap = PlotPars["patchlabelgap"],
                    tellheight = false, tellwidth = false, rowgap = 0,
                    halign = :right, valign = :top, framecolor = :black, labelsize =PlotPars["legendsize"], patchsize = PlotPars["patchsize"], framevisible = true,
                    backgroundcolor=PlotPars["backgroundcolor"])
                else
                    hideydecorations!(ax, ticks = true, ticklabels = true,  label = true)
                end
            end
        end
    end
    resize_to_layout!(figure)
    trim!(figure.layout)
    return figure
end

function drawgaussian!(canvas, σ, i, central, color)
    ax = canvas[i,i]
    x = Array(LinRange(-4σ+central,4σ+central, 200))
    Makie.lines!(ax, x, gaussian.(central, σ, x)./gaussian.(central, σ, central), color = color, linewidth = 4)
    Makie.band!(ax, x, 0, gaussian.(central, σ, x)./gaussian.(central, σ, central) , color=(color, 0.2))
    x = Array(LinRange(-σ+central,σ+central, 200))
    Makie.band!(ax, x, 0, gaussian.(central, σ, x)./gaussian.(central, σ, central) , color=(color, 0.4))
end

function drawellipse!(canvas, i, j, x, y, central_values, color)
    ax = canvas[i,j]

    Makie.lines!(ax, x .+ central_values[j], y .+ central_values[i], color = color, linewidth = 4)
    Makie.lines!(ax, 3x .+ central_values[j], 3y .+ central_values[i], color = color, linewidth = 4)

    Makie.band!(ax, x .+ central_values[j], 0, y .+ central_values[i], color=(color, 0.4))
    Makie.band!(ax, 3x .+ central_values[j], 0, 3y .+ central_values[i], color=(color, 0.2))
end

function paintcorrmatrix!(canvas, central_values, corr_matrix, color)
    ciccio = canvas[1,1]
    for i in 1:length(central_values)
        for j in 1:i
            if i == j
                drawgaussian!(ciccio, sqrt(corr_matrix[i,i]), i, central_values[i], color)
            else
                σi, σj, a, b, θ = ellipseparameters(corr_matrix, i,j)
                x,y = ellipseparameterization(a, b, θ)
                drawellipse!(ciccio, i, j, x, y, central_values, color)
            end
        end
    end
end

function save(filename, obj)
    Makie.save(filename, obj)
end

end # module
