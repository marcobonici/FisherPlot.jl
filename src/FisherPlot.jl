module FisherPlot

using LaTeXStrings
using CairoMakie
using Makie
using Base: @kwdef

function EllipseParametrization(a::Float64, b::Float64, θ::Float64)
    t = LinRange(0,2π, 200)
    x = Array(a .* sin.(θ) .* cos.(t) + b .* cos.(θ) .* sin.(t))
    y = Array(a .* cos.(θ) .* cos.(t) - b .* sin.(θ) .* sin.(t))
    return x, y
end

function Gaussian(μ::Float64, σ::Float64, x)
    return 1/(sqrt(2π*σ^2))*exp(-0.5*(x-μ)^2/σ^2)
end

function EllipseParameters(covmatrix::Matrix{Float64}, i::Int64, j::Int64)
    σi = sqrt(covmatrix[i,i])
    σj = sqrt(covmatrix[j,j])
    σij = covmatrix[i,j]
    θ = (atan(2σij,(σi^2-σj^2)))/2
    a = sqrt((σi^2+σj^2)/2+sqrt(((σi^2-σj^2)^2)/4+σij^2))
    if i==j
        b = 0.
    else
        b = sqrt((σi^2+σj^2)/2-sqrt(((σi^2-σj^2)^2)/4+σij^2))
    end
    return σi, σj, a, b, θ
end

function PrepareCanvas(LaTeX_array, central_values, limits, ticks, probes, colors, PlotPars)
    matrix_dim = length(LaTeX_array)
    #TODO: add the textsize to PlotPars
    figure = Makie.Figure(textsize = 40, font = PlotPars["font"])
    for i in 1:matrix_dim
        for j in 1:i
            if i == j
                ax = Axis(figure[i,i],
                    width = PlotPars["sidesquare"], height = PlotPars["sidesquare"],
                    xticklabelsize = PlotPars["dimticklabel"],
                    yticklabelsize = PlotPars["dimticklabel"], yaxisposition = (:right),
                    xlabel = L"%$((LaTeX_array)[i])", xlabelsize = PlotPars["parslabelsize"],
                    ylabel = L"P/P_\mathrm{max}",ylabelsize = PlotPars["PPmaxlabelsize"], yticks = [0,1],
                    xticklabelrotation = PlotPars["xticklabelrotation"],
                    xticks = ([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]],
                        [string(myi) for myi in round.([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]], sigdigits = 3)]))
                Makie.ylims!(ax, (-0.0,1.05))
                Makie.xlims!(ax, (limits[i,1],limits[i,2]))
                Makie.hideydecorations!(ax, ticks = false, ticklabels = false, label = false)
                if i != matrix_dim
                    ax.alignmode = Mixed(right = MakieLayout.Protrusion(0), bottom = MakieLayout.Protrusion(0), top= MakieLayout.Protrusion(0))
                    hidexdecorations!(ax, ticks = true, ticklabels = true,  label = true)
                else
                    hidexdecorations!(ax, ticks = false, ticklabels = false,  label = false)
                end
            else
                ax = Axis(figure[i,j], width = PlotPars["sidesquare"], height = PlotPars["sidesquare"],
                    xticklabelsize = PlotPars["dimticklabel"], yticklabelsize = PlotPars["dimticklabel"],
                    ylabel = L"%$(LaTeX_array[i])", xlabel = L"%$((LaTeX_array)[j])",
                    ylabelsize = PlotPars["parslabelsize"], xlabelsize = PlotPars["parslabelsize"], xticklabelrotation = PlotPars["xticklabelrotation"],
                    yticks = ([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]],
                        [string(myi) for myi in round.([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]], sigdigits = 3)]),
                    xticks = ([ticks[j,1], 0.5*(ticks[j,1]+ticks[j,2]), ticks[j,2]],
                        [string(myi) for myi in round.([ticks[j,1], 0.5*(ticks[j,1]+ticks[j,2]), ticks[j,2]], sigdigits = 3)]),
                yticklabelpad=8)
                Makie.ylims!(ax, (limits[i,1],limits[i,2]))
                Makie.xlims!(ax, (limits[j,1],limits[j,2]))
                if i == matrix_dim
                    hidexdecorations!(ax, ticks = false, ticklabels = false,  label = false)
                else
                    hidexdecorations!(ax, ticks = true, ticklabels = true,  label = true)
                end
                if j == 1
                    hideydecorations!(ax, ticks = false, ticklabels = false,  label = false)
                    Legend(figure[1,matrix_dim],
                    [PolyElement(color = color, strokecolor = color, strokewidth = 1) for color in colors],
                    probes,
                    tellheight = false, tellwidth = false, rowgap = 10,
                    halign = :right, valign = :top, framecolor = :black, labelsize =55, patchsize = (70, 40), framevisible = true)
                else
                    hideydecorations!(ax, ticks = true, ticklabels = true,  label = true)
                    ax.alignmode = Mixed(right = MakieLayout.Protrusion(0), bottom = MakieLayout.Protrusion(0), top = MakieLayout.Protrusion(0))
                end
            end
        end
    end
    Makie.resize!(figure.scene, figure.layout.layoutobservables.reportedsize[]...)
    return figure
end

function DrawGaussian!(canvas, σ, i, central, color)
    ax = canvas[i,i]
    x = Array(LinRange(-4σ+central,4σ+central, 200))
    Makie.lines!(ax, x, Gaussian.(central, σ, x)./Gaussian.(central, σ, central), color = color, linewidth = 4)
    Makie.band!(ax, x, 0, Gaussian.(central, σ, x)./Gaussian.(central, σ, central) , color=(color, 0.2))
    x = Array(LinRange(-σ+central,σ+central, 200))
    Makie.band!(ax, x, 0, Gaussian.(central, σ, x)./Gaussian.(central, σ, central) , color=(color, 0.4))
end

function DrawEllipse!(canvas, i, j, x, y, central_values, color)
    ax = canvas[i,j]
    
    Makie.lines!(ax, x .+ central_values[j], y .+ central_values[i], color = color, linewidth = 4)
    Makie.lines!(ax, 3x .+ central_values[j], 3y .+ central_values[i], color = color, linewidth = 4)

    Makie.band!(ax, x .+ central_values[j], 0, y .+ central_values[i], color=(color, 0.4))
    Makie.band!(ax, 3x .+ central_values[j], 0, 3y .+ central_values[i], color=(color, 0.2))
end

function PaintCorrMattrix!(canvas, central_values, corr_matrix, color)
    for i in 1:length(central_values)
        for j in 1:i
            if i == j
                DrawGaussian!(canvas, sqrt(corr_matrix[i,i]), i, central_values[i], color)
            else
                σi, σj, a, b, θ = EllipseParameters(corr_matrix, i,j)
                x,y = EllipseParametrization(a, b, θ)
                DrawEllipse!(canvas, i, j, x, y, central_values, color)
            end
        end
    end
end

function save(filename, obj)
    Makie.save(filename, obj)
end

end # module
