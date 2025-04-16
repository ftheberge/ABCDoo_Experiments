using Pkg
using ABCDGraphGenerator
using Random

youtube_like_params = Dict(
    "name"  => "yt",
    "nout"  => 0,
    "d_min" => 5,
    "d_max" => 100,
    "τ₁"    => 1.87,       
    "c_min" => 50,
    "c_max" => 500,
    "τ₂"    => 2.13,
    "ρ"     => 0.37,
)

# Universal params
seed = 123
d_max_iter = 1000
c_max_iter = 1000
Random.seed!(seed)

for params in [youtube_like_params]
    for d in [5]
        @info "$(params["name"]), d=$d"

        name = params["name"]
        nout = params["nout"]
        d_min = params["d_min"]
        d_max = params["d_max"]
        τ₁ = params["τ₁"]
        c_min = params["c_min"]
        c_max = params["c_max"]
        τ₂ = params["τ₂"]
        ρ = params["ρ"]

        ## first loop - fix n=5000
        for η in [1.0,1.25,1.5,1.75,2.0,2.25,2.5]
            for ξ in [0.1,0.2,0.3,0.4,0.5,0.6]
                for n in [5000]
                    # in what follows n is number of non-outlier nodes
                    n = n - nout
                    # Actually Generate the graph
                    @info "Expected value of degree: $(ABCDGraphGenerator.get_ev(τ₁, d_min, d_max))"
                    degs = ABCDGraphGenerator.sample_degrees(τ₁, d_min, d_max, n + nout, d_max_iter)
                    @assert iseven(sum(degs))
                    @info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"
                    coms = ABCDGraphGenerator.sample_communities(τ₂, ceil(Int, c_min / η), floor(Int, c_max / η), n, c_max_iter)
                    @assert sum(coms) == n
                    pushfirst!(coms, nout)
                    @info "    Done degs and coms, generating graph."
                    p = ABCDGraphGenerator.ABCDParams(degs, coms, ξ, η, d, ρ)
                    edges, clusters = ABCDGraphGenerator.gen_graph(p)
                    open("abcdoo_$(ξ)_$(η)_$(n)_edge.dat", "w") do io
                    #open("abcdoo_$(name)_d$(d)_nocorr_edge.dat", "w") do io
                       for (a, b) in sort!(collect(edges))
                           println(io, a, "\t", b)
                       end
                    end
                    open("abcdoo_$(ξ)_$(η)_$(n)_com.dat", "w") do io
                    #open("abcdoo_$(name)_d$(d)_nocorr_com.dat", "w") do io
                        for (i, c) in enumerate(clusters)
                            println(io, i, "\t", c)
                        end
                    end
                end
            end
        end

        ## second loop - vary n only
        for η in [1.25]
            for ξ in [0.1]
                for n in [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000]
                    # in what follows n is number of non-outlier nodes
                    n = n - nout
                    # Actually Generate the graph
                    @info "Expected value of degree: $(ABCDGraphGenerator.get_ev(τ₁, d_min, d_max))"
                    degs = ABCDGraphGenerator.sample_degrees(τ₁, d_min, d_max, n + nout, d_max_iter)
                    @assert iseven(sum(degs))
                    @info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"
                    coms = ABCDGraphGenerator.sample_communities(τ₂, ceil(Int, c_min / η), floor(Int, c_max / η), n, c_max_iter)
                    @assert sum(coms) == n
                    pushfirst!(coms, nout)
                    @info "    Done degs and coms, generating graph."
                    p = ABCDGraphGenerator.ABCDParams(degs, coms, ξ, η, d, ρ)
                    edges, clusters = ABCDGraphGenerator.gen_graph(p)
                    open("abcdoo_$(ξ)_$(η)_$(n)_edge.dat", "w") do io
                    #open("abcdoo_$(name)_d$(d)_nocorr_edge.dat", "w") do io
                       for (a, b) in sort!(collect(edges))
                           println(io, a, "\t", b)
                       end
                    end
                    open("abcdoo_$(ξ)_$(η)_$(n)_com.dat", "w") do io
                    #open("abcdoo_$(name)_d$(d)_nocorr_com.dat", "w") do io
                        for (i, c) in enumerate(clusters)
                            println(io, i, "\t", c)
                        end
                    end
                end
            end
        end
        
    end
end
