using Random, BenchmarkTools

using HighVoronoi
# Note: HyperVoronoiDelaunay is assumed available but unconfirmed API

using DelaunayTriangulation
using Quickhull

# Generate random points
function main()
    ρ = 1.0
    Random.seed!(1234)

    for  N = [500, 5_000]
        L = sqrt(N/ρ)
        points = rand(2, N)*L  # 2×N array of Float64

        println("Benchmarking Voronoi workflows with N = $N points...\n")

        # HighVoronoi.jl — Periodic Voronoi
        println("--- HighVoronoi.jl (periodic) ---")
        @btime begin
            VG = VoronoiGeometry(VoronoiNodes($ points ./ $L), cuboid(2, periodic=[1, 2]), silence=true);
            VD = VoronoiData(VG);
        end


        # DelaunayTriangulation.jl — Planar Voronoi
        println("\n--- DelaunayTriangulation.jl ---")
        @btime begin
            tri = triangulate($points);
            vorn = voronoi(tri);
        end

        # Quickhull.jl — Delaunay + Voronoi via Quickhull
        println("\n--- Quickhull.jl ---")
        @btime begin
            tri_qh = delaunay($points);
            edges = voronoi_edge_points(tri_qh);
        end
    end
end

main()