
function eval_ver4O1ParquetAD100!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O1ParquetAD.so").eval_graph100(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end