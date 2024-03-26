
function eval_ver4O0ParquetAD000!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O0ParquetAD.so").eval_graph000(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end