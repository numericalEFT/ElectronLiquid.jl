
function eval_ver4O1ParquetAD000!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O1ParquetAD.so").eval_graph000(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O1ParquetAD001!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O1ParquetAD.so").eval_graph001(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O1ParquetAD100!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O1ParquetAD.so").eval_graph100(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end