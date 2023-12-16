
function eval_ver4O2ParquetAD200!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O2ParquetAD.so").eval_graph200(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O2ParquetAD101!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O2ParquetAD.so").eval_graph101(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O2ParquetAD100!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O2ParquetAD.so").eval_graph100(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end