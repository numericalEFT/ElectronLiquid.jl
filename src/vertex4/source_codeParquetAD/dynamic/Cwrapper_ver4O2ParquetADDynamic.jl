
function eval_ver4O2ParquetADDynamic200!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O2ParquetADDynamic.so").eval_graph200(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O2ParquetADDynamic101!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O2ParquetADDynamic.so").eval_graph101(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O2ParquetADDynamic100!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O2ParquetADDynamic.so").eval_graph100(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end