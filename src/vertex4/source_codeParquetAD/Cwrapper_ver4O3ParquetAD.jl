
function eval_ver4O3ParquetAD102!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O3ParquetAD.so").eval_graph102(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O3ParquetAD201!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O3ParquetAD.so").eval_graph201(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O3ParquetAD200!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O3ParquetAD.so").eval_graph200(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O3ParquetAD210!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O3ParquetAD.so").eval_graph210(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O3ParquetAD101!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O3ParquetAD.so").eval_graph101(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O3ParquetAD300!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O3ParquetAD.so").eval_graph300(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_ver4O3ParquetAD100!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "ver4O3ParquetAD.so").eval_graph100(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end