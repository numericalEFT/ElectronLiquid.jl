
function eval_freeEnergy_ParquetAD000!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph000(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD010!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph010(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD020!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph020(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD030!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph030(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD040!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph040(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD050!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph050(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD100!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph100(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD101!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph101(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD102!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph102(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD103!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph103(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD104!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph104(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD110!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph110(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD111!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph111(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD112!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph112(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD113!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph113(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD120!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph120(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD121!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph121(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD122!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph122(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD130!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph130(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD131!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph131(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD140!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph140(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD200!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph200(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD201!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph201(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD202!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph202(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD203!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph203(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD210!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph210(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD211!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph211(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD212!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph212(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD220!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph220(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD221!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph221(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD230!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph230(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD300!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph300(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD301!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph301(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD302!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph302(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD310!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph310(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD311!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph311(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD320!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph320(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD400!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph400(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD401!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph401(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD410!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph410(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_freeEnergy_ParquetAD500!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "freeEnergy_ParquetAD.so").eval_graph500(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end