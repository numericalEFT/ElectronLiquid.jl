
function eval_spinPolar_GV100!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph100(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV110!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph110(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV120!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph120(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV130!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph130(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV140!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph140(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV200!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph200(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV201!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph201(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV202!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph202(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV203!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph203(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV210!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph210(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV211!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph211(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV212!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph212(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV220!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph220(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV221!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph221(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV230!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph230(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV300!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph300(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV301!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph301(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV302!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph302(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV310!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph310(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV311!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph311(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV320!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph320(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV400!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph400(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV401!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph401(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV410!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph410(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end
function eval_spinPolar_GV500!(root::Vector{Float64}, leafVal::Vector{Float64})
    @ccall joinpath(@__DIR__, "spinPolar_GV.so").eval_graph500(root::Ptr{Cdouble}, leafVal::Ptr{Cdouble})::Cvoid
end