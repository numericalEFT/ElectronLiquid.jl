
function eval_graph100!(root::AbstractVector, leafVal::AbstractVector)
     g1 = leafVal[1]
     g2 = leafVal[2]
     g3 = leafVal[3]
     g4 = leafVal[4]
     g5 = (g1 * g2 * g3 * g4) * -1.0
     g7 = (g5 * -1.0)
     root[1] = g7
 end 