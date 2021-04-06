using LinearAlgebra

function so3ToVec(so3mat)
    if( size(so3mat)!=(3,3) )
        println("!Warning! : Input matrix should be 3x3.")
        return
    end
    if( !((so3mat[1,1]==so3mat[2,2]==so3mat[3,3]==0)) )
        println("!Warning! : Diagonals should be zero.")
        return
    end
    if( so3mat[2,1]!=-so3mat[1,2] ||
        so3mat[3,1]!=-so3mat[1,3] ||
        so3mat[3,2]!=-so3mat[2,3] )
        println("!Warning! : Entered matrix isn't a skew symmetric matrix.")
        return
    end
    ret_vec = [so3mat[3,2];
               so3mat[1,3];
               so3mat[2,1] ]
        return ret_vec
end

function vecToSo3(omg)
    if( size(omg)!=(3,) )
        println("!Warning! : Input matrix should be 3x1.")
        return
    end

    ret_ss_matrix = [0      -omg[3]   omg[2];
                     omg[3]    0     -omg[1];
                    -omg[2]  omg[1]      0];
    return ret_ss_matrix
end
