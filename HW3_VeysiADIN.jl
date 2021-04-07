#****************************************************************#
# VeysiADN 4 Apr 2021
# Advanced Robotics Homework 3 Codes
#****************************************************************#
using LinearAlgebra
#****************************************************************#
# Returns 3 x 3 rotation matrix representing
# rotation about unit vector w by angle theta.
# theta should be in radian
function Rodrigues(w,theta)
    if( length(w) != 3)
        println("!Warning! : W matrix should have 3 elements.")
        return
    end
    if(size(w)==(1,3))
        w=w';
    end
    w_hat = [0      -w[3]    w[2];
             w[3]     0     -w[1];
            -w[2]    w[1]      0];
    rotation_matrix = cos(theta)*Matrix(1I,3,3) + sin(theta)*w_hat+(1-cos(theta))*w*w';
    return rotation_matrix;
end
#****************************************************************#
# Builds rigid transformation matrix T from given
#rotation matrix R and translation t
function RpToTrans(R,t)
    if( size(R)!=(3,3) )
        println("!Warning! : Rotation matrix should be 3x3.")
        return
    end
    if( length(t) != 3 )
        println("!Warning! : Translation matrix should have 3 elements.")
        return
    end
    if(size(t)==(1,3))
        t=t';
    end
    R_B_T = hcat(R,t) ;
    R_B_T = vcat(R_B_T,[0 0 0 1]) ;
    return R_B_T
end
#****************************************************************#
#Extracts rotation matrix R and translation t from given /
# rigid transformation matrix T
function TransToRp(T)
    if( size(T)!=(4,4) )
        println("!Warning! : Transformation matrix should be 4x4.")
        return
    end
    R = T[1:3,1:3]
    t = T[1:3,4]
    return R,t;
end
#****************************************************************#
# Inverts rigid transformation matrix T
function TransInv(T)
    if( size(T)!=(4,4) )
        println("!Warning! : W matrix should be 4x4.")
        return
    end
    if( det(T)==0 )
        println("!Warning! : Inverse doesn't exist (Determinant is equal to zero).")
        return
    end
    R,t = TransToRp(T);
    r_inv =hcat(R',-(R'*t))
    r_inv = vcat(r_inv,[0 0 0 1]);
    return r_inv;
end
#****************************************************************#
