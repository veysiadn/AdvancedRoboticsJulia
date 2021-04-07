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
    if( size(w)!=(3,) )
        println("!Warning! : W matrix should be 3x1.")
        return
    end
    Id_3x3=Diagonal(ones(3,3));
    w_hat = [0      -w[3]    w[2];
             w[3]     0     -w[1];
            -w[2]    w[1]      0];
    Rotation_matrix = cos(theta)*Id_3x3 + sin(theta)*w_hat+(1-cos(theta))*w*w';
    return Rotation_matrix;
end
#****************************************************************#
# Builds rigid transformation matrix T from given
#rotation matrix R and translation t
function RpToTrans(R,t)
    if( size(R)!=(3,3) )
        println("!Warning! : Rotation matrix should be 3x3.")
        return
    end
    if( size(t)!=(3,) )
        println("!Warning! : Translation matrix should be 3x1.")
        return
    end
    temp =  hcat(R,t) ;
    row_4 = [0 0 0 1] ;
    R_B_T = vcat(temp,row_4);
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
        println("!Warning! : W matrix should be 3x1.")
        return
    end
    if( det(T)==0 )
        println("!Warning! : Inverse doesn't exist (Determinant is equal to zero).")
        return
    end
    return T^-1;
end
#****************************************************************#
