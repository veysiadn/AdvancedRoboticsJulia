#****************************************************************#
# VeysiADN 12 May 2021
# Advanced Robotics Homework 5 Codes
#****************************************************************#
using LinearAlgebra
const lAlgebra = LinearAlgebra
#****************************************************************#
#Computes exponential coordinate ‘so3mat’
#(skew symmetric matrix) of the rotation matrix R

function MatrixLog3(R)
    acosinput = (lAlgebra.tr(R) - 1) / 2
    if acosinput >= 1
        return zeros(3, 3)
    elseif acosinput <= -1
        if !NearZero(1 + R[3, 3])
            omg = (1 / sqrt(2 * (1 + R[3, 3]))) * [R[1, 3], R[2, 3], 1 + R[3, 3]]
        elseif !NearZero(1 + R[2, 2])
            omg = (1 / sqrt(2 * (1 + R[2, 2]))) * [R[1, 2], 1 + R[2, 2], R[3, 2]]
        else
            omg = (1 / sqrt(2 * (1 + R[1, 1]))) * [1 + R[1, 1], R[2, 1], R[3, 1]]
        end
        ret_ss_matrix = pi * [0      -omg[3]   omg[2];
                             omg[3]    0       -omg[1];
                             -omg[2]  omg[1]      0];
        return ret_ss_matrix
    else
        theta = acos(acosinput)
        return theta / 2 / sin(theta) * (R - R')
    end
end
#****************************************************************#
# Computes exponential coordinate ‘se3mat’
# (twist matrix) of the rigid body transformation matrix T
function MatrixLog6(T)
    if( size(T)!=(4,4) )
        println("!Warning! : Transformation matrix should be 4x4.")
        return
    end
    R = T[1:3,1:3]
    p = T[1:3,4]
    omgmat = MatrixLog3(R)
    if omgmat == zeros(3, 3)
        return vcat(hcat(zeros(3, 3), T[1:3, 4]), [0 0 0 0])
    else
        theta = acos((lAlgebra.tr(R) - 1) / 2)
        return vcat(hcat(omgmat,
                         (lAlgebra.I - omgmat / 2 +
                          (1 /theta - 1 / tan(theta / 2) / 2) *
                          omgmat * omgmat / theta) * T[1:3, 4]),
                          [0 0 0 0])
    end
end
#****************************************************************#
