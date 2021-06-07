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
    else
        theta = acos(acosinput)
        return 1 / 2 / sin(theta) * (R - R')
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
    t = T[1:3,4]
    omg_h = MatrixLog3(R)
    w = [omg_h[3,2];omg_h[1,3];omg_h[2,1]]
    if omg_h == zeros(3, 3)
        return vcat(hcat(zeros(3, 3), T[1:3, 4]), [0 0 0 0])
    else
        theta = acos((lAlgebra.tr(R) - 1) / 2)
        return vcat(hcat(omg_h,inv((I-R)*omg_h + w*w'*theta)*t),
                          [0 0 0 0])
    end
end
#****************************************************************#
