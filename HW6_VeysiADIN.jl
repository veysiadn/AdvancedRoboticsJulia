#****************************************************************#
# VeysiADN 12 May 2021
# Advanced Robotics Homework 5 Codes
#****************************************************************#
using LinearAlgebra
#****************************************************************#
# Returns Adjoint matrix with given rigid body
# transformation matrix g
function Adjoint(g)

    R=g[1:3,1:3];
    p=g[1:3,4];
    p_hat = [ 0       -p[3]     p[2];
             p[3]       0      -p[1];
            -p[2]      p[1]       0 ]

    return [R zeros(3,3);p_hat*R R]
end
