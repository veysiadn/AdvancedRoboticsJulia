#****************************************************************#
# VeysiADN 17 Apr 2021
# Advanced Robotics Homework 4 Codes
#****************************************************************#
using LinearAlgebra
const lAlgebra = LinearAlgebra
#****************************************************************#
# Returns normalized version of given matrix/array
function Normalize(N)
    return N/lAlgebra.norm(N)
end
#****************************************************************#
# Returns 4x4 twist matrix ‘se3mat’ (an element of Lie Algebra se3)
# from given twist vector (= 6x1 vector) ‘xi
function VecToSe3( xi )
    w_skew = [ 0   -xi[3]   xi[2];
              xi[3]   0    -xi[1];
             -xi[2]  xi[1]  0   ]
    return vcat(hcat(w_skew,xi[4:6]), zeros(1, 4))
end
#****************************************************************#
# Returns twist vector ‘xi’ from given twist matrix ‘se3mat’
function Se3ToVec( se3mat )
    omg = [se3mat[3, 2], se3mat[1, 3],se3mat[2, 1]]
    v = se3mat[1:3, 4]
    return vcat(omg,v)
end
#****************************************************************#
# Calculates and return rotation matrix ‘R’ by
# evaluating matrix exponential of 3x3 skew symmetric matrix ‘so3mat’

function MatrixExp3(so3mat)
        omg = [so3mat[3, 2], so3mat[1, 3],so3mat[2, 1]]
        theta = lAlgebra.norm(omg)
        omg_n = so3mat / theta
        return lAlgebra.I + sin(theta) * omg_n + (1 - cos(theta)) * omg_n * omg_n
end

#****************************************************************#
# Calculates and return transformation matrix ‘g’ by
# evaluating matrix exponential of 4x4 twist matrix ‘se3mat’
function MatrixExp6( se3mat )
    omg = [se3mat[3, 2], se3mat[1, 3],se3mat[2, 1]]
    theta = lAlgebra.norm(omg)
    omg_n = se3mat[1:3, 1:3] / theta
    temp  = hcat(MatrixExp3(se3mat[1:3, 1:3]),
                         (lAlgebra.I * theta +
                          (1 - cos(theta)) * omg_n +
                          (theta - sin(theta)) * omg_n * omg_n) *
                         se3mat[1:3, 4] / theta)
    return    vcat(temp,[0 0 0 1])
end
