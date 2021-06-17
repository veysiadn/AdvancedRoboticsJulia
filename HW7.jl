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

function VecToso3(omg)
    if( size(omg)!=(3,) )
        println("!Warning! : Input matrix should be 3x1.")
        return
    end

    ret_ss_matrix = [0      -omg[3]   omg[2];
                     omg[3]    0     -omg[1];
                    -omg[2]  omg[1]      0];
    return ret_ss_matrix
end
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
function VecTose3( xi )
    w_skew = [ 0   -xi[3]   xi[2];
              xi[3]   0    -xi[1];
             -xi[2]  xi[1]    0 ]
    return vcat(hcat(w_skew,xi[4:6]), zeros(1, 4))
end
#****************************************************************#
# Returns twist vector ‘xi’ from given twist matrix ‘se3mat’
function se3ToVec( se3mat )
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
    theta = norm(omg)
    if theta==0
        return vcat(hcat(I, se3mat[1:3, 4]), [0 0 0 1])
    else
    omg_n = se3mat[1:3, 1:3] / theta
    temp  = hcat(MatrixExp3(se3mat[1:3, 1:3]),
                         (I * theta +
                          (1 - cos(theta)) * omg_n +
                          (theta - sin(theta)) * omg_n * omg_n) *
                         se3mat[1:3, 4] / theta)
    end
    return    vcat(temp,[0 0 0 1])
end

#****************************************************************#
# VeysiADN 12 May 2021
# Advanced Robotics Homework 5 Codes
#****************************************************************#
using LinearAlgebra
const lAlgebra = LinearAlgebra
#****************************************************************#
# Computes exponential coordinate ‘so3mat’
# (skew symmetric matrix) of the rotation matrix R

function MatrixLog3(R)
    acosinput = (lAlgebra.tr(R) - 1) / 2
    if acosinput >= 1
        return zeros(3, 3)
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
    t = T[1:3,4]
    omg_h = MatrixLog3(R)
    w = [omg_h[3,2];omg_h[1,3];omg_h[2,1]]
    if omg_h == zeros(3, 3)
        return vcat(hcat(zeros(3, 3), T[1:3, 4]), [0 0 0 0])
    else
        theta = acos((lAlgebra.tr(R) - 1) / 2)
        return vcat(hcat(omg_h,
                         (I - omg_h / 2 +
                          (1 / theta - 1 / tan(theta / 2) / 2) *
                          omg_h * omg_h / theta) * T[1:3, 4]),
                    [0 0 0 0])
    end
end
#****************************************************************#

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
#  Starting from here  is homework 7.
"""
    Calculates Spatial geometric Jacobian Js, given 6 x n matrix
    ‘Twists’ whose columns are Twist parameter of each joints
    and n x 1 vector of joint angles
"""
function JacobianSpace(Twists, JointAngles)
    T = I;
    Js = Twists
    for i = 2:length(JointAngles)
        T *= MatrixExp6(VecTose3(Twists[:, i - 1] * JointAngles[i - 1]))
        Js[:, i] = Adjoint(T) * Twists[:, i]
    end
    return Js
end
"""
Calculates Forward Kinematics given initial configuration M,
6 x n matrix ‘Twists’ whose columns are twist parameter
of each joints, and n x 1 joint angles
"""
function FwdKin(M,Twists,JointAngles)
    T=M
    for i = length(JointAngles):-1:1
        T = MatrixExp6(VecTose3(Twists[:,i]*JointAngles[i])) * T
    end
    return T
end
"""
Solves Paden-Kahan subproblem 1

"""
function PadenKahanFirst(w, x, y)
    p = x
    q = y
    u = p - (w'*p)*w;
    v = q - (w'*q)*w;
    theta = atan((w'*(cross(u,v)))/(u'*v))
    return theta
end
"""
Solves Paden-Kahan subproblem 2. th will be 2 x 2 matrix
containing solution pair (th1, th2) each row
"""

function PadenKahanSecond(w1, w2, x, y)
    p=x
    q=y

    A = [1 w1'*w2;
        w1'*w2 1];
    B=[w1'*q; w2'*p]
    alfa,beta = A\B

    gamma_sq = ((norm(p)^2)-(alfa^2)-(beta^2)-(2*alfa*beta*(w1'*w2)))/((norm(cross(w1,w2)))^2)
    gamma_1 = sqrt.(gamma_sq);
    gamma_2 = -sqrt.(gamma_sq);
    c_1 = alfa*w1+beta*w2+gamma_1*(cross(w1,w2));
    c_2 = alfa*w1+beta*w2+gamma_2*(cross(w1,w2));
    th_1 = PadenKahanFirst(w2,p,c_1)
    th_2 = PadenKahanFirst(w2,p,c_2)
    th_3 = PadenKahanFirst(-w1,q,c_1)
    th_4 = PadenKahanFirst(-w1,q,c_2)
    th = [th_1 th_3;th_4 th_2];
    return th
end
