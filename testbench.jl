using LinearAlgebra
using BenchmarkTools

function invTest()
    R_inv = zeros(4,4);
    for i = 1:100
        R=rand(Int32,(3,4));
        R=vcat(R,[0 0 0 1]);
        R_inv = R^-1;
    end
return R_inv;
end
function r_invTest()
    r_inv = zeros(4,4);
    for i=1:100
        R = rand(Int32,(3,4));
        R = vcat(R,[0 0 0 1]);
        r_inv = hcat(R[1:3,1:3]',-R[1:3,1:3]' * R[4,1:3]);
        r_inv = vcat(r_inv,[0 0 0 1]);
    end
    return r_inv;
end
