include("HW7.jl")
l_0 = 0.3
l_1 = 0.6
l_2 = 0.4
## Answer for 3.b
twist_1 = [0 0 1 0 0 0]';
twist_2 = [0 0 1 l_1 0 0]';
twist_3 = [0 0 1 l_1+l_2 0 0]';
twist_4 = [0 0 0 0 0 1]';
twists = [twist_1 twist_2 twist_3 twist_4];
R = Matrix(I,3,3)
t = [0 l_1+l_2 l_0];
T=RpToTrans(R,t)
angles=[pi/3 -pi/2 pi/4 -0.35]
## Answer for 3.b
g = FwdKin(T,twists,angles)
## Answer for 3.c
Js = JacobianSpace(twists,angles)
## Answer for 3.d
l_rod = 0.15;
P_rod = [0 l_1+l_2+l_rod l_0 1]';
Q_dot = [-0.2 1 0.5 -2]';

V_s = JacobianSpace(twists,angles)*Q_dot;
V_rod = VecTose3( V_s ) * P_rod;
