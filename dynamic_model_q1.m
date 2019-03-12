% 'Robot Dynamics' HW_3 Problem 1
clc;clear all;close all;
syms m1 m2 m3 t1(t) t2(t) t3(t) l1 l2 l3 g 

%%
%Link 1
x1 = 0;
y1 = 0;
z1 = l1;
x1d = diff(x1,t); y1d = diff(y1,t); z1d = diff(z1,t);
v1_sq = x1d^2 + y1d^2 + z1d^2;
    %Kinetic energy for link1
k1 = (1/2)*m1*(v1_sq)
    %Potential energy for link1
p1 = m1*g*z1

%%
%Link 2
x2 = l2*cos(t1)*cos(t2);
y2 = l2*sin(t1)*cos(t2);
z2 = z1+l2*sin(t2);
x2d = diff(x2,t);
y2d = diff(y2,t);
z2d = diff(z2,t);
v2_sq = x2d^2 + y2d^2 + z2d^2;
    %Kinetic energy for link2
k2 = (1/2)*m2*v2_sq
    %Potential energy for link2
p2 = m2*g*z2

%%
%Link 3
x3 = l2*cos(t1)*cos(t2) + l3*cos(t1)*cos(t2+t3);
y3 = l2*sin(t1)*cos(t2) + l3*sin(t1)*cos(t2+t3);
z3 = l3*sin(t2+t3) + l2*sin(t2) +l1;
x3d = diff(x3,t);
y3d = diff(y3,t);
z3d = diff(z3,t);
v3_sq = x3d^2 + y3d^2 + z3d^2;
    %Kinetic energy for link3
k3 = (1/2)*m3*v3_sq
    %Potential energy for link3
p3 = m3*g*z3

%%
%Total kinetic energy
k = k1+k2+k3;
k = simplify(k)

%%
%Total potential energy
p=p1+p2+p3; 
p = simplify(p)

%%
%The lagrangian (L)
L=k-p

%%
%To find torque vector
syms t1d t2d t3d th1 th2 th3

%Substituting terms for allowing differentiation
L = subs(L,[diff(t1(t),t), diff(t2(t),t), diff(t3(t),t)],[t1d, t2d t3d]);
L = subs(L,[t1(t), t2(t), t3(t)],[th1, th2, th3]);

%Partial derivatives for tau1
partial_wrt_t1 = diff(L,th1);
partial_wrt_t1d = diff(L,t1d);

%Partial derivatives for tau2
partial_wrt_t2 = diff(L,th2);
partial_wrt_t2d = diff(L,t2d);

%Partial derivatives for tau3
partial_wrt_t3 = diff(L,th3);
partial_wrt_t3d = diff(L,t3d);

%Reversing substitution to go back to original form 
partial_wrt_t1 = subs(partial_wrt_t1,[th1, th2, th3],[t1(t) t2(t) t3(t)]);
partial_wrt_t1 = subs(partial_wrt_t1,[t1d, t2d, t3d],[diff(t1(t),t),...
    diff(t2(t),t), diff(t3(t),t)]);
partial_wrt_t1d = subs(partial_wrt_t1d,[th1, th2, th3],[t1(t) t2(t) t3(t)]);
partial_wrt_t1d = subs(partial_wrt_t1d,[t1d, t2d, t3d],[diff(t1(t),t),...
    diff(t2(t),t), diff(t3(t),t)]);

partial_wrt_t2 = subs(partial_wrt_t2,[th1, th2, th3],[t1(t) t2(t) t3(t)]);
partial_wrt_t2 = subs(partial_wrt_t2,[t1d, t2d, t3d],[diff(t1(t),t),...
    diff(t2(t),t), diff(t3(t),t)]);
partial_wrt_t2d = subs(partial_wrt_t2d,[th1, th2, th3],[t1(t) t2(t) t3(t)]);
partial_wrt_t2d = subs(partial_wrt_t2d,[t1d, t2d, t3d],[diff(t1(t),t),...
    diff(t2(t),t), diff(t3(t),t)]);

partial_wrt_t3 = subs(partial_wrt_t3,[th1, th2, th3],[t1(t) t2(t) t3(t)]);
partial_wrt_t3 = subs(partial_wrt_t3,[t1d, t2d, t3d],[diff(t1(t),t),...
    diff(t2(t),t), diff(t3(t),t)]);
partial_wrt_t3d = subs(partial_wrt_t3d,[th1, th2, th3],[t1(t) t2(t) t3(t)]);
partial_wrt_t3d = subs(partial_wrt_t3d,[t1d, t2d, t3d],[diff(t1(t),t),...
    diff(t2(t),t), diff(t3(t),t)]);

%Time derivative for tau1
time_derivative_tau1 = diff(partial_wrt_t1d,t);
%Time derivative for tau2
time_derivative_tau2 = diff(partial_wrt_t2d,t);
%Time derivative for tau3
time_derivative_tau3 = diff(partial_wrt_t3d,t);

%%
%Final equation for tau1
tau1 = time_derivative_tau1 - partial_wrt_t1;
tau1 = simplify(tau1)
%Final equation for tau2
tau2 = time_derivative_tau2 - partial_wrt_t2;
tau2 = simplify(tau2)
%Final equation for tau3
tau3 = time_derivative_tau3 - partial_wrt_t3;
tau3 = simplify(tau3)

%Torque vector
Tau = [tau1; tau2; tau3]

%%
%Performing substitutions again, to find inertia, gravity, and coriolis matrices.

syms t1dd t2dd t3dd
Tau = subs(Tau,[diff(t1(t),t,t), diff(t2(t),t,t), diff(t3(t),t,t)],...
    [t1dd, t2dd, t3dd]);
Tau = subs(Tau,[diff(t1(t),t), diff(t2(t),t), diff(t3(t),t)],...
    [t1d, t2d, t3d]);
disp('Your new tau')
disp(Tau)

%%
%Inertia matrix
IM = equationsToMatrix(Tau,[t1dd t2dd t3dd])

%%
%Gravity matrix
GM = equationsToMatrix(Tau,g)

%%
%Coriolis matrix
CM = simplify(expand(Tau - IM*[t1dd;t2dd;t3dd] - GM*g))
