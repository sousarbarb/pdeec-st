close all
clear all
clc

%% INITIALISATION
syms q1 q2 q3
syms m1 m2 m3
syms lc1 lc2
syms l1 l2
syms rho1 rho2 rho3 l1x l1y l1z l2x l2y l2z l3x l3y l3z
syms g

%% POTENTIAL ENERGY
P1 = -g*lc1*sin(q1)*m1;
P2 = -g*l1 *sin(q1)*m2 - g*lc2*sin(q1+q2)*m2;
P3 = -g*l1 *sin(q1)*m3 - g*l2 *sin(q1+q2)*m3;
P  = P1 + P2 + P3;

%% PHI = d P / d qi
Phi1 = diff(P,q1)
Phi2 = diff(P,q2)
Phi3 = diff(P,q3)