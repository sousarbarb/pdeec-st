close all
clear all
clc

%% INITIALIZATION
% Joints
syms q1 q2 q3
% Lengths
syms l1 l2

%% FORWARD KINEMATICS
A1 = DH(q1, 0,l1, 0)
A2 = DH(q2, 0,l2,pi)
A3 = DH( 0,q3, 0, 0)

H = A1*A2*A3