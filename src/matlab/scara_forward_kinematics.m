close all
clear all
clc

%% INITIALIZATION
% Joints
syms th1 th2 d3
% Lengths
syms l0 l1 l2

%% FORWARD KINEMATICS
A1 = DH(th1,l0,l1, 0)
A2 = DH(th2, 0,l2,pi)
A3 = DH(  0,d3, 0, 0)

H = A1*A2*A3

%% TEST
theta1 = 0;
theta2 = pi/2;
delta3 = 0;
length0 = 0.3;
length1 = 0.3;
length2 = 0.3;

Htest = double(                                        ...
  subs(H                                             , ...
       [th1    th2    d3     l0      l1      l2]     , ...
       [theta1 theta2 delta3 length0 length1 length2]) ...
)