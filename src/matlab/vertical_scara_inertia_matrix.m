close all
clear all
clc

%% INITIALISATION
syms q1 q2 q3
syms m1 m2 m3
syms lc1 lc2
syms l1 l2
syms rho1 rho2 rho3 l1x l1y l1z l2x l2y l2z l3x l3y l3z

%% "PARTIAL" JACOBIANS
Jvc1 = [
  -lc1*sin(q1) , 0 , 0 ;
   lc1*cos(q1) , 0 , 0 ;
   0           , 0 , 0
];
Jwc1 = [
  0 , 0 , 0 ;
  0 , 0 , 0 ;
  1 , 0 , 0
];
Jvc2 = [
  -l1*sin(q1) - lc2*sin(q1+q2) , -lc2*sin(q1+q2) , 0 ;
   l1*cos(q1) + lc2*cos(q1+q2) ,  lc2*cos(q1+q2) , 0 ;
   0                           ,  0              , 0
];
Jwc2 = [
  0 , 0 , 0 ;
  0 , 0 , 0 ;
  1 , 1 , 0
];
Jvc3 = [
  -l1*sin(q1) - l2*sin(q1+q2) , -l2*sin(q1+q2) ,  0 ;
   l1*cos(q1) + l2*cos(q1+q2) ,  l2*cos(q1+q2) ,  0 ;
   0                          ,  0             , -1
];
Jwc3 = Jwc2;

%% INERTIA TENSORS
Ic1 = diag([
  rho1*l1x*l1y*l1z*(l1y*l1y+l1z*l1z)/12 ,
  rho1*l1x*l1y*l1z*(l1x*l1x+l1z*l1z)/12 ,
  rho1*l1x*l1y*l1z*(l1x*l1x+l1y*l1y)/12
]);
Ic2 = diag([
  rho2*l2x*l2y*l2z*(l2y*l2y+l2z*l2z)/12 ,
  rho2*l2x*l2y*l2z*(l2x*l2x+l2z*l2z)/12 ,
  rho2*l2x*l2y*l2z*(l2x*l2x+l2y*l2y)/12
]);
Ic3 = diag([
  rho3*l3x*l3y*l3z*(l3y*l3y+l3z*l3z)/12 ,
  rho3*l3x*l3y*l3z*(l3x*l3x+l3z*l3z)/12 ,
  rho3*l3x*l3y*l3z*(l3x*l3x+l3y*l3y)/12
]);

%% ROTATION INERTIA MATRICES
Rc1 = [
  cos(q1) , -sin(q1) , 0 ;
  sin(q1) ,  cos(q1) , 0 ;
  0       ,  0       , 1
];
Rc2 = [
  cos(q1+q2) ,  sin(q1+q2) ,  0 ;
  sin(q1+q2) , -cos(q1+q2) ,  0 ;
  0          ,  0          , -1
];
Rc3 = Rc2;

%% INERTIA MATRIX
D = ...
  m1*Transpose(Jvc1)*Jvc1 + Transpose(Jwc1)*Rc1*Ic1*Transpose(Rc1)*Jwc1 + ...
  m2*Transpose(Jvc2)*Jvc2 + Transpose(Jwc2)*Rc2*Ic2*Transpose(Rc2)*Jwc2 + ...
  m3*Transpose(Jvc3)*Jvc3 + Transpose(Jwc3)*Rc3*Ic3*Transpose(Rc3)*Jwc3 
