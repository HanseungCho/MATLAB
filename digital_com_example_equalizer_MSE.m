clc; clear;
x=[[0.0110 0 0 0 0 0 0];
   [0.0227 0.0110 0 0 0 0 0];
   [-0.1749 0.0227 0.0110 0 0 0 0];
   [1.0000 -0.1749 0.0227 0.0110 0 0 0];
   [0.1617 1.0000 -0.1749 0.0227 0.0110 0 0];
   [-0.0558 0.1617 1.0000 -0.1749 0.0227 0.0110 0];
   [0.0108 -0.0558 0.1617 1.0000 -0.1749 0.0227 0.0110];
   [0 0.0108 -0.0558 0.1617 1.0000 -0.1749 0.0227];
   [0 0 0.0108 -0.0558 0.1617 1.0000 -0.1749];
   [0 0 0 0.0108 -0.0558 0.1617 1.0000];
   [0 0 0 0 0.0108 -0.0558 0.1617];
   [0 0 0 0 0 0.0108 -0.0558];
   [0 0 0 0 0 0 0.0108]];
length=size(x);
z = [0;-0.0428; 0; 1; 0; -0.0071; 0.0345];
Rxx=transpose(x)*x;
Rxz=transpose(x)*z;
c=inv(Rxx)*Rxz;
