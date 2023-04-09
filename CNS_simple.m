% CNS(MANS) encode & decode
% edited by JBR，2023年4月9日
clear;clc;close all;
hhtmqi = [50,50,50]';
xkxk1 = [0,0,1000]'; %行星坐标
xkxk2 = [1000,0,0]';
s1 = [1,2,3]'/(1^2+2^2+3^2)^0.5; %恒星向量
s2 = [2,1,3]'/(2^2+1^2+3^2)^0.5;
s3 = [3,2,1]'/(3^2+2^2+1^2)^0.5;
%% 行星1
A11_ideal = acosd((hhtmqi-xkxk1)'/norm(hhtmqi-xkxk1)*s1);
A12_ideal = acosd((hhtmqi-xkxk1)'/norm(hhtmqi-xkxk1)*s2);
A13_ideal = acosd((hhtmqi-xkxk1)'/norm(hhtmqi-xkxk1)*s3);
A1 = [A11_ideal,A12_ideal,A13_ideal]'+0.01*randn(3,1);
% (l1 l2 l3)(s11,s21,s31;s12,s22,s32;s13,s23,s33) = (cosa(A1),cosd(A2),cosd(A3))
L1_ideal = [cosd(A11_ideal),cosd(A12_ideal),cosd(A13_ideal)]*[s1,s2,s3]^-1; %求得的行星1对航天器的向量（单位向量，理想值）
L1 = [cosd(A1(1)),cosd(A1(2)),cosd(A1(3))]*[s1,s2,s3]^-1; %行星1对航天器的向量（单位向量，带噪声）
%% 行星2
A21_ideal = acosd((hhtmqi-xkxk2)'/norm(hhtmqi-xkxk2)*s1);
A22_ideal = acosd((hhtmqi-xkxk2)'/norm(hhtmqi-xkxk2)*s2);
A23_ideal = acosd((hhtmqi-xkxk2)'/norm(hhtmqi-xkxk2)*s3);
A2 = [A21_ideal,A22_ideal,A23_ideal]'+0.01*randn(3,1);
% (l1 l2 l3)(s11,s21,s31;s12,s22,s32;s13,s23,s33) = (cosa(A1),cosd(A2),cosd(A3))
L2_ideal = [cosd(A21_ideal),cosd(A22_ideal),cosd(A23_ideal)]*[s1,s2,s3]^-1;
L2 = [cosd(A2(1)),cosd(A2(2)),cosd(A2(3))]*[s1,s2,s3]^-1;
%% 求航天器位置
Ac = [L1',-L2']; %vector
rho = (Ac'*Ac)^-1*Ac'*(xkxk2-xkxk1); %the length of L1 & L2(consider noisy)
P1 = xkxk1+rho(1)*L1'
P2 = xkxk2+rho(2)*L2'