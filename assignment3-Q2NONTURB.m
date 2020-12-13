%% Clean up
clear; close all; clc;

%% Initialization
load data.mat
load blade_table.mat
load airfoildata.mat
load bemq1NONTURB.mat

R=89.17;                % rotor radius[m]
B=3;                    % number of blades[]
P_rated=10.64*10^6;     % rated power[W]
V0=8;                  % wind speed [m/s]
rho=1.225;              % air density[kg/m^3]
thetacone=(0*pi)/180;   % [rad]
thetayaw=(0*pi)/180;    % [rad]
thetatilt=(0*pi)/180;   % [rad]
A=pi*R^2;               % rotor area [m2]
M_nac=446000;           % mass of nacelle [kg]
M_tower=628442;         % mass of tower [kg]
k1=1.7*10^6;            % spring stiffness[N/m]
pitch=-3.34;            % pitch angle [deg]
omega1f=3.93;           % first flapwise[rad/s]
omega1e=6.10;           % first edgewise [rad/s]
omega2f=11.28;          % second flapwise [rad/s]
DOF=3;                  % degrees of freedon []
N1=2048;                % grid points in direction 1
deltat=0.1;             % time step [s]
n=1:1:N1;               % elements in time [s]
t=n*deltat;             % time [s]
nt=numel(t); 



%%

GM1=trapz(r',u1f_y.*m.*u1f_y)+trapz(r',u1f_z.*m.*u1f_z); % gen. mass for DOF 1
GM2=trapz(r',u1e_y.*m.*u1e_y)+trapz(r',u1e_z.*m.*u1e_z); % gen. mass for DOF 2
GM3=trapz(r',u2f_y.*m.*u2f_y)+trapz(r',u2f_z.*m.*u2f_z); % gen. mass for DOF 3

GM=[GM1 GM2 GM3];
MassMatrix=diag(GM);

eigsq=[omega1f^2, omega1e^2, omega2f^2];
KMatrix=diag(GM.*eigsq);

delta1=0.03;
delta2=delta1;
delta3=delta1;

eig=[omega1f, omega1e, omega2f];
deltavec=[delta1, delta2, delta3];
DampingMatrix=diag((eig.*deltavec.*GM)/pi);

%NOTE: pt and pn are from last time step
for it=1:nt
    
GF1=trapz(r,p_t(:,it).*u1f_y)+trapz(r,p_n(:,it).*u1f_z);
GF2=trapz(r,p_t(:,it).*u1e_y)+trapz(r,p_n(:,it).*u1e_z);
GF3=trapz(r,p_t(:,it).*u2f_y)+trapz(r,p_n(:,it).*u2f_z);
GF(:,it)=[GF1 GF2 GF3];
end

%% Preallocation 

PositionMatrix=zeros(DOF,nt);
VelocityMatrix=zeros(DOF,nt);
AccelerationMatrix=zeros(DOF,nt);

%%

for it=1:nt
     
   AccelerationMatrix(:,it)=(((GF(:,it)'-VelocityMatrix(:,it)'*DampingMatrix)-(PositionMatrix(:,it)'*KMatrix))/MassMatrix)';
   A=(deltat*AccelerationMatrix(:,it))/2;
   b=(deltat*(VelocityMatrix(:,it)+0.5*A)/2);
   
   PositionNew=PositionMatrix(:,it)+b;
   VelocityNew=VelocityMatrix(:,it)+A;
   
   AccelerationMatrix(:,it)=((GF(:,it)'-VelocityNew'*DampingMatrix-PositionNew'*KMatrix)/MassMatrix)';
   B=deltat*0.5*AccelerationMatrix(:,it);
   
   PositionNew=PositionMatrix(:,it)+b;
   VelocityNew=VelocityMatrix(:,it)+B;
   
   AccelerationMatrix(:,it)=((GF(:,it)'-VelocityNew'*DampingMatrix-PositionNew'*KMatrix)/MassMatrix)';
   C=deltat*0.5*AccelerationMatrix(:,it);
   d=deltat*(VelocityMatrix(:,it)+C);
   
   PositionNew=PositionMatrix(:,it)+d;
   VelocityNew=VelocityMatrix(:,it)+2*C;
   
   AccelerationMatrix(:,it)=((GF(:,it)'-VelocityNew'*DampingMatrix-PositionNew'*KMatrix)/MassMatrix)';
   D=deltat*0.5*AccelerationMatrix(:,it);
   
   PositionMatrix(:,it+1)=PositionMatrix(:,it)+deltat*(VelocityMatrix(:,it)+(1/3)*(A+B+C));
   VelocityMatrix(:,it+1)=VelocityMatrix(:,it)+(1/3)*(A+2*B+2*C+D); 
 
end
%%
PositionMatrix(:,end)=[];
VelocityMatrix(:,end)=[];

%% Deflection in Flapwise and Edgewise at the Tip
% 
% DeflectionTipF1_y=PositionMatrix(1,:)*u1f_y(end);
% DeflectionTipF1_z=PositionMatrix(1,:)*u1f_z(end);
% 
% DeflectionTipE1_y=PositionMatrix(2,:)*u1e_y(end);
% DeflectionTipE1_z=PositionMatrix(2,:)*u1e_z(end);
% 
% DeflectionTipF2_y=PositionMatrix(3,:)*u2f_y(end);
% DeflectionTipF2_z=PositionMatrix(3,:)*u2f_z(end);

DeflectionTipFlap=PositionMatrix(1,:)*u1f_z(end)+PositionMatrix(2,:)*u1e_z(end)+PositionMatrix(3,:)*u2f_z(end);

DeflectionTipEdge=PositionMatrix(1,:)*u1f_y(end)+PositionMatrix(2,:)*u1e_y(end)+PositionMatrix(3,:)*u2f_y(end);

meanF=mean(DeflectionTipFlap);
meanE=mean(DeflectionTipEdge);
%% Bending moment in Flapwise and Edgewise at the root

% 
% inertiaload_F1=m.*AccelerationMatrix(1,:).*u1f_z;
% inertiaload_E1=m.*AccelerationMatrix(2,:).*u1e_y;
% inertiaload_F2=m.*AccelerationMatrix(3,:).*u2f_z;
% 
inertiaload_Flap=m.*(AccelerationMatrix(1,:).*u1f_z+AccelerationMatrix(2,:).*u1e_z+AccelerationMatrix(3,:).*u2f_z);
inertiaload_Edge=m.*(AccelerationMatrix(1,:).*u1f_y+AccelerationMatrix(2,:).*u1e_y+AccelerationMatrix(3,:).*u2f_y);


% 
% bending_moment_F1=trapz(r,(p_n-inertiaload_F1).*r');
% bending_moment_E1=trapz(r,(p_n-inertiaload_E1).*r');
% bending_moment_F2=trapz(r,(p_n-inertiaload_F2).*r');
% 
bending_moment_Flap=trapz(r,(p_n-inertiaload_Flap).*r');
bending_moment_Edge=trapz(r,(p_t-inertiaload_Edge).*r');

% bending_moment_F1wo=trapz(r,(p_n.*r'));
% bending_moment_E1wo=trapz(r,(p_n.*r'));
% bending_moment_F2wo=trapz(r,(p_n.*r'));


bending_moment_Flapwo=trapz(r,(p_n.*r'));
bending_moment_Edgewo=trapz(r,(p_t.*r'));

% impact_F1=(1-bending_moment_F1wo/bending_moment_F1)*100; % the impact of the inertial load is
% impact_F2=(1-bending_moment_E1wo/bending_moment_E1)*100; % the impact of the inertial load is
% impact_F3=(1-bending_moment_F2wo/bending_moment_F2)*100; % the impact of the inertial load is

impact_Flap=((bending_moment_Flapwo-bending_moment_Flap)/bending_moment_Flapwo)*100;
impact_Edge=((bending_moment_Edgewo-bending_moment_Edge)/bending_moment_Edgewo)*100;

%% 
bending_moment_Flap_noturb=bending_moment_Flap;
bending_moment_Edge_noturb=bending_moment_Edge;

DeflectionTipFlap_noturb=DeflectionTipFlap;
DeflectionTipEdge_noturb=DeflectionTipEdge;

%%
meanBMF=mean(bending_moment_Flap);
meanBME=mean(bending_moment_Edge);

meanBMFwo=mean(bending_moment_Flapwo);
meanBMEwo=mean(bending_moment_Edgewo);
%%
save('Q2nonturb_answ.mat','bending_moment_Flap_noturb','PositionMatrix','bending_moment_Edge_noturb','DeflectionTipFlap_noturb','DeflectionTipEdge_noturb','t','AccelerationMatrix','VelocityMatrix')


