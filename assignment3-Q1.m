%% Clean up
clear; close all; clc;

%% Initialization
load data.mat
load bemq1.mat

%%

R=89.17;                % rotor radius[m]
B=3;                    % number of blades[]
P_rated=10.64*10^6;     % rated power[W]
V0=15;                  % wind speed [m/s]
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
r=modeshapes(:,1);
%%
rtemp=[r' R];
dr=diff(rtemp);

%% Mass Matrix

GM1=trapz(r',u1f_y.*m.*u1f_y)+trapz(r',u1f_z.*m.*u1f_z); % gen. mass for DOF 1
GM2=trapz(r',u1e_y.*m.*u1e_y)+trapz(r',u1e_z.*m.*u1e_z); % gen. mass for DOF 2
GM3=trapz(r',u2f_y.*m.*u2f_y)+trapz(r',u2f_z.*m.*u2f_z); % gen. mass for DOF 3
GM=repmat([GM1 GM2 GM3],1,3); % repeating gen. mass 1,2 and 3, thrice 
DIAG_GM=diag(GM);             % diagonal matrix with rep. gen. mass in diagonal
m_1=trapz(r',m.*u1f_z);       
m_2=trapz(r',m.*u1e_z);
m_3=trapz(r',m.*u2f_z);
m1line=repmat([m_1, m_2, m_3],1,3);
Mtotal=3*sum(dr'.*m)+M_nac;   % total mass for the three blades and nacelle
MassTemp=[m1line', DIAG_GM];
MassMatrix=[[Mtotal m1line]; MassTemp];

%% Stiffness Matrix

eigenfreq=repmat([omega1f^2, omega1e^2, omega2f^2],1,3); % eigenfrequencies squared repeated thrice
K_F=eigenfreq.*GM;                                       % force vector
KMatrix=diag([k1, K_F]);                                 % stiffness matrix

%% Generalized Force Vector NOTE: we assumed here pn and pt just for last time series! 
nt=2048;

for it=1:nt
    
GF1=trapz(r,p_t(:,it).*u1f_y)+trapz(r,p_n(:,it).*u1f_z);
GF2=trapz(r,p_t(:,it).*u1e_y)+trapz(r,p_n(:,it).*u1e_z);
GF3=trapz(r,p_t(:,it).*u2f_y)+trapz(r,p_n(:,it).*u2f_z);
GF(:,it)=[T*3, repmat([GF1 GF2 GF3],1,3)]; % repeating force mass 1,2 and 3, thrice
%GF=[sum([GF1 GF2 GF3]), repmat([GF1 GF2 GF3],1,3)]; % repeating force mass 1,2 and 3, thrice 
end

%% Estimation of the Tower For-aft

omega_tower_root=sqrt(k1/M_tower);   % root (simple) formula
e=eig(KMatrix,MassMatrix);           % eigen values for matrix K and M
eigenvalues=sqrt(e);                  
omega_tower_eig=sqrt(e(1));          


%% Saving 

save('Q1_ans.mat','MassMatrix','KMatrix','GF','omega_tower_root','omega_tower_eig');

