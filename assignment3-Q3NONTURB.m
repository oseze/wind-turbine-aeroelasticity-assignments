%% Clean up
clear; close all; clc;

%% Loading

load data.mat
load bemq1NONTURB.mat
load Q1_ans.mat
load blade_table.mat
load airfoildata.mat
load u.mat


%% Initialization

R=89.17;                % rotor radius[m]
nB=1;                    % number of blades[]
P_rated=10.64*10^6;     % rated power[W]
V0=8;                   % wind speed [m/s]
rho=1.225;              % air density[kg/m^3]
thetacone=(0*pi)/180;   % [rad]
thetayaw=(0*pi)/180;    % [rad]
thetatilt=(0*pi)/180;   % [rad]
A=pi*R^2;               % rotor area [m2]
M_nac=446000;           % mass of nacelle [kg]
k1=1.7*10^6;            % spring stiffness[N/m]
thetapitch=-3.34*pi/180;% pitch angle [rad]
omega1f=3.93;           % first flapwise[rad/s]
omega1e=6.10;           % first edgewise [rad/s]
omega2f=11.28;          % second flapwise [rad/s]
H=119;                  % height of tower[m]
Ls=7.1;                 % shaft length[m]
Wy0=0; Wz0=0;       % first value for induced wind in each direction
L1=1637.6;              % length in the flow direction of the turbulence field [m]
L2=204;                 % length in the horitzaontal direction of the turbulence field [m]
L3=204;                 % length in the vertical direction of the turbulence field [m]
N1=2048;                % grid points in direction 1
N2=256;                 % grid points in direction 2
N3=256;                 % grid points in direction 3
delt=0.1;               % time step [s]
n=1:1:N1;               % elements in time [s]
t=n*delt;               % time [s]
nt=numel(t);
nr=numel(r);
np=N3;
DOF=10;                 % we are using all the system now, 10 DOF
omega=6.423*pi/30;

%% Preallocation

PositionMatrix=zeros(DOF,nt);
VelocityMatrix=zeros(DOF,nt);
AccelerationMatrix=zeros(DOF,nt);
GF=zeros(DOF,nt);

%% BEM Preallocating

a23_b=zeros(3);
a34=zeros(3);
xb=zeros(1,nt);
yb=zeros(1,nt);
zb=zeros(1,nt);
V1=zeros(3,1);
normVrel=zeros(1,nt);
L=zeros(nt,nr);
D=zeros(nt,nr);
p_n=zeros(nr,nt);
p_t=zeros(nr,nt);
Wy=zeros(nt+1,nr);
Wz=zeros(nt+1,nr);
V4=zeros(3,1);
Vrely=zeros(nt,nr);
Vrelz=zeros(nt,nr);
phi=zeros(1,nr);
alpha=zeros(1,nr);
xbinary=zeros(1,np);
ybinary=zeros(1,np);
Output=struct('xyzr',[],'T',[],'P',[],'MR',[]);
xyzr=zeros(3,nr);


%% Transformation to the different coordinate systems

% tranformation matrix from system 1 to 2
a12=zeros(3); a12(1,1)=cos(thetatilt); a12(1,2)=0; a12(1,3)=-sin(thetatilt);
a12(2,1)=sin(thetatilt)*sin(thetayaw); a12(2,2)=cos(thetayaw); a12(2,3)=sin(thetayaw)*cos(thetatilt);
a12(3,1)=sin(thetatilt)*cos(thetayaw); a12(3,2)=-sin(thetayaw); a12(3,3)=cos(thetatilt)*cos(thetayaw);
a21=a12';
ib=1;

%%
for it=1:nt
    
    AccelerationMatrix(:,it)=((GF(:,it)'-(PositionMatrix(:,it)'*KMatrix))/MassMatrix)';
    A=(delt*AccelerationMatrix(:,it))/2;
    b=(delt*(VelocityMatrix(:,it)+0.5*A)/2);
    
    PositionNew=PositionMatrix(:,it)+b;
    VelocityNew=VelocityMatrix(:,it)+A;
    
    
    %now we will run the new velocity triangle with velocityNew used
    %for it=1:nt % iteration over time
    
    thetablade=[omega*t(it) omega*t(it)+(2*pi)/3 omega*t(it)+(4*pi)/3];
    for ir=1:nr-1 % iteration across blade element
%         if ir==17
%             iiii=1;
%         end
        
        V1y=V0*sin(thetayaw);
        V1x=0;
        V1z=V0*cos(thetayaw);
        % tranformation matrix from system 2 to 3
        a23_b(1,:)=[cos(thetablade(ib)),sin(thetablade(ib)),0];
        a23_b(2,:)=[-sin(thetablade(ib)),cos(thetablade(ib)),0];
        a23_b(3,:)=[0,0,1];
        
        % tranformation matrix from system 3 to 4
        a34(1,:)=[cos(thetacone),0,-sin(thetacone)];
        a34(2,:)=[0,1,0];
        a34(3,:)=[sin(thetacone),0,cos(thetacone)];
        
        % tranformation matrix from system 1 to 4 and 4 to 1
        a14_b=(a34*a23_b)*a12;   a41_b=a14_b';
        
        % posistion calculations of the elements
        rt=[H,0,0]; rs=a21*[0,0,-Ls]';
        rb=a41_b*[r(ir),0,0]'; rr=rt+rs+rb;
        
        % position in blade system
        xb=rt(1)+rs(1)+rb(1);
        yb=rt(2)+rs(2)+rb(2);
        zb=rt(3)+rs(3)+rb(3);
        
        delta1=L1/(N1-1);
        delta2=L2/(N2-1);
        delta3=L3/(N3-1);
        
        
        % incoming wind speed in system 1 and 4
        V1=[V1x,V1y,V1z]'; V4=a14_b*V1;
        
        % the velocity in each direction for each blade element in
        % system 4
        xyzr(:,ir)=V4;
        
        % define relative velocity and axial induction factor also in
        %the case when time is 1.
               
        if it==1
            a=abs(Wz0/V0);
            Vrely(it,ir)=V4(2)+Wy0-omega*r(ir)*cos(thetacone)-VelocityNew(2)*u1f_y(ir);
            Vrelz(it,ir)=V4(3)+Wz0-VelocityNew(2)*u1f_z(ir)-VelocityNew(1);
        else
            a=abs(Wz(it,ir)/V0);
            Vrely(it,ir)=V4(2)+Wy(it,ir)-omega*r(ir)*cos(thetacone)-VelocityNew(2)*u1f_y(ir);
            Vrelz(it,ir)=V4(3)+Wz(it,ir)-VelocityNew(2)*u1f_z(ir)-VelocityNew(1);
        end
        
        % flow angle and angle of attack calculation
        phi(ir)=atan(Vrelz(it,ir)/(-Vrely(it,ir)));
        alpha(ir)=(phi(ir)*180/pi)-(-beta(ir)+thetapitch*180/pi);
        
        %interpolation to get drag and lift coefficent
        thick=[100, 60, 48, 36, 30.1, 24.1];
        cl1=interp1(W3_100(:,1),W3_100(:,2),alpha(ir));
        cl2=interp1(W3_60(:,1),W3_60(:,2),alpha(ir));
        cl3=interp1(W3_48(:,1),W3_48(:,2),alpha(ir));
        cl4=interp1(W3_36(:,1),W3_36(:,2),alpha(ir));
        cl5=interp1(W3_30(:,1),W3_30(:,2),alpha(ir));
        cl6=interp1(W3_24(:,1),W3_24(:,2),alpha(ir));
        cd1=interp1(W3_100(:,1),W3_100(:,3),alpha(ir));
        cd2=interp1(W3_60(:,1),W3_60(:,3),alpha(ir));
        cd3=interp1(W3_48(:,1),W3_48(:,3),alpha(ir));
        cd4=interp1(W3_36(:,1),W3_36(:,3),alpha(ir));
        cd5=interp1(W3_30(:,1),W3_30(:,3),alpha(ir));
        cd6=interp1(W3_24(:,1),W3_24(:,3),alpha(ir));
        clvec=[cl1 cl2 cl3 cl4 cl5 cl6];
        Cl=interp1(thick,clvec,tc(ir));
        cdvec=[cd1 cd2 cd3 cd4 cd5 cd6];
        Cd=interp1(thick,cdvec,tc(ir));
        
        % norm of relative velocity
        normVrel(it,ir)=sqrt((Vrely(it,ir))^2+(Vrelz(it,ir))^2);
        
        % computation of lift, drag, tangential load and normal load
        L(ir)=0.5*rho*normVrel(it,ir)^2*c(ir)*Cl;
        D(ir)=0.5*rho*normVrel(it,ir)^2*c(ir)*Cd;
        p_n(ir,it)=L(ir)*cos(phi(ir))+D(ir)*sin(phi(ir));
        p_t(ir,it)=L(ir)*sin(phi(ir))-D(ir)*cos(phi(ir));
        
        % prandtl's tip loss correction
        F=(2/pi)*acos(exp((-nB*(R-r(ir)))/(2*r(ir)*sin(abs(phi(ir))))));
        
        % glauert correction, defined based on a-value
        if a <=(1/3)
            fg=1;
        elseif a > (1/3)
            fg=0.25*(5-3*a);
        end
        
        if it==1
            norm=sqrt((V4(2))^2+(V4(3)+fg*Wz0)^2);
        else
            norm=sqrt((V4(2))^2+(V4(3)+fg*Wz(it,ir))^2);
        end
        
        % induced wind calculation
        Wz(it+1,ir)=-(nB*L(ir)*cos(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
        Wy(it+1,ir)=-(nB*L(ir)*sin(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
    end
    
    % the normal and tangetial loads are defined to be 0 at the tip
    p_n(nr,it)=0; p_t(nr,it)=0;
    
    % computing torque, power and thrust
    MR=trapz(r,p_t(:,it).*r);
    P=omega*MR;
    T=trapz(r,p_n(:,it)');
    
    % now BEM is over, so we will continue solving the EOM
    %update the FG matrix as the loads had changed
    
    GF1=trapz(r,p_t(:,it).*u1f_y)+trapz(r,p_n(:,it).*u1f_z);
    GF2=trapz(r,p_t(:,it).*u1e_y)+trapz(r,p_n(:,it).*u1e_z);
    GF3=trapz(r,p_t(:,it).*u2f_y)+trapz(r,p_n(:,it).*u2f_z);
    GF(:,it)=[T*3, repmat([GF1 GF2 GF3],1,3)];
    
    
    AccelerationMatrix(:,it)=((GF(:,it)'-PositionNew'*KMatrix)/MassMatrix)';
    B=delt*0.5*AccelerationMatrix(:,it);
    
    PositionNew=PositionMatrix(:,it)+b;
    VelocityNew=VelocityMatrix(:,it)+B;
    
    % Again, new BEM to fins new loads for new velocity
    
    for ir=1:nr-1 % iteration across blade element
        
        V1y=V0*sin(thetayaw);
        V1x=0;
        V1z=V0*cos(thetayaw);
        % tranformation matrix from system 2 to 3
        a23_b(1,:)=[cos(thetablade(ib)),sin(thetablade(ib)),0];
        a23_b(2,:)=[-sin(thetablade(ib)),cos(thetablade(ib)),0];
        a23_b(3,:)=[0,0,1];
        
        % tranformation matrix from system 3 to 4
        a34(1,:)=[cos(thetacone),0,-sin(thetacone)];
        a34(2,:)=[0,1,0];
        a34(3,:)=[sin(thetacone),0,cos(thetacone)];
        
        % tranformation matrix from system 1 to 4 and 4 to 1
        a14_b=(a34*a23_b)*a12;   a41_b=a14_b';
        
        % posistion calculations of the elements
        rt=[H,0,0]; rs=a21*[0,0,-Ls]';
        rb=a41_b*[r(ir),0,0]'; rr=rt+rs+rb;
        
        % position in blade system
        xb=rt(1)+rs(1)+rb(1);
        yb=rt(2)+rs(2)+rb(2);
        zb=rt(3)+rs(3)+rb(3);
        
        delta1=L1/(N1-1);
        delta2=L2/(N2-1);
        delta3=L3/(N3-1);
        
                
        % incoming wind speed in system 1 and 4
        V1=[V1x,V1y,V1z]'; V4=a14_b*V1;
        
        % the velocity in each direction for each blade element in
        % system 4
        xyzr(:,ir)=V4;
        
        % define relative velocity and axial induction factor also in
        %the case when time is 1.
        
        if it==1
            a=abs(Wz0/V0);
            Vrely(it,ir)=V4(2)+Wy0-omega*r(ir)*cos(thetacone)-VelocityNew(2)*u1f_y(ir);
            Vrelz(it,ir)=V4(3)+Wz0-VelocityNew(2)*u1f_z(ir)-VelocityNew(1);
        else
            a=abs(Wz(it,ir)/V0);
            Vrely(it,ir)=V4(2)+Wy(it,ir)-omega*r(ir)*cos(thetacone)-VelocityNew(2)*u1f_y(ir);
            Vrelz(it,ir)=V4(3)+Wz(it,ir)-VelocityNew(2)*u1f_z(ir)-VelocityNew(1);
        end
        
        % flow angle and angle of attack calculation
        phi(ir)=atan(Vrelz(it,ir)/(-Vrely(it,ir)));
        alpha(ir)=(phi(ir)*180/pi)-(-beta(ir)+thetapitch*180/pi);
        
        %interpolation to get drag and lift coefficent
        thick=[100, 60, 48, 36, 30.1, 24.1];
        cl1=interp1(W3_100(:,1),W3_100(:,2),alpha(ir));
        cl2=interp1(W3_60(:,1),W3_60(:,2),alpha(ir));
        cl3=interp1(W3_48(:,1),W3_48(:,2),alpha(ir));
        cl4=interp1(W3_36(:,1),W3_36(:,2),alpha(ir));
        cl5=interp1(W3_30(:,1),W3_30(:,2),alpha(ir));
        cl6=interp1(W3_24(:,1),W3_24(:,2),alpha(ir));
        cd1=interp1(W3_100(:,1),W3_100(:,3),alpha(ir));
        cd2=interp1(W3_60(:,1),W3_60(:,3),alpha(ir));
        cd3=interp1(W3_48(:,1),W3_48(:,3),alpha(ir));
        cd4=interp1(W3_36(:,1),W3_36(:,3),alpha(ir));
        cd5=interp1(W3_30(:,1),W3_30(:,3),alpha(ir));
        cd6=interp1(W3_24(:,1),W3_24(:,3),alpha(ir));
        clvec=[cl1 cl2 cl3 cl4 cl5 cl6];
        Cl=interp1(thick,clvec,tc(ir));
        cdvec=[cd1 cd2 cd3 cd4 cd5 cd6];
        Cd=interp1(thick,cdvec,tc(ir));
        
        % norm of relative velocity
        normVrel(it,ir)=sqrt((Vrely(it,ir))^2+(Vrelz(it,ir))^2);
        
        % computation of lift, drag, tangential load and normal load
        L(ir)=0.5*rho*normVrel(it,ir)^2*c(ir)*Cl;
        D(ir)=0.5*rho*normVrel(it,ir)^2*c(ir)*Cd;
        p_n(ir,it)=L(ir)*cos(phi(ir))+D(ir)*sin(phi(ir));
        p_t(ir,it)=L(ir)*sin(phi(ir))-D(ir)*cos(phi(ir));
        
        %     % prandtl's tip loss correction
        %     F=(2/pi)*acos(exp((-nB*(R-r(ir)))/(2*r(ir)*sin(phi(ir)))));
        %
        %     % glauert correction, defined based on a-value
        %     if a <=(1/3)
        %         fg=1;
        %     elseif a > (1/3)
        %         fg=0.25*(5-3*a);
        %     end
        %
        %     if it==1
        %         norm=sqrt((V4(2))^2+(V4(3)+fg*Wz0)^2);
        %     else
        %         norm=sqrt((V4(2))^2+(V4(3)+fg*Wz(it,ir))^2);
        %     end
        
        % induced wind calculation
        %         Wz(it+1,ir)=-(nB*L(ir)*cos(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
        %         Wy(it+1,ir)=-(nB*L(ir)*sin(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
        
    end
    % the normal and tangetial loads are defined to be 0 at the tip
    p_n(nr,it)=0; p_t(nr,it)=0;
    
    % computing torque, power and thrust
    MR=trapz(r,p_t(:,it).*r);
    P=omega*MR;
    T=trapz(r,p_n(:,it)');
    
    %update the FG matrix as the loads had changed
    
    GF1=trapz(r,p_t(:,it).*u1f_y)+trapz(r,p_n(:,it).*u1f_z);
    GF2=trapz(r,p_t(:,it).*u1e_y)+trapz(r,p_n(:,it).*u1e_z);
    GF3=trapz(r,p_t(:,it).*u2f_y)+trapz(r,p_n(:,it).*u2f_z);
    GF(:,it)=[T*3, repmat([GF1 GF2 GF3],1,3)];
    
    %now we solve EOM again
    AccelerationMatrix(:,it)=((GF(:,it)'-PositionNew'*KMatrix)/MassMatrix)';
    C=delt*0.5*AccelerationMatrix(:,it);
    d=delt*(VelocityMatrix(:,it)+C);
    
    PositionNew=PositionMatrix(:,it)+d;
    VelocityNew=VelocityMatrix(:,it)+2*C;
    
    %BEM again
    
    for ir=1:nr-1 % iteration across blade element
        
        V1y=V0*sin(thetayaw);
        V1x=0;
        V1z=V0*cos(thetayaw);
        % tranformation matrix from system 2 to 3
        a23_b(1,:)=[cos(thetablade(ib)),sin(thetablade(ib)),0];
        a23_b(2,:)=[-sin(thetablade(ib)),cos(thetablade(ib)),0];
        a23_b(3,:)=[0,0,1];
        
        % tranformation matrix from system 3 to 4
        a34(1,:)=[cos(thetacone),0,-sin(thetacone)];
        a34(2,:)=[0,1,0];
        a34(3,:)=[sin(thetacone),0,cos(thetacone)];
        
        % tranformation matrix from system 1 to 4 and 4 to 1
        a14_b=(a34*a23_b)*a12;   a41_b=a14_b';
        
        % posistion calculations of the elements
        rt=[H,0,0]; rs=a21*[0,0,-Ls]';
        rb=a41_b*[r(ir),0,0]'; rr=rt+rs+rb;
        
        % position in blade system
        xb=rt(1)+rs(1)+rb(1);
        yb=rt(2)+rs(2)+rb(2);
        zb=rt(3)+rs(3)+rb(3);
        
        delta1=L1/(N1-1);
        delta2=L2/(N2-1);
        delta3=L3/(N3-1);
        
       
        % incoming wind speed in system 1 and 4
        V1=[V1x,V1y,V1z]'; V4=a14_b*V1;
        
        % the velocity in each direction for each blade element in
        % system 4
        xyzr(:,ir)=V4;
        
        % define relative velocity and axial induction factor also in
        %the case when time is 1.
        
        if it==1
            a=abs(Wz0/V0);
            Vrely(it,ir)=V4(2)+Wy0-omega*r(ir)*cos(thetacone)-VelocityNew(2)*u1f_y(ir);
            Vrelz(it,ir)=V4(3)+Wz0-VelocityNew(2)*u1f_z(ir)-VelocityNew(1);
        else
            a=abs(Wz(it,ir)/V0);
            Vrely(it,ir)=V4(2)+Wy(it,ir)-omega*r(ir)*cos(thetacone)-VelocityNew(2)*u1f_y(ir);
            Vrelz(it,ir)=V4(3)+Wz(it,ir)-VelocityNew(2)*u1f_z(ir)-VelocityNew(1);
        end
        
        % flow angle and angle of attack calculation
        phi(ir)=atan(Vrelz(it,ir)/(-Vrely(it,ir)));
        alpha(ir)=(phi(ir)*180/pi)-(-beta(ir)+thetapitch*180/pi);
        
        %interpolation to get drag and lift coefficent
        thick=[100, 60, 48, 36, 30.1, 24.1];
        cl1=interp1(W3_100(:,1),W3_100(:,2),alpha(ir));
        cl2=interp1(W3_60(:,1),W3_60(:,2),alpha(ir));
        cl3=interp1(W3_48(:,1),W3_48(:,2),alpha(ir));
        cl4=interp1(W3_36(:,1),W3_36(:,2),alpha(ir));
        cl5=interp1(W3_30(:,1),W3_30(:,2),alpha(ir));
        cl6=interp1(W3_24(:,1),W3_24(:,2),alpha(ir));
        cd1=interp1(W3_100(:,1),W3_100(:,3),alpha(ir));
        cd2=interp1(W3_60(:,1),W3_60(:,3),alpha(ir));
        cd3=interp1(W3_48(:,1),W3_48(:,3),alpha(ir));
        cd4=interp1(W3_36(:,1),W3_36(:,3),alpha(ir));
        cd5=interp1(W3_30(:,1),W3_30(:,3),alpha(ir));
        cd6=interp1(W3_24(:,1),W3_24(:,3),alpha(ir));
        clvec=[cl1 cl2 cl3 cl4 cl5 cl6];
        Cl=interp1(thick,clvec,tc(ir));
        cdvec=[cd1 cd2 cd3 cd4 cd5 cd6];
        Cd=interp1(thick,cdvec,tc(ir));
        
        % norm of relative velocity
        normVrel(it,ir)=sqrt((Vrely(it,ir))^2+(Vrelz(it,ir))^2);
        
        % computation of lift, drag, tangential load and normal load
        L(ir)=0.5*rho*normVrel(it,ir)^2*c(ir)*Cl;
        D(ir)=0.5*rho*normVrel(it,ir)^2*c(ir)*Cd;
        p_n(ir,it)=L(ir)*cos(phi(ir))+D(ir)*sin(phi(ir));
        p_t(ir,it)=L(ir)*sin(phi(ir))-D(ir)*cos(phi(ir));
        
        %     % prandtl's tip loss correction
        %     F=(2/pi)*acos(exp((-nB*(R-r(ir)))/(2*r(ir)*sin(phi(ir)))));
        %
        %     % glauert correction, defined based on a-value
        %     if a <=(1/3)
        %         fg=1;
        %     elseif a > (1/3)
        %         fg=0.25*(5-3*a);
        %     end
        %
        %
        %         norm=sqrt((V4(2))^2+(V4(3)+fg*Wz(it,ir))^2);
        
        
        % induced wind calculation
        %         Wz(it+1,ir)=-(nB*L(ir)*cos(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
        %         Wy(it+1,ir)=-(nB*L(ir)*sin(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
        
    end
    % the normal and tangetial loads are defined to be 0 at the tip
    p_n(nr,it)=0; p_t(nr,it)=0;
    
    % computing torque, power and thrust
    MR=trapz(r,p_t(:,it).*r);
    P=omega*MR;
    T=trapz(r,p_n(:,it)');
    
    %update the FG matrix as the loads had changed
    
    GF1=trapz(r,p_t(:,it).*u1f_y)+trapz(r,p_n(:,it).*u1f_z);
    GF2=trapz(r,p_t(:,it).*u1e_y)+trapz(r,p_n(:,it).*u1e_z);
    GF3=trapz(r,p_t(:,it).*u2f_y)+trapz(r,p_n(:,it).*u2f_z);
    GF(:,it)=[T*3, repmat([GF1 GF2 GF3],1,3)];
    
    %final calculation
    AccelerationMatrix(:,it)=((GF(:,it)'-PositionNew'*KMatrix)/MassMatrix)';
    D=delt*0.5*AccelerationMatrix(:,it);
    
    PositionMatrix(:,it+1)=PositionMatrix(:,it)+delt*(VelocityMatrix(:,it)+(1/3)*(A+B+C));
    VelocityMatrix(:,it+1)=VelocityMatrix(:,it)+(1/3)*(A+2*B+2*C+D);
    
end

PositionMatrix(:,end)=[];
VelocityMatrix(:,end)=[];

%% Deflection in Flapwise and Edgewise at the Tip

DeflectionTipFlap=PositionMatrix(2,:)*u1f_z(end)+PositionMatrix(3,:)*u1e_z(end)+PositionMatrix(4,:)*u2f_z(end);

DeflectionTipEdge=PositionMatrix(2,:)*u1f_y(end)+PositionMatrix(3,:)*u1e_y(end)+PositionMatrix(4,:)*u2f_y(end);

%% Bending moment in Flapwise and Edgewise at the root


inertia_tower=M_nac*AccelerationMatrix(1,:);

inertiaload_Flap=m.*(AccelerationMatrix(2,:).*u1f_z+AccelerationMatrix(3,:).*u1e_z+AccelerationMatrix(4,:).*u2f_z);
inertiaload_Edge=m.*(AccelerationMatrix(2,:).*u1f_y+AccelerationMatrix(3,:).*u1e_y+AccelerationMatrix(4,:).*u2f_y);

bending_moment_Flap=trapz(r,(p_n-inertiaload_Flap).*r');
bending_moment_Edge=trapz(r,(p_t-inertiaload_Edge).*r');% tower inertia doesn't conribute in edge direction

%%

bending_moment_Flapwo=trapz(r,(p_n.*r'));
bending_moment_Edgewo=trapz(r,(p_t.*r'));

%%
impact_Flap=((bending_moment_Flapwo-bending_moment_Flap)/bending_moment_Flapwo)*100;
impact_Edge=((bending_moment_Edgewo-bending_moment_Edge)/bending_moment_Edgewo)*100;



%%

save('Q3nonturb_answ.mat','t','PositionMatrix','AccelerationMatrix','VelocityMatrix','bending_moment_Flap','bending_moment_Edge','DeflectionTipFlap','DeflectionTipEdge','bending_moment_Flapwo','bending_moment_Edgewo','impact_Flap','impact_Edge')