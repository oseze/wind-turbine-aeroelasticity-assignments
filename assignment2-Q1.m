%% Clean up
clear; close all; clc;

%% Initialization
load blade_table.mat
load airfoildata.mat
load Cp_opt.mat

R=89.17;                %rotor radius[m]
B=3;                    %number of blades[]
P_rated=10.64*10^6;     %rated power[W]
V0_cutin=4;             %cut-in wind speed [m/s]
V0_cutout=25;           %cut-out wind speed[m/s]
V0=15;                  %wind speed [m/s]
rho=1.225;              %air density[kg/m^3]
thetacone=(0*pi)/180;   %[rad]
thetayaw=(0*pi)/180;    %[rad]
thetatilt=(0*pi)/180;   %[rad]
lambda=8;               %tip speed ratio
A=pi*R^2;               %rotor area [m2]

omega_cutin=lambda*V0_cutin/R;          %cut-in value of omega [rad/s]
omega_rated=1.01;                       %omega rated [rad/s]
omega_cutout=lambda*V0_cutout/R;        %cut-out value of omega [rad/s]

V0_rated=omega_rated*R/lambda;          %rated velocity

K=0.5*rho*A*CP_opt*((R^3)/lambda^3);    %constant value for the generator torque characteristics
I_rot=1.6*10^8;                         %inertia moment of the drivetrain [kgm2]

omega_rated_new=(P_rated/K)^(1/3);      %rated value of omega calculated with P_rated, slightly higher than omega_rated
omega_reference=1.01*omega_rated_new;   %reference value of omega, calculated 1% higher than omega_rated to make sure we are in the right area
%omega_reference is the value we want to stabilize at.


H=119;                                  %height of tower[m]
Ls=7.1;                                 %shaft length[m]
delt=0.1;                               %time step [s]
simtime=60;                             %simulation time [s]
n=1:delt:simtime/delt;                  %element in time []
t=n*delt;                               %time [s]
nt=numel(t); nb=3; nr=numel(r);         %number of elements in: time, blade and blade elements
Wy0=0; Wz0=0;                           %first value for induced wind in each direction


KK=(14*pi)/180;                         % Gain scheduling factor [rad]
Kp=1.5;                                 % Proportional gain [rad/(rad/s)]
Ki=0.64;                                % Integral gain [rad/rad]
thetamaxdot=(8*pi)/180;                 % Maximum pitch velocity [rad/s]
thetapitchmax=(40*pi)/180;              % Maximum pitch angle [rad]
thetapitchmin=0;                        % Minimum pitch angle [rad]

%% Preallocating

Mg=zeros(1,nt);
Mg_nopi=zeros(1,nt);

a23_b=zeros(3);
a34=zeros(3);
xb=zeros(1,nt);
yb=zeros(1,nt);
zb=zeros(1,nt);
V1=zeros(3,1);
normVrel=zeros(1,nt);
L=zeros(nt,nr);
D=zeros(nt,nr);
p_n=zeros(1,nr);
p_t=zeros(1,nr);
Wy=zeros(nt+1,nr);
Wz=zeros(nt+1,nr);
V4=zeros(3,1);
Vrely=zeros(nt,nr);
Vrelz=zeros(nt,nr);

Output=struct('xyzr',[],'T',[],'P',[],'MR',[]);
xyzr=zeros(3,nr);
phi=zeros(1,nr);
alpha=zeros(1,nr);
Power1=zeros(1,nt);
Thrust1=zeros(1,nt);

thetapitch=zeros(1,nt);     % Total pitch angle [rad]
thetapitch_p=zeros(1,nt);   % Proportional part of the pitch angle [rad]
thetapitch_i=zeros(1,nt);   % Integral part of the pitch angle [rad]

totalP=zeros(1,nt);
totalT=zeros(1,nt);
totalMR=zeros(1,nt);

omega=zeros(1,nt);
omega(1)=0.1; % Here, we are setting the first value of omega in time
%% Transformation to the different coordinate systems

% wind speed in coordinate system 1 in each direction
V1z=V0*cos(thetayaw);
V1y=V0*sin(thetayaw);
V1x=0;

% tranformation matrix from system 1 to 2
a12=zeros(3); a12(1,1)=cos(thetatilt); a12(1,2)=0; a12(1,3)=-sin(thetatilt);
a12(2,1)=sin(thetatilt)*sin(thetayaw); a12(2,2)=cos(thetayaw); a12(2,3)=sin(thetayaw)*cos(thetatilt);
a12(3,1)=sin(thetatilt)*cos(thetayaw); a12(3,2)=-sin(thetayaw); a12(3,3)=cos(thetatilt)*cos(thetayaw);
a21=a12';


%% main loop

for it=1:nt % iteration over time
    
    thetablade=[omega(it)*t(it) omega(it)*t(it)+(2*pi)/3 omega(it)*t(it)+(4*pi)/3];
    for ib=1:3 % iteration across blade number
        
        for ir=1:nr-1 % iteration across blade element
            
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
            
            % incoming wind speed in system 1 and 4
            V1=[V1x,V1y,V1z]'; V4=a14_b*V1;
            
            % the velocity in each direction for each blade element in
            % system 4
            xyzr(:,ir)=V4;
            
            % define relative velocity and axial induction factor also in
            %the case when time is 1.
            if it==1
                a=abs(Wz0/V0);
                Vrely(it,ir)=V4(2)+Wy0-omega(it)*r(ir)*cos(thetacone);
                Vrelz(it,ir)=V4(3)+Wz0;
            else
                a=abs(Wz(it,ir)/V0);
                Vrely(it,ir)=V4(2)+Wy(it,ir)-omega(it)*r(ir)*cos(thetacone);
                Vrelz(it,ir)=V4(3)+Wz(it,ir);
            end
            
            % flow angle and angle of attack calculation
            phi(ir)=atan(Vrelz(it,ir)/(-Vrely(it,ir)));
            alpha(ir)=(phi(ir)*180/pi)-(-beta(ir)+thetapitch(it)*180/pi);
            
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
            p_n(ir)=L(ir)*cos(phi(ir))+D(ir)*sin(phi(ir));
            p_t(ir)=L(ir)*sin(phi(ir))-D(ir)*cos(phi(ir));
            
            % prandtl's tip loss correction
            F=(2/pi)*acos(exp((-B*(R-r(ir)))/(2*r(ir)*sin(phi(ir)))));
            
            % glauert correction, defined based on a-value
            if a <=(1/3)
                fg=1;
            elseif a > (1/3)
                fg=0.25*(5-3*a);
            end
            
            %
            if it==1
                norm=sqrt((V4(2))^2+(V4(3)+fg*Wz0)^2);
            else
                norm=sqrt((V4(2))^2+(V4(3)+fg*Wz(it,ir))^2);
            end
            
            % induced wind calculation
            Wz(it+1,ir)=-(B*L(ir)*cos(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
            Wy(it+1,ir)=-(B*L(ir)*sin(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
        end
        
        % the normal and tangetial loads are defined to be 0 at the tip
        p_n(nr)=0; p_t(nr)=0;
        
        % computing torque, power and thrust
        MR=trapz(r,p_t.*r);
        P=omega(it)*MR;
        T=trapz(r,p_n');
        
        % saving torque, power, thrust and velocity matrix onto struct
        Output(it,ib).MR=MR;
        Output(it,ib).P=P;
        Output(it,ib).T=T;
        Output(it,ib).xyzr=xyzr;
        
    end
    
    %Calculation of the total power and thrust
    
    totalP(it)=Output(it,1).P+Output(it,2).P+Output(it,3).P;
    totalT(it)=Output(it,1).T+Output(it,2).T+Output(it,3).T;
    totalMR(it)=Output(it,1).MR+Output(it,2).MR+Output(it,3).MR;
    
    
    %     % Calculation of the generator torque without the PI controller
    %      if omega_nopi(it)<omega_rated
    %
    %         Mg_nopi(it)=K*omega_nopi(it)^(2);
    %     else
    %
    %         Mg_nopi(it)=P_rated/omega_rated;
    %      end
    
    %calculation of generator torque MG and implementation of the PI
    %controller
    
    if omega(it)<omega_rated_new
        Mg(it)=K*omega(it)^(2);
        
    else
        Mg(it)=K*omega_rated_new^2;  %We use this expression instead of P_rated/omega_rated because we want to avoid region 2
    end
    GK=1/(1+thetapitch(it)/KK); % Gain scheduling constant
    thetapitch_p(it+1)=GK*Kp*(omega(it)-omega_reference);
    thetapitch_i(it+1)=thetapitch_i(it)+GK*Ki*(omega(it)-omega_reference)*delt;
    thetapitch(it+1)=thetapitch_p(it+1)+thetapitch_i(it+1);
    if thetapitch(it+1)>(thetapitch(it)+thetamaxdot*delt) && thetapitch(it+1)<thetapitchmax
        thetapitch(it+1)=thetapitch(it)+thetamaxdot*delt;
    elseif thetapitch(it+1)<(thetapitch(it)-thetamaxdot*delt) && thetapitch(it+1)>thetapitchmin
        thetapitch(it+1)=thetapitch(it)-thetamaxdot*delt;
    elseif thetapitch(it+1)>=thetapitchmax
        thetapitch(it+1)=thetapitchmax;
    elseif thetapitch(it+1)<=thetapitchmin
        thetapitch(it+1)=thetapitchmin;
    end
    omega(it+1)=omega(it)+((totalMR(it)-Mg(it))/I_rot)*delt;
end

%% Delete the last element of omega and thetapitch
omega_final=omega(1:end-1);
thetapitch_final=thetapitch(1:end-1);

%% Save power and thrust

save('Q1_answers.mat', 'Output','Mg','totalP','totalT','omega_final','thetapitch_final','t');

%% Plot for the generator torque without PI control

% figure;
% sizef=16;
% plot(omega_nopi,Mg_nopi,'linewidth',1);
% ylabel({'Generator Torque [Nm]'},'Interpreter','latex','fontsize',sizef);
% xlabel({'Rotational Speed [rad/s]'},'Interpreter','latex','fontsize',sizef);
% print('-depsc','Mg_noPI');

%% Plot the generator torque, the pitch angle, the rotational velocity and the power with PI control

figure;
sizef=16;
plot(omega_final,Mg,'linewidth',1);
ylabel({'Generator Torque [Nm]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Rotational Speed [rad/s]'},'Interpreter','latex','fontsize',sizef);
print('-depsc','Mg_PI');

figure;
sizef=16;
plot(t,omega_final,'linewidth',1);
ylabel({'Rotational speed [rad/s]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
print('-depsc','omega_PI');

figure;
sizef=16;
plot(t,(thetapitch_final*180)/pi,'linewidth',1);
ylabel({'Pitch angle [deg]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
print('-depsc','pitch_PI');

figure;
sizef=16;
plot(t,totalP,'linewidth',1);
ylabel({'Power [W]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
print('-depsc','power_PI');

