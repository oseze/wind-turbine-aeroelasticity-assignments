%% clean up
clear; close all; clc;

%% initialize
load blade_table.mat    
load airfoildata.mat    
R=89.17;                %rotor radius[m]
B=3;                    %number of blade[]
P_rated=10000;          %rated power[kW]
V0_cutin=4;             %cut in speed [m/s]
V0_cutout=25;           %cut out speed[m/s]
V0=8;                   %incoming wind speed[m/s]
omega=0.673;            %rotational speed[rad/s]
rho=1.225;              %air density[kg/m^3]
thetacone=(0*pi)/180;   %[rad]
thetayaw=(0*pi)/180;    %[rad]
thetapitch=(0*pi)/180;  %[rad]
thetatilt=(0*pi)/180;   %[rad]
H=119;                  %height of tower[m]
Ls=7.1;                 %shaft length[m]
delt=0.02;              %time step [s] 
n=1:1:10000;            %element in time []   
t=n*delt;               %time [s]
nt=numel(t); nb=3; nr=numel(r); %number of elements in: time, blade and blade elements
A=4;                    %coefficient  (dynamic stall)      
k=0.6;                  %differential equation coefficient
fs0=0;                  %first value for flow seperation
Wy0=0; Wz0=0;           %first value for induced wind in each direction

%% preallocating
Wy=zeros(nt+1,nr); 
Wz=zeros(nt+1,nr); 
Wyquasi=zeros(nt+1,nr); 
Wzquasi=zeros(nt+1,nr);
Winty=zeros(nt+1,nr); 
Wintz=zeros(nt+1,nr);
Vrely=zeros(nt,nr);
Vrelz=zeros(nt,nr);
thetablade=[omega*t omega*t+(2*pi)/3 omega*t+(4*pi)/3];
Output=struct('xyzr',[],'T',[],'P',[],'MR',[]);
xyzr=zeros(3,nr);
normVrel=zeros(1,nt);
L=zeros(nt,nr);
D=zeros(nt,nr);
p_n=zeros(1,nr);
p_t=zeros(1,nr);
Wy=zeros(nt+1,nr);
Wz=zeros(nt+1,nr);
a23_b=zeros(3); 
a34=zeros(3);
xb=zeros(1,nt);
yb=zeros(1,nt);
zb=zeros(1,nt);
V4=zeros(3,1);
V1=zeros(3,1);
phi=zeros(1,nr);
alpha=zeros(1,nr);
fs=zeros(1,nt);
Cl=zeros(1,nt);
Power2=zeros(1,nt);
Thrust2=zeros(1,nt);
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

%% loops

for it=1:nt % iteration over time
    
    % compuation of the pitch step evolution
    if t(it) < 100 || t(it) > 150
        thetapitch=0*pi/180; 
    else
        thetapitch=2*pi/180;
    end
    
    for ib=1:nb % iteration across blade number
        
        for ir=1:nr-1 % iteration across blade element
            
            % tranformation matrix from system 2 to 3
            a23_b(1,:)=[cos(thetablade(ib)),sin(thetablade(ib)),0];
            a23_b(2,:)=[-sin(thetablade(ib)),cos(thetablade(ib)),0];
            a23_b(3,:)=[0,0,1];
            
            % tranformation matrix from system 3 to 4
            a34(1,:)=[cos(thetacone),0,-sin(thetacone)];
            a34(2,:)=[0,1,0];
            a34(3,:)=[sin(thetacone),0,cos(thetacone)];
            
            % tranformation matrix from system 1 to 4, and 4 to 1    
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
                Vrely(it,ir)=V4(2)+Wy0-omega*r(ir)*cos(thetacone);
                Vrelz(it,ir)=V4(3)+Wz0;
            else
                a=abs(Wz(it,ir)/V0);
                Vrely(it,ir)=V4(2)+Wy(it,ir)-omega*r(ir)*cos(thetacone);
                Vrelz(it,ir)=V4(3)+Wz(it,ir);
            end
            
            % flow angle and angle of attack calculation
            phi(ir)=atan(Vrelz(it,ir)/(-Vrely(it,ir))); %in radians
            alpha(ir)=(phi(ir)*180/pi)-(-beta(ir)+thetapitch*180/pi); %in degrees
            
            % interpolation to get flow seperation function
            thick=[100, 60, 48, 36, 30.1, 24.1];
            fs1=interp1(W3_100(:,1),W3_100(:,5),alpha(ir));
            fs2=interp1(W3_60(:,1),W3_60(:,5),alpha(ir));
            fs3=interp1(W3_48(:,1),W3_48(:,5),alpha(ir));
            fs4=interp1(W3_36(:,1),W3_36(:,5),alpha(ir));
            fs5=interp1(W3_30(:,1),W3_30(:,5),alpha(ir));
            fs6=interp1(W3_24(:,1),W3_24(:,5),alpha(ir));
            fsvec=[fs1 fs2 fs3 fs4 fs5 fs6];
            fss=interp1(thick,fsvec,tc(ir));
            
            % interpolation to get the inviscid lift coefficient
            cli1=interp1(W3_100(:,1),W3_100(:,6),alpha(ir));
            cli2=interp1(W3_60(:,1),W3_60(:,6),alpha(ir));
            cli3=interp1(W3_48(:,1),W3_48(:,6),alpha(ir));
            cli4=interp1(W3_36(:,1),W3_36(:,6),alpha(ir));
            cli5=interp1(W3_30(:,1),W3_30(:,6),alpha(ir));
            cli6=interp1(W3_24(:,1),W3_24(:,6),alpha(ir));
            clivec=[cli1 cli2 cli3 cli4 cli5 cli6];
            Cli=interp1(thick,clivec,tc(ir));
            
            % interpolaion to get the lift coefficient with full seperation
            clfs1=interp1(W3_100(:,1),W3_100(:,7),alpha(ir));
            clfs2=interp1(W3_60(:,1),W3_60(:,7),alpha(ir));
            clfs3=interp1(W3_48(:,1),W3_48(:,7),alpha(ir));
            clfs4=interp1(W3_36(:,1),W3_36(:,7),alpha(ir));
            clfs5=interp1(W3_30(:,1),W3_30(:,7),alpha(ir));
            clfs6=interp1(W3_24(:,1),W3_24(:,7),alpha(ir));
            clfsvec=[clfs1 clfs2 clfs3 clfs4 clfs5 clfs6];
            Clfs=interp1(thick,clfsvec,tc(ir));
            
            % norm of relative velocity 
            normVrel(it,ir)=sqrt((Vrely(it,ir))^2+(Vrelz(it,ir))^2);
            
            % tau - characteristic time used in dynamic stall model
            tau=A*c(ir)/normVrel(it,ir);
            
            % compute flow seperation function, if statement to ensure the
            %function is defined also for t=1.
            if it==1
                fs(it)=fss+(fs0-fss)*exp(-delt/tau); 
            else
                fs(it)=fss+(fs(it-1)-fss)*exp(-delt/tau);
            end
            
            % lift coefficient after dynamic stall
            Cl(it)=fs(it)*Cli+(1-fs(it))*Clfs;
            
            % interpolaion to get the drag coefficient
            cd1=interp1(W3_100(:,1),W3_100(:,3),alpha(ir));
            cd2=interp1(W3_60(:,1),W3_60(:,3),alpha(ir));
            cd3=interp1(W3_48(:,1),W3_48(:,3),alpha(ir));
            cd4=interp1(W3_36(:,1),W3_36(:,3),alpha(ir));
            cd5=interp1(W3_30(:,1),W3_30(:,3),alpha(ir));
            cd6=interp1(W3_24(:,1),W3_24(:,3),alpha(ir));
            cdvec=[cd1 cd2 cd3 cd4 cd5 cd6];
            Cd=interp1(thick,cdvec,tc(ir));
            
            % computation of lift, drag, tangential load and normal load
            L(ir)=0.5*rho*normVrel(it,ir)^2*c(ir)*Cl(it);
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
            
            % norm of
            if it==1
                norm=sqrt((V4(2))^2+(V4(3)+fg*Wz0)^2);
            else
                norm=sqrt((V4(2))^2+(V4(3)+fg*Wz(it,ir))^2);
            end
            
            % in each time loop the value of W is being calculated for the
            % next loop - the solution converges in time. Dynamic flow
            % implementation
            Wzquasi(it+1,ir)=-(B*L(ir)*cos(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
            Wyquasi(it+1,ir)=-(B*L(ir)*sin(phi(ir)))/(4*pi*rho*r(ir)*F*norm);
         
            tau1=(1.1/(1-1.3*a))*R/V0;
            tau2=(0.39-0.26*(r(ir)/R)^2)*tau1;

            Hz=Wzquasi(it+1,ir)+k*tau1*((Wzquasi(it+1,ir)-Wzquasi(it,ir))/delt);
            Hy=Wyquasi(it+1,ir)+k*tau1*((Wyquasi(it+1,ir)-Wyquasi(it,ir))/delt);  
            Wintz(it+1,ir)=Hz+(Wintz(it,ir)-Hz)*exp(-delt/tau1);
            Winty(it+1,ir)=Hy+(Winty(it,ir)-Hy)*exp(-delt/tau1);
            Wz(it+1,ir)=Wintz(it+1,ir)+(Wz(it,ir)-Wintz(it+1,ir))*exp(-delt/tau2);
            Wy(it+1,ir)=Winty(it+1,ir)+(Wy(it,ir)-Winty(it+1,ir))*exp(-delt/tau2);   
        end
        % the normal and tangetial loads are defined to be 0 at the tip
        p_n(nr)=0; p_t(nr)=0;
        
        % computing torque, power and thrust 
        MR=trapz(r,p_t.*r);
        P=omega*MR;
        T=trapz(r,p_n');
        
        % saving torque, power, thrust and velocity matrix onto struct
        Output(it,ib).MR=MR;
        Output(it,ib).P=P;
        Output(it,ib).T=T;
        Output(it,ib).xyzr=xyzr;
        
        % saving power and thrust for each timestep as vector
        Power2(it)=Output(it).P;
        Thrust2(it)=Output(it).T;
        
    end
end
%remove last element in induced wind 
Wy(10001,:)=[];
Wz(10001,:)=[];
%% save

save('Q2_answers.mat', 'Power2','Thrust2','Wy','Wz','p_n','p_t'); 

%% load and define

load Q1_answers.mat
load Q2_answers.mat 
%% Parameters for the plots
m=90/delt;
n=160/delt;
sizef=16;

%% plot thrust

figure;
plot(t,Thrust2,'linewidth',1)

% set(legend,'Interpreter','latex','fontsize',sizef);
ylabel({'Thrust [N]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
% hold off
print('-depsc','thrust2');

%% plot power

figure;
plot(t,Power2,'linewidth',1)

% set(legend,'Interpreter','latex','fontsize',sizef);
ylabel({'Power [N]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
% hold off
print('-depsc','Power2');


%% plot induced velocity and power

W=sqrt(Wz(:,9).^2+Wy(:,9).^2); %norm of induced wind velocity
figure;
subplot(2,1,1) 
plot(t,W,'linewidth',1)
ylabel({'Induced Wind [m/s]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
hold on
subplot(2,1,2)
plot(t,Power2,'linewidth',1)
ylabel({'Power [W]'},'Interpreter','latex','fontsize',sizef);
hold off
print('-depsc','WandP');
%% plot induced velocity zoomed

figure;
plot(t(m:n),W(m:n),'linewidth',1);
ylabel({'Induced Wind [m/s]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
print('-depsc','WQ2');


%% plot loads 

figure;
plot(r,p_n,'linewidth',1.2)
ylabel({'Normal load [N]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Radius [m]'},'Interpreter','latex','fontsize',sizef);
print('-depsc','pnQ2');

figure;
plot(r,p_t,'linewidth',1.2)
ylabel({'Tangential load [N]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Radius [m]'},'Interpreter','latex','fontsize',sizef);
%xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);

print('-depsc','ptQ2');
%% plot thrust and induced wind speed

figure;
subplot(2,1,1) 
plot(t,W,'linewidth',1)
ylabel({'Induced Wind [m/s]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
hold on
subplot(2,1,2)
plot(t,Thrust2,'linewidth',1)
ylabel({'Thrust [N]'},'Interpreter','latex','fontsize',sizef);
hold off
print('-depsc','WandT');

