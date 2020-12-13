%% clean up
clear; close all; clc;

%% initialization
load blade_table.mat    
load airfoildata.mat 
load u.mat

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
delt=0.1;               %time step [s] 
n=1:1:2048;             %element in time []   
t=n*delt;               %time [s]
nt=numel(t); nb=3; nr=numel(r); %number of elements in: time, blade and blade elements
Wy0=0; Wz0=0;           %first value for induced wind in each direction
% L1, L2 and L3 will be defined in the box coordinate system like in the
% manual
L1=1637.6;              %[m] length in the flow direction of the turbulence field
L2=204;                 %[m] length in the horizontal direction of the turbulence field
L3=204;                 %[m] length in the vertical direction of the turbulence field

%Number of grid points in 1,2,3 directions
N1=2048;
N2=256;
N3=256;

%% preallocating
Wy=0;
Wz=0;
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
thetablade=[omega*t omega*t+(2*pi)/3 omega*t+(4*pi)/3];
Output=struct('xyzr',[],'T',[],'P',[],'MR',[]);
xyzr=zeros(3,nr);
phi=zeros(1,nr);
alpha=zeros(1,nr);
Power3=zeros(1,nt);
Thrust3=zeros(1,nt);

%% Transformation to the different coordinate systems

% wind speed in coordinate system 1 in each direction
% V1z will be calculated afterwards from the turbulence box
V1y=0;   
V1x=0;

% tranformation matrix from system 1 to 2
a12=zeros(3); a12(1,1)=cos(thetatilt); a12(1,2)=0; a12(1,3)=-sin(thetatilt);
a12(2,1)=sin(thetatilt)*sin(thetayaw); a12(2,2)=cos(thetayaw); a12(2,3)=sin(thetayaw)*cos(thetatilt);
a12(3,1)=sin(thetatilt)*cos(thetayaw); a12(3,2)=-sin(thetayaw); a12(3,3)=cos(thetatilt)*cos(thetayaw);
a21=a12';

%% loops 

for it=1:nt % iteration over time
    
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
            
            % tranformation matrix from system 1 to 4 and 4 to 1
            a14_b=(a34*a23_b)*a12;   a41_b=a14_b';
            
            % posistion calculations of the elements
            rt=[H,0,0]; rs=a21*[0,0,-Ls]';
            rb=a41_b*[r(ir),0,0]'; rr=rt+rs+rb;
            
            % position in blade system
            xb=rt(1)+rs(1)+rb(1);
            yb=rt(2)+rs(2)+rb(2);
            zb=rt(3)+rs(3)+rb(3);
            
           % here we will create the xbinary, ybinary vectors corresponding to the phyisical position
           
           delta1=L1/(N1-1);
           delta2=L2/(N2-1);
           delta3=L3/(N3-1);
           
           for ip=1:N2 %N3 can also be used because N2=N3
               %All these coordinates are in the blade coordinate system
               xbinary(ip)=(ip-1)*(delta3)+(H-L3/2); %this direction is upwards
               ybinary(ip)=(ip-1)*(delta2)-L2/2;
               
           end

           %interpolation to find the velocity in the exact xb, yb, zb 
           u_int=squeeze(u(it,:,:));
           V1z=interp2(ybinary',xbinary',u_int,yb,xb); %
           V1z=V1z+V0;
           
           
            
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
            p_n(it,ir)=L(ir)*cos(phi(ir))+D(ir)*sin(phi(ir));
            p_t(it,ir)=L(ir)*sin(phi(ir))-D(ir)*cos(phi(ir));
            
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
        
        % the normal and tangential loads are defined to be 0 at the tip
        p_n(it,nr)=0; p_t(it,nr)=0;
        
        % computing torque, power and thrust
        MR=trapz(r,p_t(it,:).*r);
        P=omega*MR;
        T=trapz(r,p_n(it,:)');
        
        % saving torque, power, thrust and velocity matrix onto struct
        Output(it,ib).MR=MR;
        Output(it,ib).P=P;
        Output(it,ib).T=T;
        Output(it,ib).xyzr=xyzr;
        
        % saving power and thrust for each timestep as vector
        Power3(it)=Output(it).P;
        Thrust3(it)=Output(it).T;
    end
end

%% save power and thrust

save('Q3_answers.mat', 'Power3','Thrust3', 'p_n','p_t','Output','t'); 

%% Total power and thrust calculation

for iii=1:2048
    totalP(iii)=Output(iii,1).P+Output(iii,2).P+Output(iii,3).P;
    totalT(iii)=Output(iii,1).T+Output(iii,2).T+Output(iii,3).T;
end

%% plot loads, thrust and power time series

figure;
sizef=16;
plot(t,p_n(:,9),'linewidth',1);
ylabel({'Normal Load [Nm]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
xlim([1 t(end)]);
print('-depsc','pnQ31');

figure;
plot(t,totalT,'linewidth',1);
ylabel({'Thrust [N]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
xlim([1 t(end)]);
print('-depsc','ThrustQ31');

figure;
plot(t,totalP,'linewidth',1);
ylabel({'Power [W]'},'Interpreter','latex','fontsize',sizef);
xlabel({'Time [s]'},'Interpreter','latex','fontsize',sizef);
xlim([1 t(end)]);
print('-depsc','PowerQ31');
%% PSD calculation for the power

pmean=mean(totalP);
NTime=2048;
delta_t=delt;
f_low=1/t(NTime);
f_high=0.5/delta_t;
fs=1/delta_t;
f=f_low:0.01:f_high;
no=2; 
figure;
[Pyy,f]=pwelch(totalP-pmean,no*128,no*64,f,fs);
semilogy(2*pi*f/omega,Pyy,'linewidth',0.8)
line([3 3],[10^0 10^(14)],'linestyle',':','linewidth',1.5,'color','r');
line([6 6],[10^0 10^(14)],'linestyle',':','linewidth',1.5,'color','r');
line([9 9],[10^0 10^(14)],'linestyle',':','linewidth',1.5,'color','r');
grid on
xlim([0 10])
ylim([10^8 10^14])
xlabel('\omega/\omega_{o}','FontSize',14)
ylabel('PSD of Power [W^{2}/Hz]','FontSize',14)
set(gca,'FontSize',14)
print('-depsc','PSD_PowerQ31_sum');

%% PSD calculation for the thrust

tmean=mean(totalT);
no=2;
figure;
[Tyy,f]=pwelch(totalT-tmean,no*128,no*64,f,fs);
semilogy(2*pi*f/(2*pi/9.3361),Tyy,'linewidth',0.9)
hold on
line([3 3],[10^0 10^(12)],'linestyle',':','linewidth',1.5,'color','r');
line([6 6],[10^0 10^(12)],'linestyle',':','linewidth',1.5,'color','r');
line([9 9],[10^0 10^(12)],'linestyle',':','linewidth',1.5,'color','r');
grid on
xlim([0 10])
ylim([10^6 10^12])
xlabel('\omega/\omega_{o}','FontSize',14)
ylabel('PSD of the thrust','FontSize',14)
hold off 

set(gca,'FontSize',14)
print('-depsc','PSD_ThrustQ31_sum');
%% PSD calculation for the normal loads

nlmean=mean(p_n(:,9)); % We are considering the 9th element of the blade

no=4;
figure;
[Nlyy,f]=pwelch(p_n(:,9)-nlmean,no*128,no*64,f,fs);
semilogy(2*pi*f/omega,Nlyy,'linewidth',0.8)
line([1 1],[10^0 10^(14)],'linestyle',':','linewidth',1.5,'color','r');
line([2 2],[10^0 10^(14)],'linestyle',':','linewidth',1.5,'color','r');
line([3 3],[10^0 10^(14)],'linestyle',':','linewidth',1.5,'color','r');
grid on
xlim([0 4])
ylim([10^3 10^8])
xlabel('\omega/\omega_{o}','FontSize',14)
ylabel('PSD of the normal load','FontSize',14)
set(gca,'FontSize',14)
print('-depsc','PSD_Normal_loadsQ31_bis');