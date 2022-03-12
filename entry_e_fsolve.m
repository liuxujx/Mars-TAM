%% 火星大气进入末段高度最大化，无约束
clc; clear all; close all;
global GM R0 g0 Vs Ts hs rho0 mass Sr CD CL Rn kQ ...
    u_min u_max c0 c1 h0 r0 V0 gamma0 Vf lambda_rf lambda_gammaf Hf ...
    qmax Qmax gmax ec eq eQ eg e

GM=4.284e13;
R0=3397000; %m
g0=GM/R0^2; %m/s^2
Vs=sqrt(R0*g0);
Ts=sqrt(R0/g0);

hs=9354.5;
rho0=0.0158;

mass=3300;
Sr=15.9;
CD=1.45;
CL=0.348;
Rn=0.6;
kQ=1.9027e-4;

u_min=cosd(120);
u_max=cosd(30);
c0=(u_max+u_min)/2;
c1=(u_max-u_min)/2;

qmax=10e3;
Qmax=70e4;
gmax=5*9.80665/g0;

%% 能量最优 J=0.5*∫u^2dτ 无约束

h0=125e3;
r0=1+h0/R0;
V0=6e3/Vs;
gamma0=-11.5*pi/180;
Vf=540/Vs;

lambda_rf=0;
lambda_gammaf=0;
Hf=0;
x0=[r0,V0,gamma0];

options=odeset('Events',@entry_u_events);
[tout,xout]=ode45(@entry_u_odes,[0,1],x0,options);
tauf=tout(end);
xf=xout(end,:);
rf=xf(1);
gammaf=xf(3);

% 由于lambda_rf=0 lambda_gammaf=0 Hf=0 导致lambda_Vf=0
% uf=0; % u=cos(sigma)==0
% hf=(rf-1)*R0;
% rhof=rho0*exp(-hf/hs);
% Lf=R0*0.5*rhof*Vf^2*Sr*CL/(2*mass);
% Df=R0*0.5*rhof*Vf^2*Sr*CL/(2*mass);
% rfdot=Vf*sin(gammaf);
% Vfdot=-Df-sin(gammaf)/rf^2;
% gammafdot=Lf*uf/Vf + Vf*cos(gammaf)/rf - cos(gammaf)/rf^2/Vf;
% lambda_Vf=(lambda_rf*rfdot+lambda_gammaf*gammafdot+0.5*uf^2)/Vfdot;%% Hf=0
lambda_Vf=0;
yf=[xf,lambda_rf,lambda_Vf,lambda_gammaf];
[~,yout]=ode45(@entry_u_lambda,[tauf,0],yf); % Yf=N×7
y0=yout(end,:);
lambda0=y0(4:6);
%% zheng yiyu 同伦法 无约束
ee=[1e-6;1e-4;1e-2;2e-2;4e-2;6e-2;8e-2;1e-1;3e-1;5e-1;7e-1;9e-1;
    9.1e-1;9.2e-1;9.3e-1;9.4e-1;9.5e-1;9.6e-1;9.7e-1;9.8e-1;
    9.9e-1;9.91e-1;9.95e-1];

h0=125e3;
r0=1+h0/R0;
V0=6e3/Vs;
gamma0=-11.5*pi/180;
Vf=540/Vs;

lambda_gammaf=0;
Hf=0;

Y0=[lambda0, tauf];
for n=1:length(ee)
    e=ee(n);
    lambda_rf=-e;
    [Y, fval, exitflag, output, jacobian] = fsolve(@entry_e_shooting,Y0); %,opt fsolve对初值非常敏感，使协态初值进行猜测
    Y0=Y;
end
y0=[r0,V0,gamma0,Y0(1:3)];
[tout,yout]=ode45(@entry_e_lambda,[0,Y0(4)],y0); % Yf=N×7
yf=yout(end,:);

%%
figure;
subplot(3,1,1); plot(tout,yout(:,1));
subplot(3,1,2); plot(tout,yout(:,2));
subplot(3,1,3); plot(tout,yout(:,3));
figure;
subplot(3,1,1); plot(tout,yout(:,4));
subplot(3,1,2); plot(tout,yout(:,5));
subplot(3,1,3); plot(tout,yout(:,6));
