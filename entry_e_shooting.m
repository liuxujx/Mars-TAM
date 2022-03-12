function outp=entry_e_shooting(inp)

global GM R0 g0 Vs Ts hs rho0 mass Sr CD CL Rn kQ ...
    u_min u_max c0 c1 h0 r0 V0 gamma0 Vf lambda_rf lambda_gammaf Hf ...
    qmax Qmax gmax ec eq eQ eg e

lambda0=inp(1:3);
tauf=inp(4);
Y0=[r0,V0,gamma0, lambda0];

[~,Yout]=ode45(@entry_e_lambda,[0,tauf],Y0); % Yf=N×7
Yf=Yout(end,:);
% 求H(tauf)，但忽略下标tauf
r=Yf(1);
V=Yf(2);
gamma=Yf(3);
lambda_r=Yf(4);
lambda_V=Yf(5);
lambda_gamma=Yf(6);

h=(r-1)*R0;
rho=rho0*exp(-h/hs);
q=rho*(V*Vs)^2/2;
aq=q/qmax;
Q=kQ*sqrt(rho/Rn)*(V*Vs)^3.15;
aQ=Q/Qmax;
g=R0*rho*V^2*Sr*sqrt(CL^2+CD^2)/(2*mass);
ag=g/gmax;

L=R0*rho*V^2*Sr*CL/(2*mass);
D=R0*rho*V^2*Sr*CD/(2*mass);

u=-lambda_gamma*L/V/(1-e); %arctan
if u>u_max
    u=u_max;
elseif u<u_min
    u=u_min;
end

rdot=V*sin(gamma);
Vdot=-D-sin(gamma)/r^2;
gammadot=L*u/V + V*cos(gamma)/r - cos(gamma)/r^2/V;

%H(tauf)表达式
H = lambda_r*rdot + lambda_V*Vdot + lambda_gamma*gammadot + (1-e)/2*u^2;

outp=[Yf(2)-Vf;
    Yf(4)-lambda_rf;
    Yf(6)-lambda_gammaf;
    H-Hf];
end