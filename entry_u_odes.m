function outp=entry_u_odes(t,inp) % ,tauf_init
global GM R0 g0 Vs Ts hs rho0 mass Sr CD CL Rn kQ ...
    u_min u_max c0 c1 h0 r0 V0 gamma0 Vf lambda_rf lambda_gammaf Hf ...
    qmax Qmax gmax ec eq eQ eg

%%
r=inp(1);
V=inp(2);
gamma=inp(3);

h=(r-1)*R0;
rho=rho0*exp(-h/hs);
% q=rho*(V*Vs)^2/2;
% aq=q/qmax;
% Q=kQ*sqrt(rho/Rn)*(V*Vs)^3;
% aQ=Q/Qmax;
% g=R0*rho*V^2*Sr*sqrt(CL^2+CD^2)/(2*mass);
% ag=g/gmax;

L=R0*rho*V^2*Sr*CL/(2*mass);
D=R0*rho*V^2*Sr*CD/(2*mass);

u=0;

rdot=V*sin(gamma);
Vdot=-D-sin(gamma)/r^2;
gammadot=L*u/V + V*cos(gamma)/r - cos(gamma)/r^2/V;


outp=[rdot, Vdot, gammadot].'; %fsolve
end