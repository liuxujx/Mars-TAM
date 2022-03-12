function outp=entry_u_lambda(t,inp) % ,tauf_init
global d2r r2d GM R0 g0 Vs Ts hs rho0 mass Sr CD CL Rn kQ ...
    u_min u_max c0 c1 h0 r0 V0 gamma0 Vf lambda_rf lambda_gammaf Hf ...
    qmax Qmax gmax ec eq eQ eg

%%
r=inp(1);
V=inp(2);
gamma=inp(3);
lambda_r=inp(4);
lambda_V=inp(5);
lambda_gamma=inp(6);

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

u=0; % u=cos(sigma)==0
% u=atan(lambda_gamma*L*c1/V/ec); %arctan

rdot=V*sin(gamma);
Vdot=-D-sin(gamma)/r^2;
gammadot=L*u/V + V*cos(gamma)/r - cos(gamma)/r^2/V;
lambda_rdot=lambda_gamma*((V*cos(gamma))/r^2 - (2*cos(gamma))/(V*r^3) + (CL*R0^2*Sr*V*rho0*u*exp(-(R0*(r - 1))/hs))/(2*hs*mass)) - lambda_V*((2*sin(gamma))/r^3 + (CD*R0^2*Sr*V^2*rho0*exp(-(R0*(r - 1))/hs))/(2*hs*mass));

lambda_Vdot=(CD*R0*Sr*V*lambda_V*rho0*exp(-(R0*(r - 1))/hs))/mass - lambda_r*sin(gamma) - lambda_gamma*(cos(gamma)/r + cos(gamma)/(V^2*r^2) + (CL*R0*Sr*rho0*u*exp(-(R0*(r - 1))/hs))/(2*mass));

lambda_gammadot=lambda_gamma*((V*sin(gamma))/r - sin(gamma)/(V*r^2)) + (lambda_V*cos(gamma))/r^2 - V*lambda_r*cos(gamma);

outp=[rdot, Vdot, gammadot, lambda_rdot, lambda_Vdot, lambda_gammadot].'; %fsolve
end