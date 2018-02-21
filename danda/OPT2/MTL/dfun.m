function [dstate] = dfun(t,state,control,X)

T = 4e-3;
g0 = 1.6e-3;
rho0 = 0.1;
Sref = 1;
Cd = 1;
Isp = 937.5;

dstate = zeros(5,1);

u = control(t,X);
rho = rho0*exp(-state(3));

dstate(1) = state(2);
dstate(2) = (T*cos(u)-0.5*rho*Sref*Cd*sqrt(state(2)^2 + state(4)^2)*state(2))/state(5);
dstate(3) = state(4)    ;
dstate(4) = (T*sin(u)-0.5*rho*Sref*Cd*sqrt(state(2)^2 + state(4)^2)*state(4))/state(5)-g0;
dstate(5) = -T/(g0*Isp) ;


    
end
