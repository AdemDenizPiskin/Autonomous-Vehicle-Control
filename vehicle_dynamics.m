function out = vehicle_dynamics(t,x,params,u,beta)
out = zeros(7,1);
m = params(1);
g = params(2);
mu = params(3);
a = params(4);
b = params(5);
e = params(6);
c = a+b;
B = params(7);
C = params(8);
I = params(9);
psi_r = params(10);
rho = 1/psi_r;
kd = params(11);


y_dot = x(1);
psi_dot = x(2);
e_psi = x(3);
e_y =x(4);
delta = x(5);
x_dot = x(6);
s = x(7);

Fx = m*g*mu*beta;
alpha_f = (y_dot+a*psi_dot)/x_dot-delta;
alpha_r = (y_dot-b*psi_dot)/x_dot;

Fyf = mu*(b-e*beta)*m*g/(a+b)*sqrt(1-beta^2)*sin(C*atan(B*alpha_f));
Fyr = mu*(a+e*beta)*m*g/(a+b)*sqrt(1-beta^2)*sin(C*atan(B*alpha_r));
Fy = Fyr+Fyf;
Fz = m*g;

Fzfl = (b*Fz-e*Fx)/(2*a+2*b)-e*Fy/(2*c);
Fzfr = (b*Fz-e*Fx)/(2*a+2*b)+e*Fy/(2*c);
Fzrl = (a*Fz+e*Fx)/(2*a+2*b)-e*Fy/(2*c);
Fzrr = (a*Fz+e*Fx)/(2*a+2*b)+e*Fy/(2*c);
fxfl = mu*beta*Fzfl;    
fxfr = mu*beta*Fzfr;
fxrl = mu*beta*Fzrl;
fxrr = mu*beta*Fzrr;
fyfl = sqrt((mu*Fzfl)^2-fxfl^2)*sin(C*atan(B*alpha_f));
fyfr = sqrt((mu*Fzfr)^2-fxfr^2)*sin(C*atan(B*alpha_f));
fyrl = sqrt((mu*Fzrl)^2-fxrl^2)*sin(C*atan(B*alpha_r));
fyrr = sqrt((mu*Fzrr)^2-fxrr^2)*sin(C*atan(B*alpha_r));

Fxrl = fxrl;
Fxrr = fxrr;
Fyrl = fyrl;
Fyrr = fyrr;
Fxfl = fxfl*cos(delta)-fyfl*sin(delta);  
Fxfr = fxfr*cos(delta)-fyfr*sin(delta);
Fyfl = fxfl*sin(delta)+fyfl*cos(delta);
Fyfr = fxfr*sin(delta)+fyfr*cos(delta);

out(1) =( -m*x_dot*psi_dot+Fyfl+Fyfr+Fyrl+Fyrr)/m;
out(2) = (a*(Fyfl+Fyfr)-b*(Fyrl+Fyrr))/I;
out(3) = psi_dot-x_dot*psi_r;
out(4) = y_dot*cos(e_psi)+x_dot*sin(e_psi);
out(5) = u;
out(6) = (Fxfl+Fxfr+Fxrl+Fxrr-kd*x_dot^2)/m;
if psi_r == 0
    out(7) = x_dot*cos(e_psi)-y_dot*sin(e_psi);
else
    out(7) = rho/(rho-e_y)*(x_dot*cos(e_psi)-y_dot*sin(e_psi));
end
end




