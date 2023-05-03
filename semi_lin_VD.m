    function out = semi_lin_VD(t,x,params,u)
out = zeros(5,1);
m = params(1);
g = params(2);
mu = params(3);
a = params(4);
b = params(5);
e = params(6);

B = params(7);
C = params(8);
I = params(9);
psi_r = params(10);
x_dot = params(11);
kd = params(12);
beta =  kd*x_dot^2/(m*g*mu);


y_dot = x(1);
psi_dot = x(2);
e_psi = x(3);
e_y =x(4);
delta = x(5);


alpha_f = (y_dot+a*psi_dot)/x_dot-delta;
alpha_r = (y_dot-b*psi_dot)/x_dot;

Fyf = mu*(b-e*beta)*m*g/(a+b)*sqrt(1-beta^2)*sin(C*atan(B*alpha_f));
Fyr = mu*(a+e*beta)*m*g/(a+b)*sqrt(1-beta^2)*sin(C*atan(B*alpha_r));


out(1) =( -m*x_dot*psi_dot+Fyf+Fyr)/m;
out(2) = (a*Fyf-b*Fyr)/I;
out(3) = psi_dot-x_dot*psi_r;
out(4) = y_dot+x_dot*e_psi;
out(5) = u;

end




