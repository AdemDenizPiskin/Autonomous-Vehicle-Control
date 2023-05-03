m = 2050;
g = 9.81; 
I = 3344;
mu = 0.3;
b = 1.52;
a = 0.92;
e = 1.112;
beta =1;
rho = 1.225;
S = 1;
Cd = 0.3;
kd = 1/2*rho*S*Cd;

B_sep = 2500*0.0018;
C_sep = 1*2.27;
psi_r =1/10000;
u = -0*pi/180;
tspan = [0 5];
y0 = [0 0 0 0 0 20 0];
params = [m g mu a b e B_sep C_sep I psi_r kd]';
%z = vehicle_dynamics(0,y0,params,u,beta)

[t,x] = ode45(@(t,x)vehicle_dynamics(t,x,params,u,beta), tspan, y0);
figure,
plot(t,x(:,1))
title('y_dot')
figure,
plot(t,x(:,2)*180/pi)
title('psi_dot')
figure,
plot(t,x(:,3)*180/pi)
title('e_psi')
figure,
plot(t,x(:,4))
title('e_y')
figure,
plot(t,x(:,5)*180/pi)
title('delta')
figure,
plot(t,x(:,6))
title('x_dot')
figure,
plot(t,x(:,7))
title('s')
figure,
%y_road = road_cord_generator(road_curve,lamda,x(:,7)');
[x_road , y_road] = generate_circular_road(10000,x(:,7)');
plot(x_road,y_road)
hold on
plot(x_road,y_road+lw/2)
hold on
plot(x_road,y_road-lw/2)
hold on
plot(x_road,x(:,4)'+y_road)
title('Road')
legend('Center lane','left','right','car')
function [x,y] = generate_circular_road(radius,arclength)
x = zeros(1,length(arclength));
y = zeros(1,length(arclength));
for i=1:length(arclength)
    theta = arclength(i)/radius;
    alpha = (pi-theta)/2;
    l = radius*sqrt(2*(1-cos(theta)));
    x(i) = l*sin(alpha);
    y(i) = l*cos(alpha);
   
end

end
