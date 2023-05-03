
delta_t = 0.1;
x_dot_ref = 80;
x_dot0 = 20;
m = 2050;
g=9.81; 
mu = 0.3;
rho = 1.225;
S = 1;
Cd = 0.3;
kd = 1/2*rho*S*Cd;

t_end = 100;
N = t_end/delta_t;
t = linspace(0,t_end,N);

error = zeros(1,N);
beta = zeros(1,N);
x_dot = zeros(1,N+1);
x_dot(1)= x_dot0;
Kp = 0.3;
Ki = 0.00003;
error_integral = 0;
for i=1:N
   error(i) = x_dot_ref-x_dot(i);
   error_integral = error_integral+error(i);
   beta(i) = Kp*error(i)+Ki*error_integral;
   beta(i) = max( min(beta(i),1),-1);
   x_dot(i+1) = x_dot(i)+delta_t*(mu*g*beta(i)-kd/m*(x_dot(i))^2);
end

figure(1)
plot(t,x_dot(1:end-1),'LineWidth',2)
hold on
plot(t,x_dot_ref*ones(1,N),'LineWidth',2)
hold on
plot(t,error,'LineWidth',2)
legend('x_dot','reference','error')
grid on
figure(2),
plot(t,beta,'LineWidth',2)
grid minor