B_sep = 6000*0.0018;
C_sep = 0.4*2.27;
N = 1000;
alpha = linspace(-pi/18,pi/18,N);

beta = 0.4;
% From http://www.thecartech.com/subjects/auto_eng/Center_of_Gravity.htm
b = 1.52;
a = 0.92;
e = 1.112;

mu = 0.3;
g = 9.81;
m = 2050;

%Linearize eqn 8

Fyf = mu*(b-e*beta)*m*g/(a+b)*sqrt(1-beta^2)*sin(C_sep*atan(B_sep*alpha));
Fyr = mu*(a+e*beta)*m*g/(a+b)*sqrt(1-beta^2)*sin(C_sep*atan(B_sep*alpha));
n = 100;
Lf = polyfit(alpha,Fyf,1);
Uf = polyfit(alpha(N/2-n:N/2+n),Fyf(N/2-n:N/2+n),1);
CLf = Lf(1)*alpha+Lf(2);
CUf = Uf(1)*alpha+Uf(2);

Lr = polyfit(alpha,Fyr,1);
Ur = polyfit(alpha(N/2-n:N/2+n),Fyr(N/2-n:N/2+n),1);
CLr = Lr(1)*alpha+Lr(2);
CUr = Ur(1)*alpha+Ur(2);

figure(1),
plot(180/pi*alpha,Fyf,'k','LineWidth',3)
hold on
plot(180*alpha/pi,CLf,'r--','LineWidth',2)
hold on
plot(180*alpha/pi,CUf,'b:','LineWidth',2)
grid on
xlabel('alpha')
ylabel('F_y front')
legend('Nonlinear','Lower','Upper')

figure(2),
plot(180/pi*alpha,Fyr,'k','LineWidth',3)
hold on
plot(180*alpha/pi,CLr,'r--','LineWidth',2)
hold on
plot(180*alpha/pi,CUr,'b:','LineWidth',2)
grid on
xlabel('alpha')
ylabel('F_y reverse')
legend('Nonlinear','Lower','Upper')
grid on