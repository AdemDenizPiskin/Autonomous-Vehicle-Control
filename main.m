%% Vehicle Parameters

%Vehicle Dims From http://www.thecartech.com/subjects/auto_eng/Center_of_Gravity.htm
b = 1.52;
a = 0.92;
e = 1.112;

%Standart parameters from the paper
m = 2050;
g = 9.81; 
I = 3344;
mu = 0.3;
%Drag Force Cd is from wikipedia
rho = 1.225;
S = 1;
Cd = 0.3;
kd = 1/2*rho*S*Cd;

%Coefficients from  semi-empirical Pacejka formula modified to fit the
%paper
B_sep = 2500*0.0018;
C_sep = 1*2.27;

%Coefficients from the linearization of the Pacejka magic formula
CLf = 1.737251926751985e+04;
CLr = 2.205172460594410e+04;
CUf = 2.477806798654727e+04;
CUr = 3.145192260792386e+04;

%Actuator Limits
delta_lim =pi/3; %60deg
delta_dot_lim = pi;
af_lim = pi/36;
ar_lim = pi/36;

%Initializaitons
x_dot0 = 20;
x0 = 0;
y_dot0 = 0;
e_y0 = 0.05;
psi_dot0 = 0;
e_psi0 = 0; %5deg
delta0 = 0;

ksi0 = [y_dot0 psi_dot0 e_psi0 e_y0 delta0]';

%% Controller Parameters

delta_t = 0.01;
%PI controller
x_dot_ref = 40;
Kp = 0.3;
Ki = 0.00003;
%MPC
%Horizon
Hp = 45;
Hp2 = 20;
%Weights
Q = diag([0 0 1000000 1000000 0]);% 1000000
Qf = 10*Q;
R = 1000000;    
Q_hat_cm = Q;
R_hat_cm = R;
for i=1:Hp
    if i<Hp
    Q_hat_cm = blkdiag(Q_hat_cm,Q);
    R_hat_cm = blkdiag(R_hat_cm,R);
    else
    Q_hat_cm = blkdiag(Q_hat_cm,Qf);
    end
    
end

%Braking points
n_beta = 3;
%Post COmputation Weights
Qbeta = 1;
B = delta_t*[0 0 0 0 1]';
%% Road Parameters
%lane width
lw = 4.6; %says wikipedia
lamda = 90;
road_curve = 0.05;    %*cos(2*pi/50*x);
%Use those to find y to draw the road and also calculating psi_r
%y = road_cord_generator(road_curve,lamda,x); 
%psi_r = road_curvature_generator(road_curve,lamda,x);
%% Simulation Initialization
t_dur = 10;
n_dur = t_dur/delta_t;
t = linspace(0,t_dur,n_dur);
%states
s = zeros(1,n_dur+1);
x_dot = zeros(1,n_dur+1);
ksi = zeros(5,n_dur+1);
beta_ref = zeros(1,n_dur);
beta_opt = zeros(1,n_dur);
beta_vec = zeros(n_beta,n_dur);
delta_dot_vec = zeros(n_beta,n_dur);
delta_dot_opt = zeros(1,n_dur);
f_cost_vec = zeros(n_beta,n_dur);
error_x_dot = zeros(1,n_dur);
error_x_dot_integrator = 0;
x_dot_predictions = zeros(1,Hp+1);
s_predict = zeros(1,Hp+1);
%initial conditions
s(1) = x0;
x_dot(1) = x_dot0;
ksi(:,1) = ksi0;

for i=1:n_dur
   %%%Longtiudal Profiles Generator
   %Speed Control
   error_x_dot(i) = x_dot_ref-x_dot(i);
   error_x_dot_integrator = error_x_dot_integrator+error_x_dot(i);
   beta_ref(i) = Kp*error_x_dot(i)+Ki*error_x_dot_integrator;
   beta_ref(i) = max( min(beta_ref(i),1),-1);
   %Braking beta
   if i>2
   beta_p = beta_opt(i-1);
   beta_pp = beta_opt(i-2);
   else
       beta_p = 0;
       beta_pp = 0;
   end
   beta_vec(:,i) = three_braking_ratios(beta_p,beta_pp,0.1);
   %%%MPC Problem
   
   for j=1:n_beta
       
       x_dot_predictions(1) = x_dot(i);
       s_predict(1) = s(i);
       for k=2:(Hp+1)
           x_dot_predictions(k) = x_dot_predictions(k-1)+delta_t*(mu*g*beta_vec(j,i)-kd/m*(x_dot_predictions(k-1))^2);
           s_predict(k) = s_predict(k-1)+delta_t*x_dot_predictions(k-1);
       end
       A_tens_cm = generate_A_mtxs(CLf,CLr,CUr,I,m,x_dot_predictions,a,b,delta_t,Hp,1);
       A_tens_om = generate_A_mtxs(CUf,CUr,CLr,I,m,x_dot_predictions,a,b,delta_t,Hp2,0);
       G_cm = generate_G_dynamic(A_tens_cm,B,Hp);
       G_om = generate_G_dynamic(A_tens_om,B,Hp2);
       H_cm = generate_H_dynamic(A_tens_cm,Hp);
       H_om = generate_H_dynamic(A_tens_om,Hp2);
       K_cm = generate_G_dynamic(A_tens_cm,[0 0 -1 0 0]',Hp);
       K_om = generate_G_dynamic(A_tens_om,[0 0 -1 0 0]',Hp2);
     
       psi_r = 1/1000000*ones(length(s_predict),1);%road_curvature_generator(road_curve,lamda,s_predict)';
       E_cm = x_dot_predictions(1:Hp)'.*psi_r(1:Hp);
       E_om = E_cm(1:Hp2);
       Ref = kron(x_dot_predictions'.*psi_r,[0;1;0;0;0]);
       Eyl_cm = lw/2*ones(Hp+1,1);
       Eyl_om = lw/2*ones(Hp2+1,1);
       Eyr_cm = -lw/2*ones(Hp+1,1);
       Eyr_om = -lw/2*ones(Hp2+1,1);
       F_cm = generate_F(G_cm,Hp,x_dot_predictions,a,b);
       F_om_temp = generate_F(G_om,Hp2,x_dot_predictions,a,b);
       V_cm =  generate_V(Hp,delta_lim,delta_dot_lim,x_dot_predictions,ar_lim,af_lim,H_cm,ksi(:,i),K_cm,E_cm,Eyl_cm,Eyr_cm,a,b);
       V_om_temp =  generate_V(Hp2,delta_lim,delta_dot_lim,x_dot_predictions,ar_lim,af_lim,H_om,ksi(:,i),K_om,E_om,Eyl_om,Eyr_om,a,b);
       F_om = zeros(10*Hp2+8,Hp);
       V_om = zeros(10*Hp2+8,1);
       F_om(1:10*Hp2+8,1:Hp2) = F_om_temp;
       V_om(1:10*Hp2+8,1) = V_om_temp;
       F = [F_cm;F_om];
       V = [V_cm;V_om];
       M = G_cm'*Q_hat_cm*G_cm+R_hat_cm;
       alpha = G_cm'*Q_hat_cm'*(H_cm*ksi(:,i)+K_cm*E_cm-0*Ref);
       con = 182;
       %[u,fval] = quadprog(M,alpha);%,F,V);
       u = -M^-1*alpha;
       fval = u'*M*u+2*alpha'*u;
       delta_dot_vec(j,i) = u(1);
       f_cost_vec(j,i) = fval;
   end
   %%%Post Computation   
   [beta_temp, delta_dot_temp] = Post_Computation(beta_ref(i),beta_vec(:,i),delta_dot_vec(:,i),f_cost_vec(:,i),Qbeta);
   beta_opt(i) = beta_temp;
   delta_dot_opt(i) = delta_dot_temp;
  
   %Calculate the Torque/Forces
   %%%Simulation
   params = [m g mu a b e B_sep C_sep I psi_r(1) kd];
   y0 = [ksi(:,i);x_dot(i);s(i)];
   
   %[t,x] = ode45(@(t,x)vehicle_dynamics(t,x,params,delta_dot_opt(i),beta_opt(i)), [0 delta_t], y0);
%    ksi(:,i+1)=x(end,1:5);
%    x_dot(i+1) = x(end,6);
%    s(i+1) = x(end,7);
     A = generate_A_CM(CLf,CLr,CUr,I,m,x_dot(i),a,b,delta_t);
     ksi(:,i+1) = A*ksi(:,i)+B*delta_dot_opt(i);
     x_dot(i+1) = x_dot(i)+delta_t*g*beta_opt(i)-delta_t*kd/m*x_dot(i)^2;
     s(i+1) = s(i)+x_dot(i)*delta_t;
end

%% Plots
close all
t = linspace(0,t_dur,n_dur);
figure,
plot(t,ksi(1,1:end-1))
title('y_dot')
figure,
plot(t,ksi(2,1:end-1)*180/pi)
title('psi_dot')
figure,
plot(t,ksi(3,1:end-1)*180/pi)
title('e_psi')
figure,
plot(t,ksi(4,1:end-1))
title('e_y')
figure,
plot(t,ksi(5,1:end-1)*180/pi)
title('delta')
figure,
plot(t,x_dot(1:end-1))
hold on
plot(t,x_dot_ref*ones(1,n_dur))
title('x_dot')
figure,
plot(t,s(1:end-1))
title('s')
figure,
plot(t,beta_opt)    
title('Beta')
figure,
%y_road = road_cord_generator(road_curve,lamda,s(1:i));
[x_road , y_road] = generate_circular_road(1000000,s(1:i));
plot(x_road,y_road)
hold on
plot(x_road,y_road+lw/2)
hold on
plot(x_road,y_road-lw/2)
hold on
plot(x_road,ksi(4,1:i)+y_road)
title('Road')
legend('Center lane','left','right','car')

%% Functions
%Generate the new A matrix given the parameters
function A = generate_A_CM(Cfl,Crl,Cru,I,m,x_dot,a,b,delta_t)
    Ac = [(Cfl+Crl)/(m*x_dot) (-x_dot+(Cfl*a-Crl*b)/(m*x_dot)) 0 0 -Cfl/m;
        (a*Cfl-b*Cru)/(I-x_dot) (a^2*Cfl+b^2*Cru)/(I*x_dot) 0 0 -a*Cfl/I;
        0 1 0 0 0;1 0 x_dot 0 0; 0 0 0 0 0];
    A = eye(5)+Ac*delta_t;
end
function A = generate_A_OM(Cfu,Cru,Crl,I,m,x_dot,a,b,delta_t)
    Ac = [(Cfu+Cru)/(m*x_dot) (-m*x_dot+(Cfu*a-Cru*b)/(m*x_dot)) 0 0 -Cfu/m;
        (a*Cfu-b*Crl)/(I*x_dot) (a^2*Cfu+b^2*Crl)/(I*x_dot) 0 0 -a*Cfu/I;
        0 1 0 0 0; 1 0 x_dot 0 0; 0 0 0 0 0];
    A = eye(5)+Ac*delta_t;
end
%Generates N+1 A matricies according to given x_Dot and puts them on the A_tensor 
function A_tensor = generate_A_mtxs(Cfu,Cru,Crl,I,m,x_dot,a,b,delta_t,N,CM)
A_tensor = zeros(5,5,N+1);    
if CM==1
    for i = 1:(N+1)
    A_tensor(:,:,i) =   generate_A_CM(Cfu,Cru,Crl,I,m,x_dot(i),a,b,delta_t);
    end
else
    for i = 1:(N+1)
     A_tensor(:,:,i) =   generate_A_OM(Cfu,Cru,Crl,I,m,x_dot(i),a,b,delta_t);
    end
end
             
end


%Multiply A matricies from start to k
function AA = multiply_A(A_tensor,k)
AA = eye(5);
for i = 1:k
    AA =A_tensor(:,:,i)*AA;
end
end



%A,B are from eqns of motion, N is horizon
function G = generate_G(A,B,N)
G = zeros(size(B,1),N*size(B,2));
for j= 2:(N+1)
    temp = zeros(size(B,1),N*size(B,2));
    for i=1:N
         if j>i
             temp(:,size(B,2)*(i-1)+1:size(B,2)*(i)) = A^(j-i-1)*B;%[temp A^(j-i-1)*B];
         else
             temp(:,size(B,2)*(i-1)+1:size(B,2)*(i))=0*B;% = [temp 0*B];
         end
    end 
    G = [G;temp];
end

end
function H = generate_H(A,N)
    H = zeros((N+1)*5,5);
    for i=1:(N+1)
        H(5*(i-1)+1:5*i,:) = A^(i-1);
    end
end


%A,B are from eqns of motion, N is horizon
function G = generate_G_dynamic(A_tensor,B,N)
G = zeros(size(B,1),N*size(B,2));
for j= 2:(N+1)
    temp = zeros(size(B,1),N*size(B,2));
    for i=1:N
         if j>i
             temp(:,size(B,2)*(i-1)+1:size(B,2)*(i)) = multiply_A(A_tensor,j-i-1)*B;%[temp A^(j-i-1)*B];
         else
             temp(:,size(B,2)*(i-1)+1:size(B,2)*(i))=0*B;% = [temp 0*B];
         end
    end 

    G = [G;temp];
end

end
function H = generate_H_dynamic(A_tensor,N)
    H = zeros((N+1)*5,5);
    for i=1:(N+1)
        H(5*(i-1)+1:5*i,:) = multiply_A(A_tensor,i-1);
    end
end


% G is from the eqns of motion N is horizon, a&b are vehicle sizes, asusmes
% x is 5D
function F = generate_F(G,N,x_dot,a,b)
f_delta = [0 0 0 0 1];
ff_delta = [0 0 0 0 1];
ff_ey = [0 0 0 1 0];
for i=1:N
    ff_delta = blkdiag(ff_delta,f_delta);
    ff_ey = blkdiag(ff_ey,[0 0 0 1 0]);
end
ff_ar = [1/x_dot(1) -b/x_dot(1) 0 0 0];
ff_al = [1/x_dot(1) a/x_dot(1) 0 0 -1];
for i=2:(N+1)
    ff_ar = blkdiag(ff_ar,[1/x_dot(i) -b/x_dot(i) 0 0 0]);
    ff_al = blkdiag(ff_al,[1/x_dot(i) a/x_dot(i) 0 0 -1]);
end
F = [ff_delta*G;-ff_delta*G;eye(N);-eye(N);ff_ar*G;-ff_ar*G; ff_al*G;-ff_al*G;ff_ey*G;-ff_ey*G];
end
function V = generate_V(N,d_lim,d_dot_lim,x_dot,ar_lim,al_lim,H,x0,K,E,Eyl,Eyr,a,b)
f_delta = [0 0 0 0 1];
ff_delta = [0 0 0 0 1];
ff_ey = [0 0 0 1 0];
for i=1:N
    ff_delta = blkdiag(ff_delta,f_delta);
    ff_ey = blkdiag(ff_ey,[0 0 0 1 0]);
end
ff_ar = [1/x_dot(1) -b/x_dot(1) 0 0 0];

ff_al = [1/x_dot(1) a/x_dot(1) 0 0 -1];
for i=2:(N+1)
    ff_ar = blkdiag(ff_ar,[1/x_dot(i) -b/x_dot(i) 0 0 0]);
    ff_al = blkdiag(ff_al,[1/x_dot(i) a/x_dot(i) 0 0 -1]);
end
dd_lim = d_lim*ones(N+1,1);
dd_dot_lim = d_dot_lim*ones(N,1);
aar_lim = ar_lim*ones(N+1,1);
aal_lim = al_lim*ones(N+1,1);
T = (H*x0+K*E);

V=[dd_lim-ff_delta*T;dd_lim+ff_delta*T;dd_dot_lim;dd_dot_lim;
    aar_lim-ff_ar*T;aar_lim+ff_ar*T;aal_lim-ff_al*T;aal_lim+ff_al*T;Eyl-ff_ey*T;-Eyr+ff_ey*T];
end
function beta_vec = n_braking_ratios(n_beta,beta_ref)

    if beta_ref<0
        beta_end=0;
    else
        beta_end = beta_ref;
    end
    beta_vec = linspace(-1,beta_end,n_beta);
    beta_vec(beta_vec>1) = 1;
beta_vec(beta_vec<-1) = -1;
end

function beta_vec = three_braking_ratios(beta_opt_p,beta_opt_pp,delta_beta)
if beta_opt_p == beta_opt_pp
    beta_vec = [(beta_opt_p -delta_beta) beta_opt_p (beta_opt_p+delta_beta)];
elseif  beta_opt_p>beta_opt_pp
    beta_vec = [beta_opt_p (beta_opt_p+delta_beta) (beta_opt_p+2*delta_beta)];
else
    beta_vec = [(beta_opt_p-2*delta_beta) (beta_opt_p-delta_beta) (beta_opt_p)];
end
beta_vec(beta_vec>1) = 1;
beta_vec(beta_vec<-1) = -1;
end
function [beta, delta_dot] = Post_Computation(beta_ref,beta_vec,delta_dot_vec,f_vec,Qf)
i_best = 1;
best_cost = inf;
	
for i=1:length(beta_vec)
    current_cost = (beta_ref-beta_vec(i))*Qf*(beta_ref-beta_vec(i))'+f_vec(i);
    
    if current_cost< best_cost
        best_cost = current_cost;
        i_best = i;
    end
end
beta = beta_vec(i_best);
delta_dot = delta_dot_vec(i_best);
end
function y = road_cord_generator(amplitude,wavelength,pos_arr)
    y = amplitude.*sin(2*pi/wavelength*pos_arr);
end
%psi_r = y''*sqrt(y'^2+1)/cos^2(y')
function psi_r = road_curvature_generator(amplitude,wavelength,position)
%y = amplitude.*sin(2*pi/wavelength*position);

y_d = amplitude.*(2*pi/wavelength).*cos(2*pi/wavelength*position);
y_dd = -amplitude.*(2*pi/wavelength)^2.*sin(2*pi/wavelength*position);
psi_r = y_dd.*sqrt(y_d.^2+1)./(cos(y_d).^2);

end
%psi_r = 1/r
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