N=5;
A = 5*eye(5);
B = [0 0 0 0 1]';
G = generate_G(A,B,N);
x_dot = [3 1 2 3 4 5]';
a = 2;
b=3;
H = generate_H(A,N);
beta_ref = 0.7;
beta_vec = three_braking_ratios(0.5,0.3,0.1);
f_arr = [1 1.1 1.1];
Qf = 3;
delta_dot_vec = [1 2 3];
[beta,delta_dot] = Post_Computation(beta_ref,beta_vec,delta_dot_vec,f_arr,Qf) ;

x = linspace(0,20,1000);
%A = 0.1;
lamda = 90;
A = 0.5;%*cos(2*pi/50*x);
rad = 10000; 
s = linspace(0,100,1000);
[x,y] = generate_circular_road(rad,s);
figure
% y = road_cord_generator(A,lamda,x);
figure(1),
plot(x,y,'LineWidth',2)
hold on
plot(x,y+2.3,'LineWidth',2)
hold on
plot(x,y-2.3,'LineWidth',2)
% psi_r = road_curvature_generator(A,lamda,x);
% % psi_r = zeros(1,1000);
% % for i=1:1000
% %     psi_r(i) = road_curvature_generator(0.1,4,x(i));
% % end
% figure(2),
% plot(x,psi_r)

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
    disp(temp)
    disp(G)
    G = [G;temp];
end

end
function H = generate_H(A,N)
    H = zeros((N+1)*5,5);
    for i=1:(N+1)
        H(5*(i-1)+1:5*i,:) = A^(i-1);
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
function V = generate_V(N,d_lim,d_dot_lim,ar_lim,al_lim,H,x0,K,E,Eyl,Eyr)
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
V=[dd_lim-ff_delta*T;dd_lim+ff_delta*T;dd_dot_lim;-dd_dot_lim;
    aar_lim-ff_ar*T;aar_lim+ff_ar*T;aal_lim-ff_al*T;aal_lim+ff_al*T;Eyl-ff_ey*T;Eyr+sff_ey*T];
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
    disp(current_cost)
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
disp(size(amplitude))
disp(size(position))
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