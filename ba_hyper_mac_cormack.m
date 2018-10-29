function W = ba_hyper_mac_cormack(Nz, Nt, dz, dt, w0, d)

W = zeros(size(w0,1),Nz,Nt);
w = w0;

for n = 1:Nt
    
    w_star = w(:,1:end-1) - dt/dz*(flux(w(:,2:end)) - flux(w(:,1:end-1))) + dt*source(w(:,1:end-1));    
    w(:,2:end-1) = 0.5*(w(:,2:end-1) + w_star(:,2:end)) - 0.5*dt/dz*(flux(w_star(:,2:end)) - flux(w_star(:,1:end-1))) + 0.5*dt*source(w_star(:,2:end));
    
    w(:,1) = boundary_l(w, d, n, dt);
    w(:,end) = boundary_r(w);
    
    W(:,:,n) = w;
    
    disp(['time ' num2str(n) ' / ' num2str(Nt)])
    
end

function g = flux(w)
global Gf
g = Gf*w;

function h = source(w)    
global Hs
h = Hs*w;

function w_l = boundary_l(w, d, n, dt)
global ode_sysc x x_vec
% w_l = [d*w(2,2); w(2,2)]; % no ODE
x_vec(:,n) = x; 
u = (2*w(2,2)-w(2,3))*ones(1,2);
t = 0:dt:dt;
[y,~,x_aux] = lsim(ode_sysc,u,t,x); 
x = x_aux(end,:); 
% x = Ad*x + bd*(2*w(2,2)-w(2,3)); % discrete
% global A b dt 
% x = (A*x + b*w(2,2))*dt + x; % Euler forward
w_l = [y(end,:) + d*(2*w(2,2)-w(2,3)); 2*w(2,2)-w(2,3)];

function w_r = boundary_r(w)
global flag_contr
if flag_contr == 1
    w_r = [w(1,end-1); ca_control(w)]; % coupled states, only one is influenced
else
    w_r = [2*w(1,end-1)-w(1,end-2); 0];
end

