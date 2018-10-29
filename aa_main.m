clear all
clc
% close all

global flag_contr Gf Hs bp Ad bd N Nz Nt x A b c z dz dt Kvu Kvv Kuu Kvu1 Kvv1 Pvv1 n10 n20 n21 c1 c4 lam1 lam2 x_vec ode_sysc

% control on flag_contr = 1: on, 0: off
flag_contr = 1;
scaleval = 1;

% PDE
lam1 = 1; lam2 = 2; c1 = -1; c2 = -1; c3 = -2; c4 = -2; d = 0.5;
% lam1 = 1; lam2 = 2; c1 = -1; c2 = -1; c3 = -2; c4 = -2; d = -0.5;
% cpar = 7.5e1;
% lam1 = 1; lam2 = 2; c1 = -cpar; c2 = -cpar; c3 = -cpar; c4 = -cpar; d = 0.5;
Gf = - [-lam1, 0; 0, lam2]*scaleval; % flux (change sign for MacCormack)
Hs = [c1, c2; c3, c4]*scaleval; % sink/source
bp = d; % boundary parameters

% ODE
N = 1; % dimension
A = 1*scaleval; % dynamic matrix
b = 1*scaleval; % input vector
c = 1; % output vector
lamD = -0.1; % desired eigenvalues
% N = 2; % dimension
% A = [1 2; 4 2]; % dynamic matrix
% b = [1; 2]; % input vector
% c = [0; 1]; % output vector
% lamD = [-5 -4]; % desired eigenvalues

% simulation parameters
dz = 0.01; 
dt = 0.99*dz/max(abs(Gf(:)));  
z = 0:dz:1;
xi = z;
t = 0:dt:5/scaleval;
Nz = length(z);
Nt = length(t);

% initial conditions
u0 = max((1-20*(z-.3).^2),0).^2;
v0 = 0.5*u0;
w0 = [u0; v0]; 
x = zeros(N,1);
x(1) = 1;
if N>1
    x(end) =  (u0(1) - d*v0(1) - c(1:end-1)'*x(1:end-1))/c(end); 
else
    x(end) =  (u0(1) - d*v0(1) - c(1)'*x(1))/c(end); 
end

% convert to discrete state space
ode_sysc = ss(A,b,c',0);
ode_sysd = c2d(ode_sysc,dt);
% Ad = ode_sysd.A;
% bd = ode_sysd.B;
x_vec = zeros(length(x),Nt);

% controller of ODE
n20 = -place(A, b, lamD)'; % continuous
% n20 = -place(Ad, bd, exp(lamD*dt))'; % discrete

% calculate kernels
Kvu = zeros(length(z));
Kvv = zeros(length(z));
for I = 1:length(z)
    Kvu(I,1:I) = da_kernel_Kvu(z(I), xi(1:I), d, lam1, lam2, c1, c2, c3, c4);
    Kvv(I,1:I) = db_kernel_Kvv(z(I), xi(1:I), d, lam1, lam2, c1, c2, c3, c4);
end

Kvu1 = Kvu(end,:);
Kvv1 = Kvv(end,:);
[n2, Pvv] = dc_n_Pvv(lam1, lam2);
n21 = n2(:,end);
Pvv1 = 1/lam2*b'*fliplr(n2);

% simulation
%W = ba_hyper_mac_cormack(Nz, Nt, dz, dt, w0, d);
W = bb_hyper_lax_friedrich(Nz, Nt, dz, dt, w0, d);
W = real(W);
x_vec = real(x_vec);

x1 = squeeze(W(1,:,:));
x2 = squeeze(W(2,:,:));
x1t = numdiff(x1, 2, dt);
x1z = numdiff(x1,  1, dz);

% figure
% surf(x1)
% return
% 
% auxterm2 = (x1t + eps1*x1z - c1*x1 - c2*x2);
% auxterm2(isnan(auxterm2)) = 0;
% figure
% 
% surf(auxterm2(1:end,1:end))
% return

% plot
figure
for nn = 1:length(t)
    subplot(2,1,1)
    plot(z,w0,'b:',z,W(1,:,nn),'k.-',z,W(2,:,nn),'r.-')
%     plot(z,w0,'b:',z,W(2,:,nn),'r.-')
    axis([min(z) max(z) min(min(min(W))) max(max(max(W)))])
%     axis([min(z) max(z) min(min(min(W(2,:,:)))) max(max(max(W(2,:,:))))])
    title(['t = ' num2str(t(nn),'%.2f')])
    subplot(2,1,2)
    plot(t,x_vec, [t(nn) t(nn)], [min(min(x_vec)) max(max(x_vec))], 'k')
    axis([min(t) max(t) min(min(x_vec)) max(max(x_vec))])
    drawnow
%     pause(0.5)
end
%return
%% calculate v_til
% Kuv = zeros(length(z));
% Kuu = zeros(length(z));
% for I = 1:length(z)
%     Kuv(I,1:I) = da_kernel_Kvu(z(I), xi(1:I), 1/d, lam2, lam1, -c4, -c3, -c2, -c1);
%     Kuu(I,1:I) = db_kernel_Kvv(z(I), xi(1:I), 1/d, lam2, lam1, -c4, -c3, -c2, -c1);
% end
% 
% n10 = c + d*n20;
% [n1, Puu] = dc_n_Puu(lam1);
% 
% u = exp(-c1/lam1*z'*ones(1,Nt)).*squeeze(W(1,:,:));
% v = exp(c4/lam2*z'*ones(1,Nt)).*squeeze(W(2,:,:));
% 
% u_bar = zeros(Nz,Nt);
% v_bar = zeros(Nz,Nt);
% for nn = 1:Nt
%     for mm = 1:Nz
%         u_bar(mm,nn) = u(mm,nn) - dz*trapz(Kuu(mm,1:mm).*u(1:mm,nn)') - dz*trapz(Kuv(mm,1:mm).*v(1:mm,nn)'); 
%         v_bar(mm,nn) = v(mm,nn) - dz*trapz(Kvu(mm,1:mm).*u(1:mm,nn)') - dz*trapz(Kvv(mm,1:mm).*v(1:mm,nn)');
%     end
% end
% 
% u_til = zeros(Nz,Nt);
% v_til = zeros(Nz,Nt);
% for nn = 1:Nt
%     for mm = 1:Nz
%         u_til(mm,nn) = u_bar(mm,nn) - dz*trapz(Puu(mm,1:mm).*u_bar(1:mm,nn)') - n1(:,mm)'*x_vec(:,nn);
%         v_til(mm,nn) = v_bar(mm,nn) - dz*trapz(Pvv(mm,1:mm).*v_bar(1:mm,nn)') - n2(:,mm)'*x_vec(:,nn);
%     end
% end

% u_til = real(u_til);
% v_til = real(v_til);
% 
% % plot
% figure
% for nn = 1:length(t)
%     subplot(2,1,1)
%     plot(z,w0,'b:',z,u_til(:,nn),'k.-',z,v_til(:,nn),'r.-')
%     axis([min(z) max(z) min(min(min([u_til, v_til]))) max(max(max([u_til, v_til])))])
%     if t(nn) >= 1/lam1 + 1/lam2
%         title(['t = ' num2str(t(nn),'%.2f') '\color{red} min cont time = ' num2str(1/lam1 + 1/lam2)])
%     else
%         title(['t = ' num2str(t(nn),'%.2f') ' min cont time = ' num2str(1/lam1 + 1/lam2)]) 
%     end
%     subplot(2,1,2)
%     plot(t,x_vec, [t(nn) t(nn)], [min(min(x_vec)) max(max(x_vec))], 'k', [1/lam1 + 1/lam2 1/lam1 + 1/lam2], [min(min(x_vec)) max(max(x_vec))], 'r')
%     axis([min(t) max(t) min(min(x_vec)) max(max(x_vec))])
%     drawnow
% %     pause(0.5)
% end

