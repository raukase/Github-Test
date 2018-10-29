function U = ca_control(w)
global x z dz Kvu Kvv Pvv1 Kvu1 Kvv1 n21 c1 c4 lam1 lam2

w(1,:) = exp(-c1/lam1*z).*w(1,:);  % Transformation aus Gehring, für antidiagonale PDE
w(2,:) = exp(c4/lam2*z).*w(2,:);

Kvu_u_int = zeros(1,length(z));     % Preallocate Kernels
Kvv_v_int = zeros(1,length(z));
for nn = 1:length(z)
    Kvu_u_int(nn) = dz*trapz(Kvu(nn,1:nn).*w(1,1:nn));
    Kvv_v_int(nn) = dz*trapz(Kvv(nn,1:nn).*w(2,1:nn));
end

U = dz*trapz(Kvu1.*w(1,:)) + dz*trapz(Kvv1.*w(2,:)) ...
    + dz*trapz(Pvv1.*w(2,:)) - dz*trapz(Pvv1.*Kvu_u_int) ...
    - dz*trapz(Pvv1.*Kvv_v_int) + n21'*x;

U = exp(-c4/lam2)*U;




