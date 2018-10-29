function [n2, Pvv] = dc_n_Pvv(lam1, lam2) % (12) Nicoles Unterlagen
global n20 A b c dz Kvu Nz N

b2 = -lam1*c*Kvu(:,1)';

% Compute n2
f = 1/lam2*A'*n20*ones(1,Nz) - 1/lam2*b2 + 1/lam2^2*dz*cumtrapz(b2,2)*(b'*n20); 

Aaux = A';
G = zeros(Nz, N^2); 
for kk = 1:N^2
    for m = 1:Nz 
        for n = 1:Nz 
            G(m,n,kk) = 1/lam2*Aaux(kk) + 1/lam2^2*dz*trapz(b2(ceil(kk/N),1:m-n+1)*b(ceil(kk/N))); 
        end 
    end
end

y = zeros(N, Nz);
for kk = 1:N
    for nn = 1:50
        y_aux = y;
        for n = 1:Nz 
            Gaux = 0;
            for jj = 1:N
                Gaux = Gaux + G(n,1:n,(kk-1)*N+jj).*y_aux(jj,1:n);
            end
            y(kk,n) = f(kk,n) + dz*trapz(Gaux); 
        end
        disp(num2str(norm(y - y_aux)))    
    end
end

n2 = n20*ones(1,Nz) + dz*cumtrapz(y,2);

Pvv = zeros(Nz);
for nn = 1:Nz
    Pvv(nn,1:nn) = 1/lam2*b'*fliplr(n2(:,1:nn));
end
n=2;




