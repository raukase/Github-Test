function [n1, Puu] = dc_n_Puu(lam1) % (12) Nicoles Unterlagen
global n10 A b c dz Kuu bp Nz N
d = bp;

b1 = -lam1*c*Kuu(:,1)';

% Compute n2
f = -1/lam1*A'*n10*ones(1,Nz) + 1/(lam1*d)*b'*n10*c*ones(1,Nz) + 1/lam1*b1 + 1/(lam1^2*d)*dz*cumtrapz(b1,2)*(b'*n10); 

Aaux = A';
CBaux = c*b';
G = zeros(Nz, Nz, N^2); 
for kk = 1:N^2
    for m = 1:Nz 
        for n = 1:Nz 
            G(m,n,kk)= -1/lam1*Aaux(kk) + 1/(d*lam1)*CBaux(kk) +  1/(lam1^2*d)*dz*trapz(b1(ceil(kk/N),1:m-n+1)*b(ceil(kk/N))); 
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

n1 = n10*ones(1,Nz) + dz*cumtrapz(y,2);

Puu = zeros(Nz);
for nn = 1:Nz
    Puu(nn,1:nn) = - 1/(d*lam1)*b'*fliplr(n1(:,1:nn));
end

