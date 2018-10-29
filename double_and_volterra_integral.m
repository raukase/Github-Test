function [] = double_and_volterra_integral()
% doubleint1()
% doubleint2()
% sucapprox()

function [] = doubleint1()
% Doubleintegration 1 -> see Maple
clc
close all
clear
% f = sin(xi)*z
% g = z^2
% h = 1/xi
% int(h(xi)*int(f(xi, tau)*g(tau), tau=0..xi), xi=0..z)
dz = 0.1;
z = 0.1:dz:1.1;
xi = z;
f = zeros(length(z));
for n = 1:length(z)
    f(n,:) = sin(xi).*z(n);
end
g = z.^2;
h = 1./xi;
fgint = zeros(1,length(z));
for n = 1:length(z)
    fgint(n) = dz*trapz(f(n,1:n).*g(1:n));
end
lsg_num = dz*cumtrapz(h.*fgint);
% keyboard
lsg_ana = -sin(z).*z.^2 + 6*sin(z) - 4*cos(z).*z - 2*z;
plot(z,lsg_ana,z, lsg_num)


function [] = doubleint2()
% Doubleintegration 1 -> see Maple
% y(z) = int(int(f(xi), xi=0..z-sigma)*h(sigma), sigma=0..z) 
clc
close all
clear all

dz = 0.01; 
z = 0:dz:1;
f = z.^2.*exp(-z); 
h = sin(z); 
% f = @(z) z.^2.*exp(-z); 
% h = @(z) sin(z); 

y = zeros(size(z)); 
% for I = 1:numel(z) 
%     tau = @(sigma) z(I) - sigma;
%     G =  @(sigma) integral(f,0,tau(sigma));
%     intfun = @(z) G(z)*h(z);
%     y(I) = integral(intfun, 0, z(I), 'ArrayValued',true); 
% end

% for I = 1:numel(z) 
%     y(I) = integral2(@(xi, sigma) f(xi) .* h(sigma), 0, z(I), 0, @(xi) z(I)-xi); 
% end
 
G = zeros(length(z)); 
for m = 1:length(z) 
    for n = 1:length(z) 
        G(m,n) = dz*trapz(f(1:m-n+1)); 
    end 
end 
for n = 1:length(z) 
    y(n) = dz*trapz(G(n,1:n).*h(1:n)); 
end 

y_ana = 1/2*(exp(z).*cos(z) - exp(z).*sin(z) - z.^2 + 4*exp(z) - 4*z - 5).*exp(-z); 
plot(z,y,z,y_ana)

function [] = sucapprox()
% Method of successive approximation -> see Nicole notes Volterra integral example
clc
close all
clear
a = 2;
c = 0.33;
x0 = 52;
dz = 0.01;
z = (0:dz:1)';
xi = z;
% T = toeplitz(z, zeros(size(xi)));
% A1 = zeros(size(z));
% for n = 1:100
%     A1 = exp(a*z)*x0 + c*dz*trapz(exp(a*T).*tril(ones(length(z))).*repmat(A1', length(A1),1),2);
% end
A = zeros(size(z));
for n = 1:5
    A_aux = A; % dadurch wird in einer Iteration nur das alte A verwendet
    for m = 1:length(z)
        A(m) = exp(a*z(m))*x0 + c*dz*trapz(exp(a*(z(m) - z(1:m))).*A_aux(1:m));
    end
end
plot(z,exp((a + c)*z)*x0)
hold on
plot(z,A)


