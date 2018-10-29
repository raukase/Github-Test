function val = P(x,y)
val_mod = zeros(1,length(x));
for n = 1:20
    val_mod = val_mod + x.^(n-1).*exp(-x)/factorial(n-1).*gammainc(y, 1 + n-1,'lower');
end
Q_mod = 1 - val_mod;
val = exp(x+y).*Q_mod;