function X = DFTAnaReal_TDL(x_in,k,K,p)

Lp = length(p);

% integrated MakeTDL for performance
[m,n] = size(x_in);

if k >= Lp,
   x = x_in(:,k:-1:k-Lp+1);
elseif k <= 0,
   x = zeros(m,Lp);
else
   x = [x_in(:,k:-1:1) zeros(m,Lp-k)];
end

% multiply by prototype coeffs and accumulate into K components
U = zeros(K,ceil(Lp/K));
U(1:Lp) = x.*p;
V = sum(U.').';

X = fft(V);
X = X(1:K/2+1);