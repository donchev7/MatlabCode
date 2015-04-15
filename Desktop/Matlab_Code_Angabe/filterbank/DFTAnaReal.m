function X = DFTAnaReal(x,K,p);

Lp = length(p);

% multiply by prototype coeffs and accumulate into K components
U = zeros(K,ceil(Lp/K));
U(1:Lp) = x.*p;
V = sum(U.').';

X = fft(V);
X = X(1:K/2+1);