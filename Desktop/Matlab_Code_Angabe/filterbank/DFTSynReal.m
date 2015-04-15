function [x_hat,yout] = DFTSynReal(Y,N,p,yin)

K = 2*(length(Y)-1);
Lp = length(p);

Y = real(ifft([Y; conj(Y(end-1:-1:2))]));

% multiply with prototype coeffs and accumulate onto updated TDL
V = Y*ones(1,ceil(Lp/K));
yout = [zeros(1,N) yin(1:Lp-N)] + V(1:Lp).*p;
x_hat = (K*N)*yout(Lp:-1:Lp-N+1);

