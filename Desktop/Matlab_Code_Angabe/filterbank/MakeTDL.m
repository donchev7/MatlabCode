function X = MakeTDL(x,k,N);
% X = MakeTDL(x,k,N);
%
% Extract a tap delay line of length N from a signal x at time k in
% reverse order:
%    X = [x(k) x(k-1) x(k-2) ... x(k-N+1)];
% For transients or negative indices, zero padding takes effect.
% 
% If x is a matrix, them TDLs are extracted with rows as temporal
% dimension.
%
% Input parameters:
%       x      input signal (row vector)
%       k      time index
%       N      TDL length
%
% Output parameters: 
%       X      samples in TDL
%
% S Weiss, 16/10/2000

[m,n] = size(x);

if k >= N,
   X = x(:,k:-1:k-N+1);
elseif k <= 0,
   X = zeros(m,N);
else
   X = [x(:,k:-1:1) zeros(m,N-k)];
end;   