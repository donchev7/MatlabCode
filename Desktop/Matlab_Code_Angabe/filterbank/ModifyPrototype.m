function pmod = ModifyPrototype(p,K);
% pmod = ModifyPrototype(p,K);
%
% This function flips the sign of blocks of K coefficients in the
% row vector p. This is here used to modify a prototype filter for
% a GDFT modulated filter bank in order to yield a period of K for
% the modulating transform.
%
% Input parameters:
%    p    prototype filter coefficients
%    K    block length of modification (= number of subbands)
%
% Output parameters:
%    pmod modified prorotype filter
%
% S. Weiss, Univ of Southampton, 21/1/2001

% modify prototype filter
Lp = length(p);
I = [ceil(Lp/K), K];
c = [ones(I) -ones(I)]';
pmod = c(1:Lp).*p;


