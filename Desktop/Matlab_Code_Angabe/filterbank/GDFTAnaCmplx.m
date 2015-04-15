function X = GDFTAnaCmplx(x,K,p);
% X = GDFTAnaCmplx(x,K,p);
%
% Performs a single step of a GDFT analysis filter bank operation on
% a complex vector x extracted from a tapped delay line(TDL). The length
% of x is equal to the length of the prototype filter whose (modified) 
% coefficients are contained in p. K is the number of subbands produced.
% The decimation factor N depends on how often the TDL is updated in 
% between successive calls to this function.
%
% Input Parameters:
%	x		values in TDL (newest datum in top position)
%	K		number of subbands
%	p		(modified) prototype filter
%
% Output Parameters:
%	X		column vector containing subband samples
%
% See also:
%	GDFTSynCmplex()
%
% Reference:
%   [1] S Weiss, RW Stewart: On Adaptive Filtering in Oversampled Subbands,
%       Shaker Verlag, Aachen, Germany, 1998.
%   [2] S Weiss, RW Stewart: Efficient Implementation of Oversampled
%       Modulated Filter Banks, Electronics Letters, August 2000.
%
% S. Weiss, Univ of Southampton, 21/1/2001

% parameters for GDFT modulations offsets
k0 = 0.5;
n0 = -(K-1)/2;
Lp = length(p);

% multiply by prototype coeffs and accumulate into K components
U = zeros(K,ceil(Lp/K));
U(1:Lp) = x.*p;
V = sum(U.').';

% perform rotation by modulation matrix
d1 = exp(sqrt(-1)*2*pi/K*n0*(0:K-1)');
d2 = exp(sqrt(-1)*2*pi/K*k0*((0:K-1)'+n0));
X = d1.*ifft(d2.*V);


