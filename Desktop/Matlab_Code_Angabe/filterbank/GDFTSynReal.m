function [x_hat,yout] = GDFTSynReal(Y,N,p,yin);
% [x_hat,yout] = GDFTSynReal(Y,N,p,yin);
%
% Performs a single step of a GDFT synthesis filter bank operation on
% a complex vector Y containing K/2 subband samples, whereby it is assumed 
% that the underlying fullband signal is real valued. The output from the
% synthesis bank is latched into a tapped delay line (TDL), which is shifted
% by N samples with every call of this function. The initial values in the
% TDL must be submitted as yin, after execution of a synthesis step, the
% values in the updated TDL are returned in yout. The N samples in x_hat
% are the fullband output samples of the filter bank (newest datum in last
% position).
%
% Input Parameters:
%	Y		column vector with K/2 subband samples
%	N		decimation factor
%	p		prototype filter
%	yin		current state of output TDL
%
% Output Parameters:
%	x_hat	N reconstructed output samples
%	yout	updated TDL of state values
%
% See also:
%	GDFAnaReal(), GDFTAnaCmplex(), GDFTSynCmplx()
%
% Reference:
%   [1] S Weiss, RW Stewart: On Adaptive Filtering in Oversampled Subbands,
%       Shaker Verlag, Aachen, Germany, 1998.
%   [2] S Weiss, RW Stewart: Efficient Implementation of Oversampled
%       Modulated Filter Banks, Electronics Letters, August 2000.
%
% S. Weiss, Univ of Southampton, 21/1/2001

% parameters for GDFT modulation offsets
K = 2*length(Y);
k0 = 0.5;
n0 = -(K-1)/2;
Lp = length(p);

% derotation by modulating transform
d1 = exp(-sqrt(-1)*2*pi/K*n0*(0:K/2-1)');
d2_r = cos(2*pi/K*k0*((0:K-1)'+n0));
d2_i = sin(2*pi/K*k0*((0:K-1)'+n0));
Y = fft([d1.*Y; zeros(K/2,1)]);
Y = d2_r.*real(Y) + d2_i.*imag(Y);

% multiply with prototype coeffs and accumulate onto updated TDL
V = Y*ones(1,ceil(Lp/K));
yout = [zeros(1,N) yin(1:Lp-N)] + V(1:Lp).*p;
x_hat = (2*K*N)*yout(Lp:-1:Lp-N+1);

