% DFTRealTest.m
%
% Andreas Schwarz (schwarz@lnt.de), 2011-02-01
% based on code from S Weiss, Univ. of Southampton, 21/1/2001
%
% Demonstrates a DFT filter bank for real-valued signals.
%
% Unlike in Stefan Weiss' code, no GDFT is used -> "even-channel stacking".
%
% The DFT approach seems more practical for efficient implementation. It
% only has the slight disadvantage that K/2+1 bands are required instead of
% K/2 bands for perfect reconstruction; however, it seems reasonable to
% simply discard band K/2+1 (i.e., set to 0 before synthesis).
%
% Band 1 and K/2+1 are real-valued due to symmetry at Omega=0 and Omega=pi.

% prototype filter and parameters
%p = filter16_14_448;
K = 64; N = 16; Lp = 128;
p=IterLsDesign(Lp,K,N);

% real valued input signal
Lx = 20*Lp;
x = randn(1,Lx);

% initialisation of output TDL for GDFTSynCmplx() 
y_tdl = zeros(1,Lp);

% iteration for analysis and synthesis operations
for k = 1:N:Lx,
   U = MakeTDL(x,k,Lp);											% extract TDL
   X = DFTAnaReal(U,K,p);								% analysis
   %X(1) = 0; X(K/2+1) = 0;
   [y(k:k+N-1) y_tdl] = DFTSynReal(X,N,p,y_tdl);	% synthesis
end;  

% display of I/O to cascaded analysis and synthesis filter bank
clf; 
plot(x(1:end-Lp+1)); hold on; plot(y(Lp:end),'r:'); hold off;
legend('FiBa input','delayed FiBa output');
ylabel('input, output'); title('I/O for real valued GDFT OSFB');
xlabel('time [fullband sampling periods]');


SNR = -10*log10(var(x(1:end-Lp+1)-y(Lp:length(x)))/var(x(1:end-Lp+1)));
disp(sprintf('SNR:    %f dB',SNR));