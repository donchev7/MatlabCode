% GDFTCmplxTest.m
%
% demonstrates the use of GDFTAnaCmplx() and GDFTSynCmplx()
% in combination with ModifyPrototype() and MakeTDL().
%
% S Weiss, Univ. of Southampton, 21/1/2001

% prototype filter and parameters
p = filter16_14_448;
K = 16; N = 14; Lp = 448;
p_mod = ModifyPrototype(p,K);

% input signal
Lx = 2*Lp;
x = randn(1,Lx) + sqrt(-1)*randn(1,Lx);

% initialisation of output TDL for GDFTSynCmplx() 
y_tdl = zeros(1,Lp);

% iteration for analysis and synthesis operations
for k = 1:N:Lx,
   U = MakeTDL(x,k,Lp);											% extract TDL
   X = GDFTAnaCmplx(U,K,p_mod);								% analysis
   [y(k:k+N-1) y_tdl] = GDFTSynCmplx(X,N,p_mod,y_tdl);	% synthesis
end;  

% display of I/O to cascaded analysis and synthesis filter bank
clf; subplot(211);
plot(real(x(1:end-Lp+1))); hold on; plot(real(y(Lp:end)),'r:'); hold off;
legend('FiBa input','delayed FiBa output');
ylabel('real(.)'); title('I/O for complex valued GDFT OSFB');
subplot(212);
plot(imag(x(1:end-Lp+1))); hold on; plot(imag(y(Lp:end)),'r:'); hold off;
ylabel('imag(.)'); xlabel('time [fullband sampling periods]');
