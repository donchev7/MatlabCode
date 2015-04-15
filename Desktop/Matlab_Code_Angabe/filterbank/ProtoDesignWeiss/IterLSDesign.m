function h = IterLSDesign(Lp,K,N,trade_off);
% h = IterLSDesign(Lp,K,N,trade_off);
%
% Design the prototype filter for a K channel DFT filter bank
% decimated by a factor N using an iterative least squares method.
% An initial guess is calculated using the Parks McClellan algorithm.
%
% Input parameters:
%      Lp     length of prototype (1+order) -- must be even!
%      K      number of channels
%      N      decimation ratio
%      trade_off weighting of stopband attenuation (optional)
%
% Output parameters:
%      h      coefficients of prototype filter
%
% See also: FreqGridMat(), PolyPhaseMat() 
%
% St.Weiss, University of Strathclyde, 3.9.1997

% constants
limit = 0.0001;  % stop iteration if improvement is less than 0.1 permille
overdet = 1.1;  % degree of overdeterminism 
disp('===== iterative least squares design =====');

if ~exist('trade_off'),
  trade_off = 1;
end;  

% initial guess for prototype
% f = [0 0.85/K 0.92/N 1];
f = [0 0.95/K 0.98/N 1];
a = [1 1 0 0];
h = remez(Lp-1,f,a);
epsilon = 0.000000001

% degrees of freedom / determinism
F_pc = 2*Lp/K-1;             % no. equations for power complementarity 
F_free = Lp - F_pc;          % degrees of freedom
Nw = ceil(Lp*overdet-F_pc);  % no. of frequency points
disp(sprintf('filter length         %d',Lp));
disp(sprintf('channel number        %d',K));
disp(sprintf('decimation ratio      %d',N));
disp('  ---  ');
disp(sprintf('degrees of freedom    %d',F_free));
disp(sprintf('no. frequency points  %d',Nw));
disp('  ---  ');

h = h(:)';

% optimum performance vector
p = zeros(2*Lp/K-1+Nw,1);
p(round(Lp/K)) = 1/K;

time1 = cputime;
%calc1 = flops;

% matrix for stopband critierion 
% T = FreqGridMat(Lp,Nw,0.92*pi/N);
T = trade_off*FreqGridMat(Lp,Nw,0.99*pi/N);
E = PolyPhaseMat(h,K);
A = [E; T];
b = A*h(1:end/2)';
e(1) = norm(b-p,2);
%disp(sprintf('initial error %f',e(1)));

r = 1;
i = 1;
pr = 1; snr = 1;
% iterate
while (r>limit),
  i = i+1;
  R = triu(qr([A p]));
  h_h = inv(R(1:Lp/2,1:Lp/2))*R(1:Lp/2,Lp/2+1);
  h_new = [h_h' fliplr(h_h')];
  h = (h_new+h)/2;           % Harteneck relaxation
  E = PolyPhaseMat(h,K);

  % weighting of frequency response points
%  snr = MeasureAliasDist(h,K,N); 
%  PR = CheckPR(h,K,N,'quiet');
%  pr = PR(2);
%  disp(sprintf('     pr error  %f [dB]',10*log10(pr)));
%  disp(sprintf('     snr        %f [dB]',snr));
%  snr = 10^(-snr/10);
  ratio = F_pc/Nw*snr/8/pr;
  D = diag(ratio*ones(Nw,1));

  A = [E; D*T];

  b = A*h(1:end/2)';
  e(i) = norm(b-p,2);
  %disp(sprintf('%d iteration: residual error %f',[i-1 e(i)]));
  r = abs(1 - e(i)/e(i-1));
  if (e(i) < epsilon),
    r = 0.0;
  end;  
end;   

time2 = cputime-time1;
%calc2 = flops;

% output computer related stuff
disp('  ---  ');
%disp(sprintf('elapsed time:    %f',time2));
%disp(sprintf('number of flops: %d',calc2));

% calculate final quality 
snr = MeasureAliasDist(h,K,N);
% E = CheckPR(h,K,N,'quiet');
disp('  ---  ');
% disp(sprintf('perfect reconstruction error    %f [dB]',10*log10(E(2))));
disp(sprintf('snr due to in-band aliasing      %f [dB]',snr));

