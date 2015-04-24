function spectrum = STFTana(x,N,win)
%N         % number of filters = DFT length
L=length(x);
n=0:L-1;
h = win;  % Better DFT lowpass = Hamming window
spectrum = zeros(N,L);   % X will be the filter bank output
for k=1:N         % Loop over channels
  wk = 2*pi*(k-1)/N;
  xk = exp(-1i*wk*n).* x;  % Modulation by complex exponential
  spectrum(k,:) = filter(h,1,xk);
end

end