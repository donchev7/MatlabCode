function spectrum = STFTana(x)
N=512;           % number of filters = DFT length
fs=16000;        % sampling frequency (arbitrary)
D=10;            % duration in seconds
L = length(x);
L2 = ceil(fs*D)+1; % signal duration (samples)
n = 0:L2-1;        % discrete-time axis (samples)
t = n/fs;         % discrete-time axis (sec)
%x = chirp(t,0,D,fs/2);   % sine sweep from 0 Hz to fs/2 Hz
%x = echirp(t,0,D,fs/2); % for complex "analytic" chirp
x = x(1:L2);       % trim trailing zeros at end
%h = ones(1,N);    % Simple DFT lowpass = rectangular window
h = hamming(N);  % Better DFT lowpass = Hamming window
spectrum = zeros(N,L2);   % X will be the filter bank output
for k=1:N         % Loop over channels
  wk = 2*pi*(k-1)/N;
  xk = exp(-1i*wk*n);
  xs = xk.* x.';    % Modulation by complex exponential
  spectrum(k,:) = filter(h,1,xs);
end
% N=512;
% nframes = length(x)/R;
% %spectrum = zeros(N,nframes); % pre-allocate STFT output array
% M = length(win);           % M = window length, N = FFT length
% zp = zeros(N-M,1);       % zero padding (to be inserted)
% xoff = 0;                % current offset in input signal x
% Mo2 = (M)/2;           % Assume M even for simplicity here
% for m=1:nframes
%   xt = x(xoff+1:xoff+M); % extract frame of input data
%   xtw = win .* xt;         % apply window to current frame
%   xtwz = [xtw(Mo2+1:M); zp; xtw(1:Mo2)]; % windowed, zero padded
%   spectrum(:,m) = fft(xtwz); % STFT for frame m
%   xoff = xoff + R;       % advance in-pointer by hop-size R
% end
end