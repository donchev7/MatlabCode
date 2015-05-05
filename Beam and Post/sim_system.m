function [y_sv,varargout] = ...
sim_system(speech,noise,phi_d,beam_nr,mue,filter_nr,cn,aa,alpha,mics,N,L,fs,flow,fhigh)

% [y_sv,y_s,y_v,y_beam,y_bs,y_bv] = ...
% sim_system(speech,noise,phi_d,beam_nr,mue,filter_nr,cn,aa,alpha,mics,N,L,fs,flow,fhigh)
% Complete Simulation System for Post-Filtering techniques with
% MVDR-Beamformers in a diffuse noisefield (Master-Slave-Algorithm)
% y_sv
% y_s
% y_v
% y_beam
% y_bs
% y_bv processed noisy input
% processed speech-only signal (optional)
% processed noise-only signal (optional)
% beamformer output (optional)
% speech-only beamformer output(optional)
% noise-only signal beamformer output (optional)
% speech
% noise clean input signal matrix
% input noise matrix
% hint: if you don't need the Slave-Algorithm please declare
% speech = noisy_input and noise = 0
% desired azimuth angle (to direction of arrival)
% you can choose between the following beamformers
% 'DSB' ... Delay&Sum-Beamformer
% 'SDB' ... Superdirective Beamformer
% default beam_nr = 'SDB'
% for the regularization of Coherence Matrix Gamma
% (see thesis, section 3.4)
% mue in dB (typ. between -10 and -40dB); mue = 0 ... uses zero mue;
% default mue = -20
% name of the chosen filter
% 'ZEL88' ...... Zelinski Filter based on Welch-estimated spectral
% density functions
% 'ZEL88p'...... Zelinski Filter based on Welch-estimation plus
% post-processing method (see Zelinski 1988)
% 'SIM92' ...... Simmer 1992
% 'APAB' ....... Adaptive Postfilter with Arbitrary Beamformer
% (see book "Micorphone Arrays" by Brandstein)
% 'APES' ....... Adaptive Postfilter Extension for
% Superdirective beamformers (see Bitzer et al. 1999)
% 'MCCC' ....... McCowan Filter (see McCowan 2003)
% 'GMCC' ....... McCowan Filter, correct solution for all angles (see
% diploma thesis, section 4.3)
% default filter_nr = 'ZEL88'
% introduce Comfort Noise to reduce speech distortions
% (Minimum-Filter, see thesis section 4.6.1)
% 0 ..... no Comfort Noise (default)
% 1 ..... activate Comfort Noise
% use adaptive Welch-parameter
% (see thesis, section 4.6.2)
% 0 ..... use fixed Welch-parameter (default)
% 1 ..... use adaptive Welch-parameter
% factor for Welch estimation; default alpha = 0.8
% phi_d
% beam_nr
% mue
% filter_nr
% cn
% aa
% alpha
% mics
% N
% L
% fs
% flow
% fhigh
% set this value also for an adaptive Welch-parameter
% microphone positon matrix; default load mics_8xh.mat
% FFT length; default N = 512
% decimation factor L = N/M; default L = 4
% sampling frequency in Hz; default fs = 16000 (16kHz)
% lowest frequency in Hz; default 200 Hz
% highest frequency in Hz; default 6800 Hz
% used functions: mvdr.m, post_filter.m, welch_est.m
% Syntax:
% [y_sv] = sim_system(noisy_signal,0,...)
% [y_sv] = sim_system(speech,noise,...)
% [y_sv,y_beam] = sim_system(noisy_signal,0,...)
% [y_sv,y_s,y_v] = sim_system(speech,noise,...)
% [y_sv,y_s,y_v,y_beam] = sim_system(speech,noise,...)
% [y_sv,y_s,y_v,y_beam,y_bs,y_bv] = sim_system(speech,noise,...)
if nargin<15 fhigh = 6800; end
if nargin<14 flow = 200; end
if nargin<13 fs = 16000; end
if nargin<12 L = 4; end
if nargin<11 N = 512; end
if nargin<10 load mics_8xh.mat; end
if nargin<9 alpha = 0.8; end
if nargin<8 aa = 0; end
if nargin<7 cn = 0; end
if nargin<6 filter_nr = 'ZEL88'; end;
if nargin<5 mue = -20; end;
if nargin<4 beam_nr = 'SDB'; end
if nargin<3
    help sim_system
    return;
end
% Parameters
theta_d = 56.4303;
% elvation angle to direction of arrival (fixed)
fact = 1.3;
% additional parameter to scale the output signals
[K,Dim] = size(mics);
[Nx_orig,K_signal] = size(speech);
if K_signal ~= K
error('number of mics does not match number of speechsignals')
end
% *************************************************************************
% Zero-padding to reach a signallength to be a multiple of L
n_orig = [1:Nx_orig];
n_zeros = ceil(Nx_orig/N)*N - Nx_orig;
if noise == 0
    solo = true;
    signal = speech;
else
    signal = speech + noise;
    solo = false;
    speech = [speech;zeros(n_zeros,K)];
    noise = [noise;zeros(n_zeros,K)];
end
signal = [signal;zeros(n_zeros,K)];
Nx = length(signal);
% ************************************************************************
% Initialise Vectors and Matrices, and calculation different parameters
% init Hanning-Window for FFT
h = hanning(N);
H = h(:) * ones(1,K);
% Calculate decimationfactor and half of FFT-Length
M = N/L;
N2 = N/2 + 1;
n2 = 1:(N2);
% Calculate mue
if mue ~= 0
    mue = 10^(mue/10);
end
% Initialise output vectors
y_sv = zeros(Nx,1);
switch nargout
    case 1
    case 2
        if solo
            y_beam = zeros(Nx,1);
        else
            error('Noise input has to be Zero for this choice of output variables!!')
        end
    case 3
        if ~solo
            y_s = zeros(Nx,1);
            y_v = zeros(Nx,1);
        else
            error('Noise input has to be Zero for this choice of output variables')
        end
    case 4
        if ~solo
            y_s = zeros(Nx,1);
            y_v = zeros(Nx,1);
            y_beam = zeros(Nx,1);
        else
            error('Noise input has to be Zero for')
        end
    case 6
        if ~solo
            y_s = zeros(Nx,1);
            y_v = zeros(Nx,1);
            y_beam = zeros(Nx,1);
            y_bs = zeros(Nx,1);
            y_bv = zeros(Nx,1);
        else
            error('Noise input has to be Zero for')
        end
    otherwise
        error('This choice of output variables is')
end
% Calculate Lowpass and Highpass Filter
[bh,ah] = butter(4,flow/(fs/2),'high');
h_high = freqz(bh,ah,N2,fs);
[bl,al] = butter(4,fhigh/(fs/2),'low');
h_low = freqz(bl,al,N2,fs);
% Calculate number of all possible sensor combinations
comb = 2/(K*(K-1));
% Initialize vectors necessary to calculate the different postfilters
switch filter_nr
    case 'SIM92'
        cross_dens = zeros(N2,1/comb);
        Pyy = zeros(N2,1);
    case {'ZEL88','ZEL88p'}
        auto_dens = zeros(N2,K);
        cross_dens = zeros(N2,1/comb);
    case 'MCC03'
        auto_dens = zeros(N2,K);
        cross_dens = zeros(N2,1/comb);
    case 'GMCC'
        auto_dens = zeros(N2,K);
        cross_dens = zeros(N2,1/comb);
    case 'APAB'
        Pxx = zeros(N2,1);
        Pyy = zeros(N2,1);
    case 'APES'
        Pyy = zeros(N2,1);
        Pxx = zeros(N2,K);
        Pzz = zeros(N2,1);
    end
% Check microphone position matrix
if (Dim < 2) | (K < 1)
    error('bad matrix of microphine coordinates');
end
if Dim == 2
    rn = [mics zeros(K,1)];
else
    rn = mics;
end
%Calc. of angles in rad
theta_d = theta_d(:).' * pi / 180;
phi_d = phi_d(:).' * pi / 180;
% Calculate time alignment vector
ed = [sin(theta_d).*cos(phi_d); sin(theta_d).*sin(phi_d); cos(theta_d)];
tau = rn*ed;
% Define Matrix of the micorphone distances
xc = rn(:,1);
xc = xc(:,ones(K,1));
dxc = xc - xc.';
yc = rn(:,2);
yc = yc(:,ones(K,1));
dyc = yc - yc.';
if Dim == 2
    dR = sqrt(dxc.^2 + dyc.^2);
else
    zc = rn(:,3);
    zc = zc(:,ones(K,1));
    dzc = zc - zc.';
    dR = sqrt(dxc.^2 + dyc.^2 + dzc.^2);
end
% Calculate Coherencematrix, MVDR-Beamformer coefficients and the steering
% vector
for l = 1:N2
    beta = 2*(l-1)/N*fs/340;
    % Coherence function for a diffuse noise field
    Gamma_dum = sinc(beta*dR);
    % use constrained design for a SDB
    if strcmp(beam_nr,'SDB')
        Gamma_const = tril(Gamma_dum,-1)./(1+mue) + diag(diag(Gamma_dum)) + ...
        triu(Gamma_dum,1)./(1+mue);
    else
        Gamma_const = Gamma_dum;
    end
    Gamma(:,:,l) = Gamma_const;
    % Calculate MVDR-Beamformer Coefficients
    [W(:,l),d0(:,l)] = mvdr(tau,Gamma(:,:,l),(l-1)/N*fs,beam_nr);
end

% Change Coherencematrix to a matrix of all sensor combinations
% needed for McCowan 2003 postfilter calculation
counter = 0;
for m = 1:(K-1)
    for n = (m+1):K
        counter = counter + 1;
        Gamma_comb(:,counter) = real(Gamma(m,n,:));
    end
end
% ************************************************************************
% Main program, Master-Slave-Algorithm
for k = 1:M:(Nx-N+1)
    k1 = k:(k+N-1);
    % ********************************************************************
    % proceeded Master-Algorithm and determination of Post-Filter bins
    % FFT - Filterbank with Hanning-Windowing
    X = fft(signal(k1,:) .* H,N).';
    % Multiplicate Beamformer weighting coefficients with signal
    X_mod = conj(W) .* X(:,n2);
    % Calculate sum of all signals (beamformer output)
    Y_sum = sum(X_mod).';
    % Calculating the post filer
    switch filter_nr
        case {'ZEL88','ZEL88p'}
            % Use time-aligned for calculating
            % post-filter
            X = (conj(d0).*X(:,n2)).';
            % Calculate Postfilter-Coefficients of the Zelinski Postfilter
            % with Welch-estimation / of the Zelinski Filter with
            % Welch-estiamtion and Zelinskis post-processing method
            [H_post,alpha,auto_dens,cross_dens] = post_filter(filter_nr,cn,aa,alpha,X,auto_dens,cross_dens);
        case 'SIM92'
            % Use time-aligned for calculating
            % post-filter
            X = (conj(d0).*X(:,n2)).';
            % Calculate Postfilter-Coefficients of the Simmer92 Postfilter
            % with Welch-estimation
            [H_post,alpha,cross_dens,Pyy] = post_filter(filter_nr,cn,aa,alpha,X,Y_sum,cross_dens,Pyy);
        case 'MCC03'
            % Use time-aligned for calculating
            % post-filter
            X = (conj(d0).*X(:,n2)).';
            % Calculate Postfilter-Coefficients of the MCowan03 Postfilter
            [H_post,alpha,auto_dens,cross_dens] = post_filter(filter_nr,cn,aa,alpha,X,Gamma_comb,auto_dens,cross_dens);
        case 'GMCC'
            % Use time-aligned for calculating
            % post-filter
            X = (conj(d0).*X(:,n2)).';
            % Calculate Postfilter-Coefficients of the McCowan03
            % Postfilter, correct solution
            [H_post,alpha,auto_dens,cross_dens] = post_filter(filter_nr,cn,aa,alpha,X,Gamma_comb,d0,auto_dens,cross_dens);
        case 'APAB'
            % Use time-aligned or beamformed signal for calculating
            % post-filter
            X = (conj(d0).*X(:,n2)).';
            % Calculate Postfilter-Coefficients of the APAB Postfilter
            [H_post,alpha,Pxx,Pyy] = post_filter(filter_nr,cn,aa,alpha,X,Y_sum,Pxx,Pyy);
        case 'APES'
             % Calculate output of a Delay&Sum Beamformer
            Y = conj(d0./K).*X(:,n2);
            Y = sum(Y).';
            X = (X(:,n2)).';
            % Calculate Postfilter-Coefficients of the APES Postfilter
            [H_post,alpha,Pxx,Pyy,Pzz] = post_filter(filter_nr,cn,aa,alpha,X,Y,Y_sum,Pxx,Pyy,Pzz);
    end
    % Calculate output of the postfilter
    Y = Y_sum.*H_post;
    % Highpass and Lowpass filter to cut frequencies
    Y = h_high.*Y;
    Y = h_low.*Y;
    % IFFT and Overlap-Add (OLA)
    Y = [Y;conj(Y(end-1:-1:2))];
    yb = (ifft(Y,N));
    y_sv(k1) = y_sv(k1) + yb(:);

    % *********************************************************************
    % Slave-Algorithm for Speech- and Noise-only signals
    % calculating signal at beamformer output, speech-only at beamformer output,
    % noise-only at beamformer output, speech-only at postfilter output and
    % noise-only at postfilter output, depending on input and output arguments
    if nargout>2
        % FFT - Filterbank with Hanning-Windowing for speech- and noise-only
        S = fft(speech(k1,:) .* H,N).';
        V = fft(noise(k1,:) .* H,N).';
        % Multiplicate Beamformer weighting coefficients with signal
        S_mod = conj(W) .* S(:,n2);
        V_mod = conj(W) .* V(:,n2);
        % Calulate Beamformer outputs of speech- and noise-only
        S_sum = sum(S_mod).';
        V_sum = sum(V_mod).';
        switch nargout
            case 4
                % Calc. signal at beamformer output
                % Highpass and Lowpass filtering of beamfomer output
                Y_beam = h_low.*h_high.*Y_sum;
                % IFFT and Overlap-Add (OLA)
                Y_beam = [Y_beam;conj(Y_beam(end-1:-1:2))];
                yb = (ifft(Y_beam,N));
                y_beam(k1) = y_beam(k1) + yb(:);
            case 6
                % Calc. signal at beamformer output
                % Highpass and Lowpass filtering of beamfomer output
                Y_beam = h_low.*h_high.*Y_sum;
                % IFFT and Overlap-Add (OLA)
                Y_beam = [Y_beam;conj(Y_beam(end-1:-1:2))];
                yb = (ifft(Y_beam,N));
                y_beam(k1) = y_beam(k1) + yb(:);
                % Calc. speech-only signal at beamformer output
                % Highpass and Lowpass filtering of beamfomer output
                Y_bs = h_low.*h_high.*S_sum;
                % IFFT and Overlap-Add (OLA)
                Y_bs = [Y_bs;conj(Y_bs(end-1:-1:2))];
                yb = (ifft(Y_bs,N));
                y_bs(k1) = y_bs(k1) + yb(:);
                % Calc. noise-only signal at beamformer output
                % Highpass and Lowpass filtering of beamfomer output
                Y_bv = h_low.*h_high.*V_sum;
                % IFFT and Overlap-Add (OLA)
                Y_bv = [Y_bv;conj(Y_bv(end-1:-1:2))];
                yb = (ifft(Y_bv,N));
                y_bv(k1) = y_bv(k1) + yb(:);
        end
        % Speech-only signal
        % Calc. at output of the postfilter
        Y_s = S_sum.*(H_post);
        % Highpass and Lowpass filtering
        Y_s = h_low.*h_high.*Y_s;
        % IFFT and Overlap-Add (OLA)
        Y_s = [Y_s;conj(Y_s(end-1:-1:2))];
        yb = (ifft(Y_s,N));
        y_s(k1) = y_s(k1) + yb(:);
        % Noise-only signal
        % Calc. at output of the postfilter
        Y_v = V_sum.*(H_post);
        % Highpass and Lowpass filtering
        Y_v = h_low.*h_high.*Y_v;
        % IFFT and Overlap-Add (OLA)
        Y_v = [Y_v;conj(Y_v(end-1:-1:2))];
        yb = (ifft(Y_v,N));
        y_v(k1) = y_v(k1) + yb(:);
    elseif nargout == 2 & solo
    % Calc. signal at beamformer output
    % Highpass and Lowpass filtering of beamfomer output
    Y_beam = h_low.*h_high.*Y_sum;
    % IFFT and Overlap-Add (OLA)
    Y_beam = [Y_beam;conj(Y_beam(end-1:-1:2))];
    yb = (ifft(Y_beam,N));
    y_beam(k1) = y_beam(k1) + yb(:);
    end
end
% Output signal of the whole postfilter algorithm
% Realpart and scaled
y_sv = real(y_sv(n_orig)*1/L*2*fact);
% Calc. realpart and scale Slave-output-signals
if nargout>2
    y_s = real(y_s*1/L*2*fact);
    y_v = real(y_v*1/L*2*fact);
    varargout{1} = y_s(n_orig);
    varargout{2} = y_v(n_orig);
switch nargout
    case 4
        y_beam = real(y_beam*1/L*2*fact);
        varargout{3} = y_beam(n_orig);
    case 6
        y_beam = real(y_beam*1/L*2*fact);
        y_bs = real(y_bs*1/L*2*fact);
        y_bv = real(y_bv*1/L*2*fact);
        varargout{3} = y_beam(n_orig);
        varargout{4} = y_bs(n_orig);
        varargout{5} = y_bv(n_orig);
end
elseif nargout == 2 & solo
    y_beam = real(y_beam*1/L*2*fact);
    varargout{1} = y_beam(n_orig);
end