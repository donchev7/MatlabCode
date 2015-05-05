function beampattern(beam_nr,phi_d,mue,mics,fs,varargin)
% beampattern(beam_nr,phi_d,mue,mics,fs,sim,phi_n)
% Plots the Beampattern of a 1-dimensional Array and the frequency response
% at a defined angle, as well as the frequncy response for small
% perturbations imposed to mic positions
% beam_nr
% phi_d
% mue
% mics
% fs
% sim
% phi_n
% you can choose between the following beamformers
% 'DSB' ... Delay&Sum-Beamformer
% 'SDB' ... Superdirective Beamformer
% default beam_nr = 'SDB';
% angle to TDOA
% for the regularization of Coherence Matrix Gamma
% mue in dB; default mue = -20
% microphone positons; default load mics_8xh.mat
% sampling frequency in Hz; default fs = 16000 (16kHz)
% use simulated coherence function or sinc-function for
% Gamma; default sim = 'sinc'
% 'incoh' ............. incoherent Noisefield
% 'sinc' .............. Sinc-Function
% 'bessel' ............ Bessel-Function
% 'zero' .............. puts a Zero at the angle
% specified by phi_n
% angle for the plotted frequency response; if sim = 'zero'
% phi_n is the angle of the specified zero
% used functions: mvdr.m
if nargin>=7 phi_n = varargin{2}; end
if nargin>=6 sim = varargin{1}; end
if nargin<6 sim = 'sinc'; end
if nargin<5 fs = 16000; end
if nargin<4 load mics_8xh.mat; end
if nargin<3
    help beampattern
    return;
end
% Parameters
theta_d = 90;
% elvation angle to direction of arrival (fixed)
% Check microphone position matrix
[K,Dim] = size(mics);
if (Dim < 2) | (K < 1)
    error('bad matrix of microphine coordinates');
end
% Calculate mue
if mue ~= 0
    mue = 10^(mue/10);
end
%Calc. of angles in rad
theta_r = theta_d(:).' * pi / 180;
phi_r = phi_d(:).' * pi / 180;
phi_n = phi_n(:).' * pi / 180;
% Calculate time alignment vector
ed = [sin(theta_r).*cos(phi_r); sin(theta_r).*sin(phi_r)];
Rc = mics*ed;
% Define Matrix of the micorphone distances
xc = mics(:,1);
xc = xc(:,ones(K,1));
dR = (xc-xc.');
% Define Frequency vector
if fs <= 12000
    f = linspace(0,3400,120);
else
    f = linspace(0,6800,240);
end
Nf = length(f);
% Calculate Coherencefunction
for l = 1:Nf
    switch sim
        case 'incoh'
            % Calc. coherencefunction for an incoherent noisefield
            Gamma_const = diag(ones(1,K));
        case 'bessel'
            % Calc. coherencefunction for a zylindrical isotropic noisefield
            beta = 2*pi*f(l)/340;
            Gamma_dum = besselj(0,beta*dR);
        case 'sinc'
            % Calc. coherencefunction for a diffuse noisefield
            beta = 2*f(l)/340;
            Gamma_dum = sinc(beta*dR);
        case 'zero'
            % Calc. coherencefunction for a coherent noisefield (interferer noise from angle phi_n)
            beta = 2*pi*f(l)/340;
            Gamma_real = cos(beta*cos(phi_n)*dR);
            Gamma_imag = -sin(beta*cos(phi_n)*dR);
            Gamma_dum = (Gamma_real + j*Gamma_imag);
    end
    if strcmp(sim,'zero') | strcmp(sim,'bessel') | strcmp(sim,'sinc')
        % regularization of the coherence matrix
        Gamma_const = tril(Gamma_dum,-1)./(1 + mue) + diag(diag(Gamma_dum)) + ...
        triu(Gamma_dum,1)./(1 + mue);
    end
    Gamma(:,:,l) = Gamma_const;
end
% Calculate Beampattern
phi_wav_d = [(-180):1:(180)];
% angle vector
phi_wav = phi_wav_d(:).' * pi / 180;
ed_wav = [sin(theta_r).*cos(phi_wav); sin(theta_r).*sin(phi_wav)];
Rc_wav = mics*ed_wav;
for l = 1:Nf
    % Calc. beamformer coefficients and steering vector for each angle and frequency
    [W(:,l),dum] = mvdr(Rc,Gamma(:,:,l),f(l),beam_nr);
    beta_wav = 2*pi*f(l)/340;
    d = exp(-j * beta_wav *Rc_wav);
    % Calc. gain for each angle and frequency
    H = abs(W(:,l)'*d).^2;
    HdB = max(-25,10*log10(H + eps));
    H_log(:,l) = HdB;
end
% Plot beampattern
figure,surf(f,phi_wav_d,H_log);
axis tight
set(gca,'TickDir','out');
set(gca,'YTick',[-180:45:180]);
if fs <= 12000
    set(gca,'XTick',[100;1000;2000;3000;3400])
    set(gca,'XTickLabel',{'100';'1000';'2000';'3000';'3400'})
else
    set(gca,'XTick',[100;1000;2000;3000;4000;5000;6000;6800])
    set(gca,'XTickLabel',{'100';'1000';'2000';'3000';'4000';'5000';'6000';'6800'})
end
colorbar('YLim',[-25 max(max(H_log))])
view([0,90]);
box on
shading interp
ylabel('\theta in Grad');
xlabel('Frequenz [Hz]');
% Calculate frequency response
if strcmp(sim,'zero') | ~isstr(sim)
    phi_freq = phi_n;
else
    phi_freq = phi_r;
end
figurebackcolor = 'black';
Hf = calc_freq_resp(mics,W,theta_r,phi_n,f,fs);
HfdB = max(-100,10*log10(Hf + eps));
% repeat for small perturbations imposed to mic positions
delta = 0.001;
% standard deviation in m
r = mics + delta*randn(size(mics));
Hr = calc_freq_resp(r,W,theta_r,phi_n,f);
HrdB = max(-100,10*log10(Hr + eps));
% Plot frequency response
pos = [0.045 0.01 0.4 0.37];
fp4 = figure('numbertitle','off','name','Frequency domain',...
'Units','normal','Position',pos);
colordef(fp4,figurebackcolor);
plot(f,HfdB,f,HrdB);
grid on;
xlabel('f in Hz');
ylabel('magnitude in dB');
title('frequency responses of ideal (y), and perturbated array (m)');
% -------------------------------------------------------------------------
% Calc. frequency response
end
function Hf = calc_freq_resp(mics,W,theta_r,phi_r,f)
    beta = (2*pi*f/340);
    ed = [sin(theta_r).*cos(phi_r); sin(theta_r).*sin(phi_r)];
    % Calc. of Constraint Matrix
    d = exp(-1i*(mics*ed)*beta);
    Hf = abs(diag(W'*d)).^2;
end