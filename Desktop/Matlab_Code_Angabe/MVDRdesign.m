function sig = MVDRdesign(cfg, sig)

%need  the signal covariance matrix, scov,
%which represents the power in each signal and the correlation between signals
%also need  the spatial noise covariance matrix, ncov.
%This value represents the noise power on each sensor as well as the
%correlation of the noise between sensors

% cov = sensorcov(cfg.pos,cfg.ang,ncov,scov);
% 
% W = mvdrweights(cfg.pos,cfg.ang,cov);
% 
% Y(idx_nu,:) = W(:,idx_nu)'*squeeze(X(idx_nu,:,:)).';
% Y_des(idx_nu,:) = W(:,idx_nu)'*squeeze(X_des(idx_nu,:,:)).';
% Y_int(idx_nu,:) = W(:,idx_nu)'*squeeze(X_int(idx_nu,:,:)).';
cfg.fstep = 10; % frequency steps between 'sampling points' of set of
% lower and upper limit in order to create extended frequency range
cfg.lfreq_lim = 100; %lower limit
cfg.hfreq_lim = 200; %to create higher limit
cfg.frange = 300:cfg.fstep:(cfg.fs/2 - cfg.hfreq_lim);
cfg.frange = (cfg.frange(1)-cfg.lfreq_lim):cfg.fstep:(cfg.frange(end)+cfg.hfreq_lim);
%p = nextpow2(length(cfg.frange));
%% frequency resolution = 16000/400 = 40Hz fs/N
%%  25ms 400/16000 N/fs
%%
N=400; 
H=100;
win = window(@hamming,N);
%S [freqx  time nmic] matrix
for idx_mic=1:cfg.nmic
    %[S(:,:,idx_mic) F(:,:,idx_mic) T(:,:,idx_mic)] = spectrogram(sig.x(:,idx_mic),win,N-H,length(win),cfg.fs);
    S(:,:,idx_mic) = STFTana(sig.x(:,idx_mic),64,win);
    %[U(:,:,idx_mic) UF(:,:,idx_mic) UT(:,:,idx_mic)] = spectrogram(sig.xSrc(:,idx_mic,2),win,N-H,length(win),cfg.fs);
end
%Fcoef = (1:1:size(S,1))*cfg.fs/N;
%Tcoef = (1:1:size(S,2))*H/cfg.fs;
%S = freq x frame matrix
Phi_xx = Recursive_PSD_estimation(S,0.75);
Phi_uu = Recursive_PSD_estimation(U,0.75);
mic_position.x = [-7.5e-2 -3.375e-2 0 3.375e-2 7.5e-2];
mic_position.y = [-6e-2 -0.75e-2 0 -0.75e-2 -6e-2];
mic_position.z = [-4e-2 0 0 0 -4e-2]; 
steerV = steeringVector(cfg,mic_position,UF);
d=9;
for idx_mic=1:cfg.nmic
    for idx_freq=1:length(F)
        for idx_frame=1:length(T)
            w(idx_freq,idx_frame,idx_mic) = pinv(Phi_uu(idx_freq,idx_frame,idx_mic))*steerV(idx_freq,idx_mic)/...
            (steerV(idx_freq,idx_mic)'*pinv(Phi_uu(idx_freq,idx_frame,idx_mic))*steerV(idx_freq,idx_mic));
        end
    end
end
b=0;
