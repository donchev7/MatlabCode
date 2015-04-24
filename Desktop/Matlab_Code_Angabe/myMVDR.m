function [sig, W] = myMVDR(cfg, sig)
%MYMVDR performs beamforming using the classical MVDR beamformer
% Input parameters:
%   * sig:  struct containing the source and microphone signal(s)
%   * cfg:  struct containing configuration parameters
% 
% Output parameters:
%   * sig:  struct that will also contain the output signal
%   * W:    matrix containing beamforming coefficients for each subband

cfg.mic_pos.x = [-7.5e-2 -3.375e-2 0 3.375e-2 7.5e-2];
cfg.mic_pos.y = [-6e-2 -0.75e-2 0 -0.75e-2 -6e-2];
cfg.mic_pos.z = [-4e-2 0 0 0 -4e-2]; 
cfg.frange = 0:250:8000;
cfg.c = 343; %speed of sound in air
cfg.k_range = 2*pi*cfg.frange/cfg.c;
cfg.K = 64; %number of channels or bands
cfg.N = 16; %decimation ratio
Lp = 128; %length of prototype filter

cfg.p=IterLSDesign(Lp,cfg.K,cfg.N);
%--------------------------------------------------------------------------
%create input data blocks and create subbands
%--------------------------------------------------------------------------
for idx_mic = 1:cfg.nmic    
    %get the frequency subbands for each block -> X of dimension (#subbands, #blocks, #microphones)
    %of microphone data
    X(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.x(:,idx_mic),cfg.K,cfg.N,cfg.p);
    %of desired signal components
    X_des(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.xSrc(:,idx_mic,1),cfg.K,cfg.N,cfg.p);
    %of interference+noise
    if cfg.noise_type
        X_int(:,:,idx_mic) = DFTAnaRealEntireSignal(sum(sig.xSrc(:,idx_mic,2:end),3)...
            +sig.xnoise(:,idx_mic), cfg.K,cfg.N,cfg.p);
    else
        X_int(:,:,idx_mic) = DFTAnaRealEntireSignal(sum(sig.xSrc(:,idx_mic,2:end),3),cfg.K,cfg.N,cfg.p);
    end
end

%--------------------------------------------------------------------------
%estimate spectral cross correlation matrix for each subband
%--------------------------------------------------------------------------
Pxx = zeros(size(X_des,3),size(X_des,3),size(X_des,1));
%easier case ->Y consider only noise components
% for idx_nu = 1:size(X_des,1)
%     Pxx(:,:,idx_nu) = cov(squeeze(X_int(idx_nu,:,:)));
% end
for idx_nu = 1:size(X_des,1)
    Pxx(:,:,idx_nu) = cov(squeeze(X(idx_nu,:,:)));
end
Pxx = Pxx./size(X,2);


%load hrtfs if required
if strcmp(cfg.design,'hrtf')
    load(cfg.path_hrirs);
    HRTF = fft(hrir,cfg.fs,1);
    %determine index of source direction    
%     idx_sourceDir = cfg.idx_sourcePosition(cfg.position(1));
end 

%--------------------------------------------------------------------------
%perform MVDR beamforming for each subband and frame
%--------------------------------------------------------------------------
%matrix for output blocks
Y = zeros(size(X,1), size(X,2));
Y_des = zeros(size(X,1), size(X,2));
Y_int = zeros(size(X,1), size(X,2));
H = zeros(size(X,3), size(X,1));

for idx_nu = 1:size(X,1)
    switch cfg.design
        case 'freefield'
            %create wavevector according to van Trees (2.25)
            kvec = - cfg.k_range(idx_nu) * [sind(cfg.look_elevation)*cosd(cfg.look_azimuth);...
                sind(cfg.look_elevation)*sind(cfg.look_azimuth);...
                cosd(cfg.look_elevation)];
            %create steering vector according to van Trees (2.28) for all
            %microphone positions
            v_k = exp(-1i*kvec.'*[cfg.mic_pos.x; cfg.mic_pos.y; cfg.mic_pos.z]).';
        case 'hrtf'
            if idx_nu == 1
                v_k = HRTF(1,:,cfg.position(1)).';
            else
                v_k = HRTF(ceil(cfg.frange(idx_nu)),:,cfg.position(1)).';
            end
    end
    %compute beamforming weights (according to (6.74) in vanTrees - Optimum
    %Array Processing
    rho = 0;%0.00001; %regularization constant for diagonal loading
    %estimate filter coefficients
    W(:,idx_nu) = ((v_k'*inv(Pxx(:,:,idx_nu)+rho*eye(cfg.nmic))) / ...
        (v_k'*inv(Pxx(:,:,idx_nu)+rho*eye(cfg.nmic))*v_k)).';
    
    %perform beamforming frequency-band-wise
    Y(idx_nu,:) = W(:,idx_nu)'*squeeze(X(idx_nu,:,:)).';
    Y_des(idx_nu,:) = W(:,idx_nu)'*squeeze(X_des(idx_nu,:,:)).';
    Y_int(idx_nu,:) = W(:,idx_nu)'*squeeze(X_int(idx_nu,:,:)).';
end %for idx_nu

%create time-domain output signal
y = DFTSynRealEntireSignal(Y, cfg.K, cfg.N, cfg.p);
y_des = DFTSynRealEntireSignal(Y_des, cfg.K, cfg.N, cfg.p);
y_int = DFTSynRealEntireSignal(Y_int, cfg.K, cfg.N, cfg.p);

%--------------------------------------------------------------------------
%Set output signal y
%--------------------------------------------------------------------------
sig.y = y;
sig.y_des = y_des;
sig.y_int = y_int;
%--------------------------------------------------------------------------
end