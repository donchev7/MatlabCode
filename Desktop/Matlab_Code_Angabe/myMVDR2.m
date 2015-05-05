function [sig, W] = myMVDR2(cfg, sig)
%MYMVDR performs beamforming using the classical MVDR beamformer
% Input parameters:
%   * sig:  struct containing the source and microphone signal(s)
%   * cfg:  struct containing configuration parameters
% 
% Output parameters:
%   * sig:  struct that will also contain the output signal
%   * W:    matrix containing beamforming coefficients for each subband

%--------------------------------------------------------------------------
%create input data blocks and create subbands
%--------------------------------------------------------------------------

for idx_mic = 1:cfg.nmic    
    %get the frequency subbands for each block -> X of dimension (#subbands, #blocks, #microphones)
    %of microphone data
    %X(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.x(:,idx_mic),cfg.K,cfg.N,cfg.p);
    X(:,:,idx_mic) = stft(sig.x(:,idx_mic),cfg.K,cfg.N,cfg.K,cfg.fs);
    %of desired signal components
    %X_des(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.xSrc(:,idx_mic,1),cfg.K,cfg.N,cfg.p);
    X_des(:,:,idx_mic) = stft(sig.xSrc(:,idx_mic,1),cfg.K,cfg.N,cfg.K,cfg.fs);
    lds_prev = zeros(length(cfg.frange),size(X,2));
    Pxx(:,:,idx_mic) = welch_est(lds_prev,X(:,:,idx_mic),X(:,:,idx_mic),0.76);
    %of interference+noise
    if cfg.noise_type
        %X_int(:,:,idx_mic) = DFTAnaRealEntireSignal(sum(sig.xSrc(:,idx_mic,2:end),3)...
         %   +sig.xnoise(:,idx_mic), cfg.K,cfg.N,cfg.p);
         X_int(:,:,idx_mic) = stft(sum(sig.xSrc(:,idx_mic,2:end),3)...
            +sig.xnoise(:,idx_mic), cfg.K,cfg.N,cfg.K,cfg.fs);
         
    else
        %X_int(:,:,idx_mic) = DFTAnaRealEntireSignal(sum(sig.xSrc(:,idx_mic,2:end),3),cfg.K,cfg.N,cfg.p);
        X_int(:,:,idx_mic) = stft(sum(sig.xSrc(:,idx_mic,2:end),3),cfg.K,cfg.N,cfg.K,cfg.fs);
    end
end
%--------------------------------------------------------------------------
%estimate spectral cross correlation matrix for each subband
%--------------------------------------------------------------------------
%Pxx = zeros(size(X_des,3),size(X_des,3),size(X_des,1));
%easier case ->Y consider only noise components
for idx_nu = 1:size(cfg.iso2,1)
    Cnn2(:,:,idx_nu) = cov(squeeze(cfg.iso2(idx_nu,:,:)));
    %Thi_Y(:,:,idx_nu) = cov(squeeze(Pxx(idx_nu,:,:)));
    Thi_Y(:,:,idx_nu) = cov(squeeze(X(idx_nu,:,:)));
end
Cnn2 = Cnn2./size(cfg.iso2,2);
Cnn = cfg.iso2;
% for idx_nu = 1:size(X_des,1)
%     Pxx(:,:,idx_nu) = cov(squeeze(X(idx_nu,:,:)));
% end
% Pxx = Pxx./size(X,2);
% for idx_mic=1:cfg.nmic
%     iso2(:,:,idx_mic) = besselj(0,(2*cfg.frange*cfg.micSpacing(idx_mic,:))/cfg.c);
% end
% iso2=iso2./size(iso2,2);

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
            mic_position = [cfg.mic_pos.x; cfg.mic_pos.y; cfg.mic_pos.z];
            kvec = - cfg.k_range(idx_nu) * [sind(cfg.look_elevation)*cosd(cfg.look_azimuth);...
                sind(cfg.look_elevation)*sind(cfg.look_azimuth);...
                cosd(cfg.look_elevation)];
            %create steering vector according to van Trees (2.28) for all
            %microphone positions
            v_k = exp(-1i*kvec.'*mic_position).';
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
    mue = 10^(-20/10);
    Gamma_const = tril(squeeze(Cnn(idx_nu,:,:)),-1)./(1+mue) + diag(diag(squeeze(Cnn(idx_nu,:,:)))) + ...
        triu(squeeze(Cnn(idx_nu,:,:)),1)./(1+mue);
    W(:,idx_nu) = (inv(Gamma_const)*v_k) / ...
        (v_k'*inv(Gamma_const)*v_k);
    
    V(:,idx_nu) = 1/(cfg.nmic-1) * trace((eye(cfg.nmic,cfg.nmic)-v_k*W(:,idx_nu)')*Thi_Y(:,:,idx_nu)*inv(Gamma_const));
    
    S(:,idx_nu) = W(:,idx_nu)'*(Thi_Y(:,:,idx_nu)-V(:,idx_nu)*Gamma_const)*W(:,idx_nu);
    
    Wmwf(:,idx_nu) = (S(:,idx_nu)/(S(:,idx_nu)+(V(:,idx_nu)*inv(v_k'*inv(Gamma_const)*v_k))))*...
        W(:,idx_nu);
    
    %perform beamforming frequency-band-wise
    Y(idx_nu,:) = Wmwf(:,idx_nu)'*squeeze(X(idx_nu,:,:)).';
    Y_des(idx_nu,:) = W(:,idx_nu)'*squeeze(X_des(idx_nu,:,:)).';
    Y_int(idx_nu,:) = W(:,idx_nu)'*squeeze(X_int(idx_nu,:,:)).';
end %for idx_nu

%create time-domain output signal
y = istft(Y, cfg.N, cfg.K, cfg.fs);
y_des = istft(Y_des, 128, cfg.K, cfg.fs);
y_int = istft(Y_int, cfg.N, cfg.K, cfg.fs);

%--------------------------------------------------------------------------
%Set output signal y
%--------------------------------------------------------------------------
sig.y.mwf = y;
sig.y_des = y_des;
sig.y_int = y_int;

%--------------------------------------------------------------------------
end