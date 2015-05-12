%main function

clc;
clear;
close all;

%------------------------------------------------------------------
% Initialization
%------------------------------------------------------------------
cfg = [];
sig = [];
flt = [];


%add path of filterbank m-files (required if beamformer design is
% applied)
addpath(genpath('filterbank'));

%add path of beamformer m-files (required if FSB beamformer design is
% applied)
addpath(genpath('Generalized_RLSFI_BF'));
% 
cfg.design = 'freefield';
cfg.beamFormer = 'MVDR';
cfg.geometry = 1; %linear array
cfg.spacing = 0; % 0 = non-uniform array, 1=uniform array;
% cfg.design = 'hrtf';
cfg.wng_limit_db = -15;
cfg.look_azimuth = 90; %azimuth angle
cfg.look_elevation = repmat(90 - atand(0.73/1.1), length(cfg.look_azimuth)); %elevation
cfg.des_look_dir.azimuth = cfg.look_azimuth;
cfg.des_look_dir.elevation = cfg.look_elevation;
cfg.nsrc = 2;       % Number of sources
cfg.nmic = 5;       % Number of mics
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 9 29]; %position of sources, interferers at 40° and 140°
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 10 32]; %position of sources, interferers at 45° and 155°
cfg.position = [19 46 22]; %position of sources, interferers at 45° and 105°
%source must be first defined and than interferer e.x. cfg.position(1) must
%be the desired source
cfg.noise_type = 0;
% noise implemented in SetAcousticScenario.m
    % 0 = no noise
    % 1 = white noise
    % 2 = noise from files
    % 3 = generated diffuse noise
cfg.sig_len = 0;    % Choose 0 if the whole signal should be used (if ASR scores need to be evaluated!)
%------------------------------------------------------------------
% Set Acoustical scenario parameters (parameters are stored in cfg
% structure
%------------------------------------------------------------------
cfg = SetAcousticScenario(cfg);
cfg = BF_Array_Geometry(cfg);
%------------------------------------------------------------------
% Create microphone signals (signals are stored in sig structure, RIRs in
% flt structure
%------------------------------------------------------------------
[cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);

%------------------------------------------------------------------
%% filterbank initialization
cfg.c = 342; %speed of sound in air
cfg.fs = 16000;
% cfg.N = 512; % FFT points
% cfg.K = 128; % frame shift
% cfg.wlen = 256; %window length
cfg.N = 512; % FFT size
cfg.K = 128; % frame shift
cfg.Lp = 1024; % prototype filter length
%p2=IterLSDesign(cfg.Lp,cfg.N,cfg.K);
load('/filterbank/prototype_K512_N128_Lp1024.mat');
cfg.p = p; clear p;
cfg.frange = linspace(0,cfg.fs/2,cfg.N/2+1)'; % frequency axis
cfg.k_range = 2*pi*cfg.frange/cfg.c;
cfg.angRange.azimuth = [(-180):5:(180)];
cfg.angRange.elevation = repmat(cfg.look_elevation,size(cfg.angRange.azimuth));
%------------------------------------------------------------------

%Filterbank Analysis
for idx_mic=1:cfg.nmic
 X(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.x(:,idx_mic),cfg.N,cfg.K,cfg.p);
 if cfg.nsrc >1
    X_int(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.xSrc(:,idx_mic,2:end),cfg.N,cfg.K,cfg.p);
 end
end
%------------------------------------------------------------------------------

%------------------------------------------------------------------
%perform beamformer design (of robust FSB)
%------------------------------------------------------------------
[flt.w.RFSB, cfg, steerV, realWNG_dB] = RobustFSBdes(cfg);

%flt.w.RFSBFIR = main_robust_FSBORG(cfg, cfg.look_azimuth, cfg.look_elevation, cfg.design, cfg.wng_limit_db);


% ------------------------------------------------------------------
% Create microphone signals (signals are stored in sig structure) at
% beamformer output signal
% ------------------------------------------------------------------
%sig.y.RFSBFIR = zeros(length(sig.x)+length(flt.w.RFSBFIR)-1,1);
%ySrc = number of source signals filtered with w - signal at beamer output
% sigu = both sources mixed and filtered with w - signal at beamer output
%sig.ySrc = zeros(length(sig.x)+length(flt.w.RFSBFIR)-1,cfg.nsrc); 
%  for idx_channels = 1:cfg.nmic
%      sig.y.RFSBFIR  = sig.y.RFSBFIR + fftfilt(flt.w.RFSBFIR(:,idx_channels),[sig.x(:,idx_channels); zeros(length(flt.w.RFSBFIR)-1,1)]);
%      
% %      for idx_sources = 1:cfg.nsrc
% %          sig.ySrc(:,idx_sources) = sig.ySrc(:,idx_sources) + fftfilt(flt.w.RFSBFIR(:,idx_channels), ...
% %              [sig.xSrc(:,idx_channels,idx_sources); zeros(length(flt.w.RFSBFIR)-1,1)]);
% %      end
%  end


%------------------------------------------------------------------
% Create microphone signals (signals are stored in sig structure) at
% beamformer output signal
%------------------------------------------------------------------   

% for idx_freq=1:length(cfg.frange) 
% Y(idx_freq,:) = flt.w.RFSB(:,idx_freq)'*squeeze(X(idx_freq,:,:)).';
% end
% sig.y.RFSB = istft(Y, cfg.K, cfg.N, cfg.fs);


% %------------------------------------------------------------------
% % Potential place for postfilter
% %------------------------------------------------------------------

%Calculate the sound field as in [1] Kuklasinski

micPosition = [cfg.mic_pos.x; cfg.mic_pos.y; cfg.mic_pos.z].';
cfg.micSpacing = micSpacingCalculation(micPosition);

%%Evaluating gammaIso the cylindrical sound field from HRIR measurments
% Also evaluating the steering vector d from the RIR measurments as
% mentioned in Kuklasinski paper
HRIR=load('./Data/RIRs/HRIR_NAO_LRC_az_0_355_16kHz_norm.mat');
HRTF = fft(HRIR.imp_resp,512,1);
%L=length(HRIR.imp_resp(:,1,1));
[Nx_orig,K_signal,Ang] = size(HRIR.imp_resp);
[Nx,Ang,K] = size(flt.h);
padding_zeros = ceil(1058/cfg.N)*cfg.N - 1058;
n_zeros = ceil(Nx_orig/cfg.N)*cfg.N - Nx_orig;
for idx_mic=1:cfg.nmic
    %RIR_Padded = [flt.h(1:1058,cfg.look_azimuth/5+1,idx_mic); zeros(padding_zeros,1)];
    %d(idx_mic,:) = fft(RIR_Padded,cfg.N);
    for idx_ang=1:size(HRIR.imp_resp,3)
        HRIR_Padded = [HRIR.imp_resp(:,cfg.mic_ch(idx_mic),idx_ang); zeros(n_zeros,1)];
        d_hrtf_tmp(:,idx_mic,idx_ang) = fft(HRIR_Padded,cfg.N);
    end
end
% d=d(:,1:cfg.N/2+1)./Nx;
d_hrtf_tmp = d_hrtf_tmp(1:cfg.N/2+1,:,:);
for idx_ang=1:size(HRIR.imp_resp,3)
    d_hrtf(:,:,idx_ang) = cov(d_hrtf_tmp(:,:,idx_ang));
end
gammaIsoSum = zeros(5,5);
for idx_ang=1:size(HRIR.imp_resp,3)
     gammaIsoSum = gammaIsoSum+d_hrtf(:,:,idx_ang)*d_hrtf(:,:,idx_ang)';
end
gammaIso=gammaIsoSum./size(HRIR.imp_resp,3);


%% End of evalution of gammaIso and steering vector d

%%Estimating the PSDs with recursive averaging
% counter=0;
% for m = 1:(cfg.nmic-1)
%     for n = (m+1):cfg.nmic
%         counter = counter + 1;
%         %Thi_Ycpsd(:,:,counter) = Recursive_PSD_estimation(X_int(:,:,m),X_int(:,:,n),0.76);
%         Thi_Ycpsd(:,:,counter) = estimate_cpsd(X(:,:,m),X(:,:,n),0.76);
%     end
% end
% for m=1:cfg.nmic
%     Thi_Ypsd(:,:,m) = estimate_psd(X(:,:,m),0.76);
%     %Thi_Ypsd(:,:,m) = Recursive_PSD_estimation(X_int(:,:,m),0.76);
% end
%%End of PSD estimation


mue = 10^(-cfg.wng_limit_db/10); %mue is used regularization
%EQ 2.21 Optimum Array Processing et. L. Van Trees
tau = micPosition*[sind(cfg.look_elevation)*cosd(cfg.look_azimuth); sind(cfg.look_elevation)*sind(cfg.look_azimuth); cosd(cfg.look_elevation)];
for idx_freq=1:size(X,1)
    beta = (2*cfg.frange(idx_freq))/cfg.c;
    Gamma_tmp_cylndrical = besselj(0,beta*cfg.micSpacing); %2D cylindrical noise model
    Gamma_tmp_diffuse = sinc(beta*cfg.micSpacing);
    if strcmp(cfg.beamFormer,'MVDR')
        Gamma_const_cylindrical = tril(Gamma_tmp_cylndrical,-1)./(1+mue) + diag(diag(Gamma_tmp_cylndrical)) + ...
        triu(Gamma_tmp_cylndrical,1)./(1+mue);
        Gamma_const_diffuse = tril(Gamma_tmp_diffuse,-1)./(1+mue) + diag(diag(Gamma_tmp_diffuse)) + ...
        triu(Gamma_tmp_diffuse,1)./(1+mue);
    else
        Gamma_const_cylindrical = Gamma_tmp_cylndrical;
    end
    GammaCyldircal(idx_freq,:,:) = Gamma_const_cylindrical;
    GammaDiffuse(idx_freq,:,:) = Gamma_const_diffuse;
    
    d=d_hrtf_tmp(idx_freq,:,cfg.position(1)).';
    %EQ 2.25 Optimum Array Processing et. L. Van Trees
    k_vec = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tau;
    %v_k EQ 2.28 Optimum Array Processing et. L. Van Trees
    d0(:,idx_freq) = exp(-1i*k_vec);
    %X_aligned(idx_freq,:,:) = bsxfun(@times,d0(:,idx_freq)',squeeze(X(idx_freq,:,:)));
    %X_int_aligned(idx_freq,:,:)=bsxfun(@times,d0(:,idx_freq)',squeeze(X_int(idx_freq,:,:)));
    Pxx2(:,:,idx_freq) = cov(squeeze(X(idx_freq,:,:)));
    %Pxx(:,:,idx_freq) = Thi_Y(squeeze(Thi_Ycpsd(idx_freq,:,:)),squeeze(Thi_Ypsd(idx_freq,:,:)),cfg.nmic)./size(X_int,2);
    %Thi_Y(:,:,idx_freq) = cov(squeeze(X(idx_freq,:,:)))./size(X,2);
    
    %Wmvdr_superDirective = filter for freefield 2D cyldirical noise and freefield steering vector d0
    Wmvdr_superDirective(:,idx_freq) = myMVDR2(squeeze(GammaDiffuse(idx_freq,:,:)),d0(:,idx_freq),cfg.beamFormer);
    Wmvdr_freeField(:,idx_freq) = myMVDR2(squeeze(GammaCyldircal(idx_freq,:,:)),d0(:,idx_freq),cfg.beamFormer);
    
    %WgammaIso = filter with evaluated cylindrical sound field from HRTF and
    %hrtf evaluated steering vector d
    WgammaIso(:,idx_freq) = myMVDR2(gammaIso,d,cfg.beamFormer);
    
    %Wmwf = the multichannel Wiener filter EQ 5b in Kuklasinski paper
    %but for freefield parameters
    Wmwf_freeField(:,idx_freq) = postFilter(Wmvdr_freeField(:,idx_freq),squeeze(GammaCyldircal(idx_freq,:,:)),d0(:,idx_freq),Pxx2(:,:,idx_freq));
    
    %Wmwf as defined in Kuklasinski et. al EQ 5b 
    Wmwf(:,idx_freq) = postFilter(WgammaIso(:,idx_freq),gammaIso,d,Pxx2(:,:,idx_freq));
    %Conventional MVDR beamformer with covariance matrix of microphone signals as input
    if cfg.nsrc > 1
        Wmvdr(:,idx_freq) = myMVDR2(cov(squeeze(X_int(idx_freq,:,:))),d0(:,idx_freq),cfg.beamFormer);
    else
        Wmvdr(:,idx_freq) = myMVDR2(cov(squeeze(X(idx_freq,:,:))),d0(:,idx_freq),cfg.beamFormer);
    end
    Ymwf_freeField(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf_freeField(:,idx_freq);
    Ymwf(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf(:,idx_freq);
    Ymvdr(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr(:,idx_freq);
    Ymvdr_superDirective(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_superDirective(:,idx_freq);
    Yrfsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*flt.w.RFSB(:,idx_freq);
end

ymvdr = DFTSynRealEntireSignal(Ymvdr, cfg.N, cfg.K, cfg.p);
ymvdr_superDirective = DFTSynRealEntireSignal(Ymvdr_superDirective, cfg.N, cfg.K, cfg.p);
ymwf = DFTSynRealEntireSignal(Ymwf, cfg.N, cfg.K, cfg.p);
ymwf_freeField = DFTSynRealEntireSignal(Ymwf_freeField, cfg.N, cfg.K, cfg.p);
yrfsb = DFTSynRealEntireSignal(Yrfsb, cfg.N, cfg.K, cfg.p);
% %------------------------------------------------------------------
% % BEAMPATTERN for MVDR and RFSB
% beampattern(steerV,W,flt.w.RFSB,cfg.angRange.azimuth,cfg.frange)
% %
% %------------------------------------------------------------------


%y_sv = sim_system(sig.xSrc(:,:,1),sig.xSrc(:,:,2),cfg.look_azimuth,'SDB',-10,'GMCC',0,0,0.8,micPosition,512,4,16000,200,7600);

% %------------------------------------------------------------------
% % Evaluating the fwsegsnr, PESQ and ASR scores
% %------------------------------------------------------------------

[y_rfsb_pesq, y_rfsb_fwsegsnr, y_rfsb_ASR]=evaluateScores(sig,yrfsb,cfg);
[y_mvdr_pesq, y_mvdr_fwsegsnr, y_mvdr_ASR]=evaluateScores(sig,ymvdr,cfg);
[y_mvdrSuperDirective_pesq, y_mvdrSuperDirectiv_fwsegsnr, y_mvdrSuperDirectiv_ASR]=evaluateScores(sig,ymvdr_superDirective,cfg);
[y_mwfFreeField_pesq, y_mwfFreeField_fwsegsnr, y_mwfFreeField_ASR]=evaluateScores(sig,ymwf_freeField,cfg);
[y_mwf_pesq, y_mwf_fwsegsnr, y_mwf_ASR]=evaluateScores(sig,ymwf,cfg);



fprintf('\n');
fprintf('\n');
fprintf('                    |----PESQ scores----|-----FwSegSNR----|----ASR----------------\n');
fprintf('           RFSB     |      %g        |      %.2f      |    %g             \n',y_rfsb_pesq,y_rfsb_fwsegsnr,y_rfsb_ASR);
fprintf('           MVDR     |      %g        |      %.2f      |    %g             \n',y_mvdr_pesq,y_mvdr_fwsegsnr,y_mvdr_ASR);
fprintf('MVDR_SuperDirective |      %g        |      %.2f      |    %g             \n',y_mvdrSuperDirective_pesq,y_mvdrSuperDirectiv_fwsegsnr,y_mvdrSuperDirectiv_ASR);
fprintf('MVDR+MWF_FreeField  |      %g        |      %.2f      |    %g             \n',y_mwfFreeField_pesq,y_mwfFreeField_fwsegsnr,y_mwfFreeField_ASR);
fprintf('   MVDR +  MWF      |      %g        |      %.2f      |    %g             \n',y_mwf_pesq,y_mwf_fwsegsnr,y_mwf_ASR);
fprintf('----------------------------------------------------------------------------------')
