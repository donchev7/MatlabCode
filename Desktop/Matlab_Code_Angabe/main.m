%main function

%clc;
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
cfg.nsrc = 1;       % Number of sources
cfg.nmic = 9;       % Number of mics
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 9 29]; %position of sources, interferers at 40° and 140°
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 10 32]; %position of sources, interferers at 45° and 155°
cfg.position = [37 25 22]; %position of sources, interferers at 30
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
% sig.y.RFSBFIR = zeros(length(sig.x)+length(flt.w.RFSBFIR)-1,1);
% %ySrc = number of source signals filtered with w - signal at beamer output
% %sigu = both sources mixed and filtered with w - signal at beamer output
% sig.ySrc = zeros(length(sig.x)+length(flt.w.RFSBFIR)-1,cfg.nsrc); 
%  for idx_channels = 1:cfg.nmic
%      sig.y.RFSBFIR  = sig.y.RFSBFIR + fftfilt(flt.w.RFSBFIR(:,idx_channels),[sig.x(:,idx_channels); zeros(length(flt.w.RFSBFIR)-1,1)]);
%      
% %      for idx_sources = 1:cfg.nsrc
% %          sig.ySrc(:,idx_sources) = sig.ySrc(:,idx_sources) + fftfilt(flt.w.RFSBFIR(:,idx_channels), ...
% %              [sig.xSrc(:,idx_channels,idx_sources); zeros(length(flt.w.RFSBFIR)-1,1)]);
% %      end
%  end

% %------------------------------------------------------------------
% % Potential place for postfilter
% %------------------------------------------------------------------

d_hrtf = GetConstants(cfg);


micPosition = [cfg.mic_pos.x; cfg.mic_pos.y; cfg.mic_pos.z].';
cfg.micSpacing = micSpacingCalculation(micPosition);

mue = 10^(-cfg.wng_limit_db/10); %mue is used for regularization
%EQ 2.21 Optimum Array Processing et. L. Van Trees
tau = micPosition*[sind(cfg.look_elevation)*cosd(cfg.look_azimuth); sind(cfg.look_elevation)*sind(cfg.look_azimuth); cosd(cfg.look_elevation)];
%tauInt = micPosition*[sind(cfg.look_elevation)*cosd(cfg.position(2)*5-5); sind(cfg.look_elevation)*sind(cfg.position(2)*5-5); cosd(cfg.look_elevation)];
for idx_mic=1:cfg.nmic
 X(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.x(:,idx_mic),cfg.N,cfg.K,cfg.p);
 if cfg.nsrc >1
    X_int(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.xSrc(:,idx_mic,2:end),cfg.N,cfg.K,cfg.p);
 end
end

%Estimating the PSDs with recursive averaging
% counter=0;
% for m = 1:(cfg.nmic-1)
%     for n = (m+1):cfg.nmic
%         counter = counter + 1;
%         %Thi_Ycpsd(:,:,counter) = Recursive_PSD_estimation(X_int(:,:,m),X_int(:,:,n),0.76);
%         Thi_Ycpsd(:,:,counter) = estimate_cpsd(X(:,:,m),X(:,:,n),0.76);
%     end
% end
% for m=1:cfg.nmic
%     Thi_Ypsd(:,:,m) = estimate_psd(X(:,:,m),0.68);
%     %Thi_Ypsd(:,:,m) = Recursive_PSD_estimation(X(:,:,m),0.68);
% end
%End of PSD estimation

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
        Gamma_const_diffuse = Gamma_tmp_diffuse;
    end
    GammaCyldircal = Gamma_const_cylindrical;
    GammaDiffuse = Gamma_const_diffuse;
    gammaIso = covarianceEstimate(squeeze(d_hrtf(idx_freq,:,:)));
    d=squeeze(d_hrtf(idx_freq,cfg.position(1),:));
    %EQ 2.25 Optimum Array Processing et. L. Van Trees
    k_vec = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tau;
    %k_vecInt = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tauInt;
    %v_k EQ 2.28 Optimum Array Processing et. L. Van Trees
    d0 = exp(-1i*k_vec);
%    dint = exp(-1i*k_vecInt);
    %A = [d0, dint];
    %e1=[1,0];
    %wNull = (e1*pinv(A)).';
    %X_aligned(idx_freq,:,:) = bsxfun(@times,d0(:,idx_freq)',squeeze(X(idx_freq,:,:)));
    %X_int_aligned(idx_freq,:,:)=bsxfun(@times,d0(:,idx_freq)',squeeze(X_int(idx_freq,:,:)));
    Pxx = covarianceEstimate(squeeze(X(idx_freq,:,:)));
    
    Wmvdr_freeField = myMVDR2(GammaCyldircal,d0,cfg.beamFormer,cfg.nmic);
    Wmvdr_superDirective = myMVDR2(GammaDiffuse,d0,cfg.beamFormer,cfg.nmic);
    %Wmvdr_superDirective = filter for freefield 2D cyldirical noise and freefield steering vector d0
    WgammaIso = myMVDR2(gammaIso,d,cfg.beamFormer,cfg.nmic);
    %Wmwf = the multichannel Wiener filter EQ 5b in Kuklasinski paper
    %but for freefield parameters
    Wmwf_freeField = postFilter(Wmvdr_freeField,GammaCyldircal,d0,Pxx,cfg.nmic);

    %Wmwf as defined in Kuklasinski et. al EQ 5b 
    Wmwf = postFilter(WgammaIso,gammaIso,d,Pxx,cfg.nmic);
    %Conventional MVDR beamformer with covariance matrix of microphone signals as input
    if cfg.nsrc > 1
        Wmvdr = myMVDR2(covarianceEstimate(squeeze(X_int(idx_freq,:,:))),d0,cfg.beamFormer,cfg.nmic);
    else

        Wmvdr = Pxx*d0/(d0'*Pxx*d0);%myMVDR2(Pxx,d0,cfg.beamFormer,cfg.nmic);
    end
    Ymwf_freeField(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf_freeField;
    Ymwf(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf;
    Ymvdr(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr;
    Ymvdr_superDirective(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_superDirective;
    Yrfsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*flt.w.RFSB(:,idx_freq);
    Ydsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d0/cfg.nmic);
    %Ynull(idx_freq,:) = squeeze(X(idx_freq,:,:))*wNull;
end
Yinput = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,1)));
ymvdr = DFTSynRealEntireSignal(Ymvdr, cfg.N, cfg.K, cfg.p);
ymvdr_superDirective = DFTSynRealEntireSignal(Ymvdr_superDirective, cfg.N, cfg.K, cfg.p);
ymwf = DFTSynRealEntireSignal(Ymwf, cfg.N, cfg.K, cfg.p);
ymwf_freeField = DFTSynRealEntireSignal(Ymwf_freeField, cfg.N, cfg.K, cfg.p);
yrfsb = DFTSynRealEntireSignal(Yrfsb, cfg.N, cfg.K, cfg.p);
yDSB = DFTSynRealEntireSignal(Ydsb, cfg.N, cfg.K, cfg.p);
yinput = DFTSynRealEntireSignal(Yinput, cfg.N, cfg.K, cfg.p);
%ynull = DFTSynRealEntireSignal(Ynull, cfg.N, cfg.K, cfg.p);


% %------------------------------------------------------------------
% % BEAMPATTERN for MVDR and RFSB
% beampattern(steerV,W,flt.w.RFSB,cfg.angRange.azimuth,cfg.frange)
% %
% %------------------------------------------------------------------


%y_sv = sim_system(sig.xSrc(:,:,1),sig.xSrc(:,:,2),cfg.look_azimuth,'SDB',-10,'GMCC',0,0,0.8,micPosition,512,4,16000,200,7600);

% %------------------------------------------------------------------
% % Evaluating the fwsegsnr, PESQ and ASR scores
% %------------------------------------------------------------------

[y_rfsb_pesq, y_rfsb_fwsegsnr, y_rfsb_ASR]=evaluateScores(sig.s(:,1),yrfsb,cfg);
[y_mvdr_pesq, y_mvdr_fwsegsnr, y_mvdr_ASR]=evaluateScores(sig.s(:,1),ymvdr,cfg);
[y_mvdrSuperDirective_pesq, y_mvdrSuperDirectiv_fwsegsnr, y_mvdrSuperDirectiv_ASR]=evaluateScores(sig.s(:,1),ymvdr_superDirective,cfg);
[y_mwfFreeField_pesq, y_mwfFreeField_fwsegsnr, y_mwfFreeField_ASR]=evaluateScores(sig.s(:,1),ymwf_freeField,cfg);
%[S_pesq, S_fwsegsnr, S_ASR]=evaluateScores(sig,Smwf,cfg);
[y_mwf_pesq, y_mwf_fwsegsnr, y_mwf_ASR]=evaluateScores(sig.s(:,1),ymwf,cfg);

[pesqDSB, fwsegsnrDSB, asrDSB]=evaluateScores(sig.s(:,1),yDSB,cfg);
[pesqInput, fwsegsnrInput, asrInput]=evaluateScores(sig.s(:,1),yinput,cfg);

fprintf('\n');
fprintf('\n');
fprintf('                    |----PESQ scores----|-----FwSegSNR----|----ASR----------------\n');
fprintf('           RFSB     |      %g        |      %.2f      |    %g             \n',y_rfsb_pesq,y_rfsb_fwsegsnr,y_rfsb_ASR);
fprintf('           DSB      |      %g        |      %.2f      |    %g             \n',pesqDSB,fwsegsnrDSB,asrDSB);
fprintf('           MVDR     |      %g        |      %.2f      |    %g             \n',y_mvdr_pesq,y_mvdr_fwsegsnr,y_mvdr_ASR);
fprintf('MVDR_SuperDirective |      %g        |      %.2f      |    %g             \n',y_mvdrSuperDirective_pesq,y_mvdrSuperDirectiv_fwsegsnr,y_mvdrSuperDirectiv_ASR);
fprintf('MVDR+MWF_FreeField  |      %g        |      %.2f      |    %g             \n',y_mwfFreeField_pesq,y_mwfFreeField_fwsegsnr,y_mwfFreeField_ASR);
fprintf('   MVDR +  MWF      |      %g        |      %.2f      |    %g             \n',y_mwf_pesq,y_mwf_fwsegsnr,y_mwf_ASR);
fprintf('   Unprocessed      |      %g        |      %.2f      |    %g             \n',pesqInput,fwsegsnrInput,asrInput);
fprintf('----------------------------------------------------------------------------------\n');
