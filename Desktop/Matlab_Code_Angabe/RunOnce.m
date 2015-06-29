%main function
%------------------------------------------------------------------
% Initialization
%------------------------------------------------------------------
clear;
close all;

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
cfg.design = 'freefield'; %or freefield
cfg.beamFormer = 'MVDR';
cfg.c = 342; %speed of sound in air
cfg.geometry = 1; %linear array
cfg.spacing = 0; % 0 = non-uniform array, 1=uniform array;
cfg.wng_limit_db = -15;
cfg.alpha = 0.68; %smoothing factor
cfg.nmic=5;
cfg.sig_len = 0;

cfg.RIRcond = 'NAO190_1m';
% 'NAO190_2m';
% 'NAO190_4m';
% 'NAO600_1m';
% 'NAO600_2m';
% 'NAO600_4m';

cfg.look_azimuth = 60;
cfg.angRange.azimuth = 0:5:180;
cfg.look_elevation = repmat(90 - atand(0.73/1.1), length(cfg.look_azimuth));
cfg.angRange.elevation = repmat(cfg.look_elevation,size(cfg.angRange.azimuth));
cfg.des_look_dir.azimuth = cfg.look_azimuth;
cfg.des_look_dir.elevation = cfg.look_elevation;
cfg.nsrc=2;
cfg.position =[cfg.look_azimuth/5 +1 10]; %19 source and interferer at 10*5 -5 =45
cfg.noise_type=0;

%------------------------------------------------------------------
%% filterbank initialization
cfg.fs = 16000;
cfg.N = 512; % FFT size
cfg.K = 128; % frame shift
cfg.Lp = 1024; % prototype filter length
%p2=IterLSDesign(cfg.Lp,cfg.N,cfg.K);
load('/filterbank/prototype_K512_N128_Lp1024.mat');
cfg.p = p; clear p;
cfg.frange=0:cfg.fs/cfg.N:cfg.fs/2;
cfg.k_range = 2*pi*cfg.frange/cfg.c;
cfg = SetAcousticScenario(cfg);
[cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);
cfg = BF_Array_Geometry(cfg);
[flt.w.RFSB, cfg, steerV, realWNG_dB] = RobustFSBdes(cfg);
frange_ext = 200:100:8000;
flt.w.RFSB2 = interpolateFrequencies(frange_ext,cfg,flt);

        


[d_hrtf, d_rir] = GetConstants(cfg,flt);

for idx_mic=1:cfg.nmic
    X(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.x(:,idx_mic),cfg.N,cfg.K,cfg.p);
end

if strcmp(cfg.design, 'hrtf')
    tau=0;
else
    micPosition = [cfg.mic_pos.x; cfg.mic_pos.y; cfg.mic_pos.z].';
    cfg.micSpacing = micSpacingCalculation(micPosition);
    tau = micPosition*[sind(cfg.look_elevation)*cosd(cfg.look_azimuth); sind(cfg.look_elevation)*sind(cfg.look_azimuth); cosd(cfg.look_elevation)];
end
Ydsb=zeros(size(X,1),size(X,2));
Yrfsb=zeros(size(X,1),size(X,2));
Ymvdr=zeros(size(X,1),size(X,2));
Ymwf=zeros(size(X,1),size(X,2));

for k=10:10:45800
    if k==10
    [Yrfsb(:,1:k),Ydsb(:,1:k), Ymvdr(:,1:k),Ymwf(:,1:k)] = frequencyDomain(X(:,1:k,:),cfg,d_hrtf,d_rir,flt,tau);
    else
        if k==45800
            [Yrfsb(:,k-10:45799),Ydsb(:,k-10:45799), Ymvdr(:,k-10:45799),Ymwf(:,k-10:45799)] = frequencyDomain(X(:,k-10:45799,:),cfg,d_hrtf,d_rir,flt,tau);
        else
            [Yrfsb(:,k-10:k),Ydsb(:,k-10:k), Ymvdr(:,k-10:k),Ymwf(:,k-10:k)] = frequencyDomain(X(:,k-10:k,:),cfg,d_hrtf,d_rir,flt,tau);
        end
    end
end

ymvdr = DFTSynRealEntireSignal(Ymvdr, cfg.N, cfg.K, cfg.p);
ymwf = DFTSynRealEntireSignal(Ymwf, cfg.N, cfg.K, cfg.p);
yrfsb = DFTSynRealEntireSignal(Yrfsb, cfg.N, cfg.K, cfg.p);
ydsb = DFTSynRealEntireSignal(Ydsb, cfg.N, cfg.K, cfg.p);

[y_rfsb_pesq, y_rfsb_fwsegsnr, y_rfsb_ASR]=evaluateScores(sig.x(:,cfg.ref),yrfsb,cfg,sig);
[pesqDSB, fwsegsnrDSB, asrDSB]=evaluateScores(sig.x(:,cfg.ref),ydsb,cfg,sig);
[y_mvdr_pesq, y_mvdr_fwsegsnr, y_mvdr_ASR]=evaluateScores(sig.x(:,cfg.ref),ymvdr,cfg,sig);
[y_mwf_pesq, y_mwf_fwsegsnr, y_mwf_ASR]=evaluateScores(sig.x(:,cfg.ref),ymwf,cfg,sig);

fprintf('                    |----PESQ scores----|-----FwSegSNR----|----ASR----------------\n');
fprintf('           LSFIB    |      %g        |      %.2f      |    %g             \n',y_rfsb_pesq,y_rfsb_fwsegsnr,y_rfsb_ASR);
fprintf('           DSB      |      %g        |      %.2f      |    %g             \n',pesqDSB,fwsegsnrDSB,asrDSB);
fprintf('           MVDR     |      %g        |      %.2f      |    %g             \n',y_mvdr_pesq,y_mvdr_fwsegsnr,y_mvdr_ASR);
fprintf('           MWF      |      %g        |      %.2f      |    %g             \n',y_mwf_pesq,y_mwf_fwsegsnr,y_mwf_ASR);
fprintf('----------------------------------------------------------------------------------\n');
