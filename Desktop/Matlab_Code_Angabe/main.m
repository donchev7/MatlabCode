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
cfg.geometry = 1; %linear array
cfg.spacing = 0; % 0 = non-uniform array, 1=uniform array;
% cfg.design = 'hrtf';
cfg.wng_limit_db = -15;
cfg.look_azimuth = 45; %azimuth angle
cfg.look_elevation = repmat(90 - atand(0.73/1.1), length(cfg.look_azimuth)); %elevation
cfg.des_look_dir.azimuth = cfg.look_azimuth;
cfg.des_look_dir.elevation = cfg.look_elevation;
cfg.nsrc = 2;       % Number of sources
cfg.nmic = 5;       % Number of mics
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 9 29]; %position of sources, interferers at 40° and 140°
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 10 32]; %position of sources, interferers at 45° and 155°
cfg.position = [10 46 22]; %position of sources, interferers at 45° and 105°
%source must be first defined and than interferer e.x. cfg.position(1) must
%be the desired source
cfg.noise_type = 0;
% noise implemented in SetAcousticScenario.m
    % 0 = no noise
    % 1 = white noise
    % 2 = noise from files
    % 3 = generated diffuse noise
    
%------------------------------------------------------------------
% Set Acoustical scenario parameters (parameters are stored in cfg
% structure
%------------------------------------------------------------------
[cfg] = SetAcousticScenario(cfg);

%------------------------------------------------------------------
% Create microphone signals (signals are stored in sig structure, RIRs in
% flt structure
%------------------------------------------------------------------
[cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);

%------------------------------------------------------------------
%% filterbank initialization
cfg.c = 342; %speed of sound in air
cfg.K = 512; % FFT size
cfg.N = 128; % frame shift
%cfg.Lp = 1024; % prototype filter length
%p=IterLSDesign(cfg.Lp,cfg.K,cfg.N);
%cfg.p=IterLSDesign(cfg.Lp,cfg.K,cfg.N);
%cfg.frange = 0:250:8000;
cfg.frange = linspace(0,cfg.fs/2,cfg.K/2+1)'; % frequency axis
cfg.k_range = 2*pi*cfg.frange/cfg.c;
cfg.angRange.azimuth = 0:5:180;
cfg.angRange.elevation = repmat(cfg.look_elevation,size(cfg.angRange.azimuth));

%------------------------------------------------------------------

%------------------------------------------------------------------
%perform beamformer design (of robust FSB)
%------------------------------------------------------------------
[flt.w.RFSB, cfg, steerV2, realWNG_dB] = RobustFSBdes(cfg);
flt.w.MVDR = myMVDR2(cfg,sig);

steerV = steeringVector(cfg);
for idx_freq=1:length(cfg.frange)
    BPatternRFSB(idx_freq,:) = mag2db(abs(squeeze(steerV(idx_freq,:,:))*flt.w.RFSB(:,idx_freq)));
end
for idx_freq=1:length(cfg.frange)
    BPatternMVDR(idx_freq,:) = mag2db(abs(squeeze(steerV(idx_freq,:,:))*flt.w.MVDR(:,idx_freq)));
end
figure
subplot(1,2,1);
imagesc(cfg.angRange.azimuth,cfg.frange,BPatternMVDR)
xlabel('DOA in degrees')
ylabel('Frequency in Hz')
title('Beampattern MVDR')
colorbar

subplot(1,2,2);
imagesc(cfg.angRange.azimuth,cfg.frange,BPatternRFSB)
xlabel('DOA in degrees')
ylabel('Frequency in Hz')
title('Beampattern RFSB')
colorbar


%------------------------------------------------------------------
% Create microphone signals (signals are stored in sig structure) at
% beamformer output signal
%------------------------------------------------------------------
rown = ceil((1+cfg.K)/2);            % calculate the total number of rows
coln = 1+fix((length(sig.x(:,1)-cfg.K)/cfg.N));   
for idx_mic=1:cfg.nmic
 X(:,:,idx_mic) = stft(sig.x(:,idx_mic),cfg.K,cfg.N,cfg.K,cfg.fs);
end
for idx_freq=1:length(cfg.frange) 
Y(idx_freq,:) = flt.w.RFSB(:,idx_freq)'*squeeze(X(idx_freq,:,:)).';
end
sig.y.RFSB = istft(Y, cfg.N, cfg.K, cfg.fs);


% %------------------------------------------------------------------
% % Potential place for postfilter
% %------------------------------------------------------------------
% 
% %------------------------------------------------------------------
% % Evaluation regarding the fwsegsnr, PESQ score and ASR score
% %------------------------------------------------------------------
% switch cfg.design
%     case 'freefield'
%         [sig, cfg, res.freefield{idx_phi}] = Evaluation(sig, cfg);
%     case 'hrtf'
%         [sig, cfg, res.hrtf{idx_phi}] = Evaluation(sig, cfg);
% end
% 
% %------------------------------------------------------------------
% % save results (if needed)
% %------------------------------------------------------------------
% save(['saves/nameOfResultFile' '.mat'], 'res');