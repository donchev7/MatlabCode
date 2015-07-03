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
cfg.design = 'hrtf';
cfg.beamFormer = 'MVDR';
cfg.geometry = 1; %linear array
cfg.spacing = 0; % 0 = non-uniform array, 1=uniform array;
% cfg.design = 'hrtf';
cfg.wng_limit_db = -15;
cfg.alpha = 0.68; %smoothing factor
cfg.nmic=5;
cfg.sig_len = 0;

c = containers.Map;
c('1') = 'NAO190_1m';
c('2') = 'NAO190_2m';
%c('3') = 'NAO190_4m';
c('3') ='NAO600_1m';
c('4') ='NAO600_2m';
%c('4') ='NAO600_4m';

source=0:30:180;
interferer=[45 165];
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
%cfg.frange = linspace(0,cfg.fs/2,cfg.N/2+1)';
cfg.frange=0:cfg.fs/cfg.N:cfg.fs/2;
cfg.k_range = 2*pi*cfg.frange/cfg.c;
cfg.angRange.azimuth = 0:5:180;
frange_ext = 200:100:8000;

for i=1:length(keys(c))
    cfg.RIRcond = c(num2str(i));
    for j=1:length(source)
        cfg.look_azimuth = source(j);
        cfg.look_elevation = repmat(90 - atand(0.73/1.1), length(cfg.look_azimuth));
        cfg.angRange.elevation = repmat(cfg.look_elevation,size(cfg.angRange.azimuth));
        cfg.des_look_dir.azimuth = cfg.look_azimuth;
        cfg.des_look_dir.elevation = cfg.look_elevation;
        cfg.position =source(j)/5 + 1;
        cfg.nsrc=1;
        cfg.noise_type=0;
        cfg = SetAcousticScenario(cfg);
        [cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);
        cfg = BF_Array_Geometry(cfg);
        [flt.w.RFSB, cfg, steerV, realWNG_dB] = RobustFSBdes(cfg);
        flt.w.RFSB2 = interpolateFrequencies(frange_ext,cfg,flt);
        RunSimulation(cfg,sig,flt);
        for k=1:length(interferer)
            if k==1
                cfg.nsrc = 2; 
                cfg.position =[source(j)/5+1 interferer(k)/5+1];
                [cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);
                RunSimulation(cfg,sig,flt);
            else
                cfg.nsrc = 3; 
                cfg.position =[source(j)/5+1 interferer(k-1)/5+1 interferer(k)/5+1];
                [cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);
                RunSimulation(cfg,sig,flt);
            end
        end
        cfg.nsrc = 1; 
        cfg.noise_type = 1;
        cfg = SetAcousticScenario(cfg);
        [cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);
        RunSimulation(cfg,sig,flt);
        cfg.nsrc = 2;
        cfg.position =[source(j)/5+1 interferer(1)/5+1];
        [cfg,sig,flt] = LoadMicInputs(cfg,sig,flt);
        RunSimulation(cfg,sig,flt);
        fprintf('%s Source %g \n',cfg.RIRcond,source(j));
    end
end

