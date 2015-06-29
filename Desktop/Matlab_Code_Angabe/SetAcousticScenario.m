function [cfg] = SetAcousticScenario(cfg)
%--------------------------------------------------------------------------
%    TRINICON-based generic GSC relaization for reverberant conditions
%
%   Klaus Reindl (reindl@LNT.de)                         
%   University of Erlangen-Nuremberg                         
%   Chair of Multimedia Communications and Signal Processing 
% 
%   Date: 03.04.2013   
%--------------------------------------------------------------------------
%SETGENERALPARAMS Set parameters for GSC simulation.
%   [CFG] = SetParams(CFG) sets some parameters and stores them in the
%   structure CFG.
%
%   CFG serves as both input and output. A (typically empty) structure is provided as input argument. 
%   The parameters are stored as structure fields and given back as output argument. 
%   If the input argument CFG is not empty, the structure fields are overwritten. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Acoustical scenario parameters
%--------------------------------------------------------------------------
cfg.RIRtype = 'measured';  % synthetic -> simulated RIR; measured -> real RIR
cfg.fs = 16000;     % Sampling frequency for audio acquisition and playback


%--------------------------------------------------------------------------
% settings for the environmental conditions
%--------------------------------------------------------------------------
if strcmp(cfg.RIRtype, 'synthetic')
%     % parameter settings for synthetic RIRs
%     cfg.synth_room.en = 1;        % Synthetic room impulse responses
%     cfg.synth_room.dim = [2.5, 3, 2.5];  % room dimensions [x, y, z]
%     cfg.synth_room.t60 = 0.15;       % reverberation time
%     cfg.synth_room.order = 1;        % reflections order of RIRs
%     cfg.synth_room.Nh = 80;         % length of RIRs
%     cfg.synth_room.mloc_center = [1, 1, 1.4];
%     cfg.synth_room.mspacing = 0.2;
%     cfg.synth_room.sloc = [0.7, 70;...         %sources location: distance, angle (sources height is fixed to 1.7m
%                            0.7, -70;...
%                            0.7, 30]; 
%     cfg.synth_room.src_alloc = [1, 2, 0]; % 0- none, 1- female, 2- male, 3- white gaussian noise
%     cfg.synth_room.src_paths = {'.\Data\SourceSigs\E205A.WAV',...
%                                 '.\Data\SourceSigs\E202A.WAV',...
%                                 '.\Data\SourceSigs\E208A.WAV',...
%                                 '.\Data\SourceSigs\whitenoise3.wav'};
%     cfg.synth_room.nsrc = size(cfg.synth_room.src_paths, 2);
%     cfg.synth_room.src_gain = [0, 0, 0, 0]; % sources gain [dB]
%     cfg.M = cfg.synth_room.Nh;
%     cfg.startIR = 1;   % cut impulse response from this sample...
%     cfg.endIR = cfg.M; % ... until this sample (0 to take the whole recorded RIR)
% 
%     cfg.synth_room.mloc = repmat(cfg.synth_room.mloc_center, 2, 1)+[1, 0, 0; -1, 0, 0]*cfg.synth_room.mspacing/2;
%     cfg.synth_room.sloc_cart = zeros(size(cfg.synth_room.sloc, 1), 3);
    
% elseif strcmp(cfg.RIRtype, 'measured')
%     cfg.RIRcond = 'NAO190_1m'; %NAO190_1m -> 190ms,1m, NAO190_2m -> 190ms,2m, NAO190_4m -> 190ms,4m 
%                                 %NAO600_1m -> 600ms,1m, NAO600_2m -> 600ms,2m, NAO600_4m -> 600ms,4m
% else
%     error('Invalid choice of the variable cfg.RIRtype\n'); 
 end

%--------------------------------------------------------------------------
% Get the source signals
%   The first mic signal will always be considered as being the primary
%   signal to clean up. 
%   Other mic signals constitute the interference reference signals.
%--------------------------------------------------------------------------
cfg.input_type = 2; 
% 0 = Read mic sig from files
% 1 = Read individual microphone source components
% 2 = Read source signals and convolve them with impulse responses
switch cfg.input_type

    case 1 % Read individual microphone source components

        cfg.path_source_mic  = {'~/Samsung/2011/recordings/2011-03-16/components/ch1-male/01-16k.wav','~/Samsung/2011/recordings/2011-03-16/components/ch1-male/03-16k.wav';...
                                '~/Samsung/2011/recordings/2011-03-16/components/ch2-female/01-16k.wav','~/Samsung/2011/recordings/2011-03-16/components/ch2-female/03-16k.wav';};

    case 2 % Read speech source signals and convolve with impulse responses


        cfg.path_source  = {'./Data/SourceSigs/test_200.wav',...
                            './Data/SourceSigs/E202A.WAV', ...
                            './Data/SourceSigs/E205A.WAV', ...
                            './Data/SourceSigs/E208A.WAV'};
        %cfg.srcgain = sqrt(10.^([40; 40; 40; 40; 40]./10));
        cfg.srcgain = 0;
        % In the specified file a variable imp_resp(sample,mic_positions,source_positions) has to exist!!!
        % switch between measured and synthetic RIRs
        
        %BOBBY - imp_response is a samples x mic_ch x direction of
        %           impulse 
        %           mic_ch is the channel the impulse response was measured
        %           direction of impulse is from 0 degrees to 355 in steps
        %           of 5 (72 total)
        if strcmp(cfg.RIRtype, 'measured')
            cfg.path_imp_resp = ['../Data/RIRs/IR_' cfg.RIRcond '_circle_norm_16kHz.mat'];
            if strcmp(cfg.RIRcond, 'NAO190_1m')
                cfg.path_imp_resp = ['./Data/RIRs/RIR_NAO_AudioLab_az_0_355_T60_190ms_1m_16kHz.mat'];
                cfg.startIR = 1;
                cfg.endIR = 2000;
            elseif strcmp(cfg.RIRcond, 'NAO190_2m')
                cfg.path_imp_resp = ['./Data/RIRs/RIR_NAO_AudioLab_az_0_355_T60_190ms_2m_16kHz.mat'];
                cfg.startIR = 1;
                cfg.endIR = 2000;
            elseif strcmp(cfg.RIRcond, 'NAO190_4m')
                cfg.path_imp_resp = ['./Data/RIRs/RIR_NAO_AudioLab_az_0_355_T60_190ms_4m_16kHz.mat'];
                cfg.startIR = 1;
                cfg.endIR = 2000;                
            elseif strcmp(cfg.RIRcond, 'NAO600_1m')
                cfg.path_imp_resp = ['./Data/RIRs/RIR_NAO_AudioLab_az_0_355_T60_600ms_1m_16kHz.mat'];
                cfg.startIR = 1;
                cfg.endIR = 5461;
            elseif strcmp(cfg.RIRcond, 'NAO600_2m')
                cfg.path_imp_resp = ['./Data/RIRs/RIR_NAO_AudioLab_az_0_355_T60_600ms_2m_16kHz.mat'];
                cfg.startIR = 1;
                cfg.endIR = 5461;
            elseif strcmp(cfg.RIRcond, 'NAO600_4m')
                cfg.path_imp_resp = ['./Data/RIRs/RIR_NAO_AudioLab_az_0_355_T60_600ms_4m_16kHz.mat'];
                cfg.startIR = 1;
                cfg.endIR = 5461;
            else
                error('no such room conditions cfg.RIRcond\n');
            end
        end
        cfg.M = cfg.endIR - cfg.startIR + 1;
        
        %define parameters describing the experimental setup
        % ---------------------------------------------------------
        cfg.theta_vec=(0:5:355);
        %define source positions            

        %define microphone channels
        switch cfg.nmic
            case 3
                cfg.mic_ch = [1 5 9];
            case 5
                cfg.mic_ch = [1 3 5 7 9];
            case 7
                cfg.mic_ch = [1 2 4 5 6 8 9];
            case 9
                cfg.mic_ch = (1:9);
            otherwise
                error('Current number of microphones not supported')
        end
               
        % select reference microphone channels (always the middle one) 
        %%What exactly is the reference microphone channel?
        cfg.ref = ceil(cfg.nmic/2);
        % -----------------------------------------------------------------            
end % switch cfg.input_type

%signal length in seconds
%cfg.sig_len = 10;    % Choose 0 if the whole signal should be used (if ASR scores need to be evaluated!)
%, otherwise set the length in seconds
    
% set source activity (t in [s])
cfg.activity = [0; 0; 0; 0; 0]; 

% Add some noise
%cfg.noise_type = 0;
    % 0 = no noise
    % 1 = white noise
    % 2 = noise from files
    % 3 = generated diffuse noise
switch cfg.noise_type
    case 0 % no noise    
    case 1 % white noise
        cfg.inputsnr = 30; % input SNR in dB
        
    case 2 % noise from files
        cfg.inputsnr = 25; % input SNR in dB
        cfg.file_noise = {'../Data/SourceSigs/noise1'...
                          '../Data/SourceSigs/noise2'...
                          '../Data/SourceSigs/noise3'}; 
    case 3
        Error('Not implemented/supported yet.\n')
%         cfg.inputsnr = 20; % input SNR in dB
%         cfg.file_noise = '../Data/SourceSigs/DiffuseNoise_6channel-ULA_16kHz.mat';        
    otherwise
        error('Invalid choice of the variable cfg.noise_type\n');
end