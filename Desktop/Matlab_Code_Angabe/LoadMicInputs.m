function [cfg,sig,flt] = LoadMicInputs(cfg,sig,flt)
%---------------------------------------------------------
%Blind Source Separation (BSS) and Localization
%
%   Anthony Lombard (lombard@LNT.de)                         
%   University of Erlangen-Nuremberg                         
%   Chair of Multimedia Communications and Signal Processing 
% 
%   Date: 06.07.2007   
%---------------------------------------------------------
%LOADMICINPUTS Load the microphone signals for BSS simulation.
%   [CFG,SIG,FLT]=LoadMicInputs(CFG,SIG,FLT) uses the parameters previously set 
%   using the functions SetParams and SetExtraParams of the BSS simulation software. 
%
%   Microphone signals can be directly read from WAV files or created 
%   by convolving clean source signals with measured room impulse responses.
%   Background noise can also be added to simulate noisy environments.
%   The signals are stored in the SIG structure, while the room impulse responses 
%   are stored in the structure FLT if required. 
%   Additional parameters are also set while loading the microphone signals 
%   like, e.g., the effective signal length, and stored in the structure CFG.
%
%   CFG         -> structure of parameters   
%   SIG         -> structure of signals
%   SIG.x       -> microphone signals
%   SIG.xsrc    -> containing all cfg.nsrc components at the microphones 
%                  structure: xsrc(samples,cfg.nmic,cfg.nsrc)
%   SIG.xclean  -> clean microphone signals
%   SIG.xnoise  -> noise at microphone signals
%   SIG.s       -> source signals
%   FLT         -> structure of filters
%   FLT.h       -> mixing matrix h(samples,cfg.nsrc,cfg.nmic)



s = [];
x = [];
xsrc = [];  
xclean = [];
xnoise = [];
h = [];

%-------------------------------------------------------
% Read the microphone signals
%-------------------------------------------------------        
% fprintf('\t* Read the source signals from wav-files and convolve them with RIRs\n');

if strcmp(cfg.RIRtype, 'synthetic')
    %Not necessary for this thesis
    
elseif strcmp(cfg.RIRtype, 'measured')
    % Read source signals
    for q = 1:cfg.nsrc
        [tmp,f] = audioread(cfg.path_source{q});
        tmp = resample(tmp,cfg.fs,f);
        if cfg.sig_len ~= 0
            while length(tmp)<(cfg.sig_len*cfg.fs)
                tmp = cat(1,tmp,tmp);
            end
            s(:,q) = tmp(1:cfg.sig_len*cfg.fs);
        else
            if q == 1
                s(:,q)=tmp;
            else
                s = s(1:min(length(tmp),length(s)),:);
                s(:,q)=tmp(1:min(length(tmp),length(s)));
            end
            cfg.sig_len = length(s)/cfg.fs;
        end
        %normalize to amplitude 1
        %why normalize to amplitude 1?
        s(:,q) = s(:,q)/(1.1*max(abs(s(:,q))));
    end
    
    % Get inputs of equal power
    %why? because if one speaker speakes louder we cant judge the
    %peroformance of the interference cancellation
    s = s./repmat(sqrt(var(s)/min(var(s))),size(s,1),1);
    
    % Apply a gain on the source signals
    % BOBBY - distorts the signal!!
%     if length(find(cfg.srcgain~=1))>0
% %         fprintf('\t* Apply a gain on the source signals\n');
%         for q = 1:cfg.nsrc
%             s(:,q) = s(:,q)*cfg.srcgain(q);
%         end
%     end
    
    % Convolve with impulse responses
%     fprintf('\t* Position the sources...\n')
    load(cfg.path_imp_resp);
    if exist('fs_RIR')
        if cfg.fs~=fs_RIR
            imp_resp_tmp = imp_resp;
            imp_resp = [];
            for q = 1:size(imp_resp_tmp,3)
                imp_resp(:,:,q) = resample(imp_resp_tmp(:,:,q),cfg.fs,fs_RIR);
            end
            clear imp_resp_tmp
        end
    else
%         fprintf('*************************************************************************\n');
%         fprintf('*** Are you sure you use RIRs recorded at the sampling rate %d Hz? ***\n',cfg.fs);
%         fprintf('*************************************************************************\n');
    end
    for p = 1:cfg.nmic
        hTemp = squeeze(imp_resp(:,cfg.mic_ch(p),:));
        h(:,:,p) = hTemp(cfg.startIR:cfg.startIR+cfg.endIR-1,:);
    end
    %normalize impulse responses to absolut maximum of all required impulse
    %responses from desired source direction
    %    
    h = h/max(max(abs(h(:,cfg.position(1),:))));    
    
    xsrc = zeros(size(s,1)+cfg.endIR-1,cfg.nmic,cfg.nsrc); % structure: xsrc(samples,mics,#speaker)
    for q = 1:cfg.nsrc
        for p = 1:cfg.nmic
%             hTemp = squeeze(imp_resp(:,cfg.mic_ch(p),cfg.position(q)));
%             hTemp = hTemp(cfg.startIR:cfg.startIR+cfg.endIR-1);           
            hTemp = h(:,cfg.position(q),p);
            xsrc(:,p,q) = fftfilt(hTemp,[s(:,q); zeros(length(hTemp)-1,1)]);
        end
    end
            
    % Sum up the individual source components to obtain the mic signals
    x = sum(xsrc,3);
    
end  %if strcmp(cfg.RIRtype, 'synthetic')     

%--------------------------------------------------------------------------
% Add some additional noise to the microphones
%--------------------------------------------------------------------------
xclean = x;

switch cfg.noise_type
    
    case 0 % Add no oise  
        
        xnoise = zeros(size(x));
    
    case 1 % Add white noise to mic
        %cfg.inputsnr defined in SetAcousticScenario.m
        fprintf('\t* Add some white noise with SNR %d dB\n',cfg.inputsnr);
        
        xnoise = 2*rand(size(x))-1;
        % adjust SNR to desired value (same scale factor for both channels)
        scalefac = sqrt( mean( mean(sum(x.^2,1)) ./ mean(sum(xnoise.^2,1)) ) ./ 10.^(cfg.inputsnr/10) );
        xnoise = xnoise.*scalefac;
        x = x + xnoise;

    case 2 % Add noise from files to mic
        
        fprintf('\t* Add noise from files with SNR %d dB\n',cfg.inputsnr);            
        
        for p=1:cfg.nmic
            [tmp_noise,noise_rate] = wavread(cfg.file_noise{p});  
            tmp_noise = resample(tmp_noise,cfg.fs,noise_rate);        
            if p==1
                xnoise(:,p) = tmp_noise;
            else
                xnoise = xnoise(1:min(length(tmp_noise),length(xnoise)),:);
                xnoise(:,p)=tmp_noise(1:min(length(tmp_noise),length(xnoise)));            
            end
        end    
        while size(xnoise,1)<size(x,1)
            xnoise = cat(1,xnoise,xnoise);
        end
        xnoise = xnoise(1:size(x,1),:);  
        % adjust SNR to desired value (same scale factor for all channels)
        scalefac = sqrt( mean( mean(sum(x.^2,1)) ./ mean(sum(xnoise.^2,1)) ) ./ 10.^(cfg.inputsnr/10) );
        xnoise = xnoise.*scalefac;
        x = x + xnoise;
        
    case 3
        tmp_noise = load(cfg.file_noise);
        %get required channels
        xnoise = tmp_noise.z([cfg.mic_ch],:).'; 
        while size(xnoise,1)<size(x,1)
            xnoise = cat(1,xnoise,tmp_noise.z([cfg.mic_ch],:).');
        end
        xnoise = xnoise(1:size(x,1),:); 
        
        % adjust SNR to desired value (same scale factor for all channels)
        scalefac = sqrt( mean( mean(sum(x.^2,1)) ./ mean(sum(xnoise.^2,1)) ) ./ 10.^(cfg.inputsnr/10) );
        xnoise = xnoise.*scalefac;
        x = x + xnoise;
        
end %switch cfg.noise_type

%--------------------------------------------------------------------------
% Output variables
%--------------------------------------------------------------------------
sig.s = s;
sig.xSrc = xsrc;
sig.xClean = xclean;  
if cfg.noise_type
    sig.xNoise = xnoise;
end
sig.s = s; 
sig.x = x;
flt.h = h;
%--------------------------------------------------------------------------