function RunSimulation(cfg,sig,flt)


%cfg.look_azimuth = 0; %azimuth angle

%cfg.nsrc = 1;       % Number of sources
%cfg.nmic = 5;       % Number of mics
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 9 29]; %position of sources, interferers at 40° and 140°
%             cfg.position = [find(cfg.theta_vec==cfg.look_azimuth) 10 32]; %position of sources, interferers at 45° and 155°
%cfg.position = [1 7 22]; %position of sources, interferers at 30
%source must be first defined and than interferer e.x. cfg.position(1) must
%be the desired source
%cfg.noise_type = 0;
% noise implemented in SetAcousticScenario.m
    % 0 = no noise
    % 1 = white noise
    % 2 = noise from files
    % 3 = generated diffuse noise
    % Choose 0 if the whole signal should be used (if ASR scores need to be evaluated!)
%------------------------------------------------------------------
% Set Acoustical scenario parameters (parameters are stored in cfg
% structure
%------------------------------------------------------------------

%------------------------------------------------------------------
% Create microphone signals (signals are stored in sig structure, RIRs in
% flt structure
%------------------------------------------------------------------

%------------------------------------------------------------------
%perform beamformer design (of robust FSB)
%------------------------------------------------------------------

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

[d_hrtf, d_rir] = GetConstants(cfg,flt);


% micPosition = [cfg.mic_pos.x; cfg.mic_pos.y; cfg.mic_pos.z].';
% cfg.micSpacing = micSpacingCalculation(micPosition);
% 
% %EQ 2.21 Optimum Array Processing et. L. Van Trees
% tau = micPosition*[sind(cfg.look_elevation)*cosd(cfg.look_azimuth); sind(cfg.look_elevation)*sind(cfg.look_azimuth); cosd(cfg.look_elevation)];
for idx_mic=1:cfg.nmic
 X(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.x(:,idx_mic),cfg.N,cfg.K,cfg.p);
%  if cfg.nsrc >1
%     X_int(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.xSrc(:,idx_mic,2:end),cfg.N,cfg.K,cfg.p);
%  end
end

% Estimating the PSDs with recursive averaging
% counter=0;
% for m = 1:(cfg.nmic-1)
%     for n = (m+1):cfg.nmic
%         counter = counter + 1;
%         %Thi_Ycpsd(:,:,counter) = Recursive_PSD_estimation(X_int(:,:,m),X_int(:,:,n),0.76);
%         Thi_Ycpsd(:,:,counter) = estimate_cpsd(X(:,:,m),X(:,:,n),0.76);
%     end
% end
% Thi_Ypsd = zeros(size(X,1),size(X,2),size(X,3));
% for m=1:cfg.nmic
%     Thi_Ypsd(:,:,m) = estimate_psd(X(:,:,m),0.1);
%     %Thi_Ypsd(:,:,m) = Recursive_PSD_estimation(X(:,:,m),0.68);
% end
% End of PSD estimation
% 
    %[yrfsb,ydsb, ymvdr,ymvdr_superDirective,ymwf_freeField,ymwf]= frequencyDomain2(X,cfg,d_hrtf,d_rir,flt,tau);
    tau=0;
[yrfsb, ymvdr, ydsb] = frequencyDomain2(X,cfg,d_hrtf,d_rir,flt,tau);
% CDR = zeros(size(X,1),size(X,2));
% Pxx = estimate_psd(X,cfg.alpha);
% counter=0;
% for m = 1:(cfg.nmic-1)
%     for n = (m+1):cfg.nmic
%         counter = counter + 1;
%         TDOA = abs(finddelay(sig.x(:,m),sig.x(:,n)))/cfg.fs;
%         Cxx = estimate_cpsd(X(:,:,m),X(:,:,n),cfg.alpha)./sqrt(Pxx(:,:,m).*Pxx(:,:,n));
%         CDR = CDR + postFilterCDR(Cxx,cfg.micSpacing(m,n),TDOA,cfg);
%     end
% end
% CDR = CDR./counter;
% mu = 1.3;     % noise overestimation factor
% Gmin = 0.1; % minimum Gain
% W = max(1 - (mu./(CDR + 1)).^0.5, 0).^2;
% W = max(W,0);
% W = max(sqrt(W),Gmin);
% W = min(W,1);
% Y = sqrt(mean(abs(X).^2,3)) .* exp(1j*angle(X(:,:,cfg.ref)));
% Ycdr = W.*Y;
%ycdr = DFTSynRealEntireSignal(Ycdr, cfg.N, cfg.K, cfg.p);
%Yinput = X(:,:,cfg.ref);
%yinput = DFTSynRealEntireSignal(Yinput, cfg.N, cfg.K, cfg.p);
% %ynull = DFTSynRealEntireSignal(Ynull, cfg.N, cfg.K, cfg.p);
% 
% 
% % % ------------------------------------------------------------------
% % % BEAMPATTERN for MVDR and RFSB
% % % beampattern(steerV,W,flt.w.RFSB,cfg.angRange.azimuth,cfg.frange)
% % % 
% % % ------------------------------------------------------------------
% % 
% % 
% % % y_sv = sim_system(sig.x(:,:),0,cfg.look_azimuth,'SDB',-10,'GMCC',0,0,0.68,micPosition,512,4,16000,200,7600);
% % 
% % % ------------------------------------------------------------------
% % % Evaluating the fwsegsnr, PESQ and ASR scores
% % % ------------------------------------------------------------------
% % 
[y_rfsb_pesq, y_rfsb_fwsegsnr, y_rfsb_ASR]=evaluateScores(sig.x(:,cfg.ref),yrfsb,cfg,sig);
[pesqDSB, fwsegsnrDSB, asrDSB]=evaluateScores(sig.x(:,cfg.ref),ydsb,cfg,sig);
[y_mvdr_HRTF_pesq, y_mvdr_HRTF_fwsegsnr, y_mvdr_HRTF_ASR]=evaluateScores(sig.x(:,cfg.ref),ymvdr,cfg,sig);
% [y_mvdr_Cyldircal_pesq, y_mvdr_Cyldircal_fwsegsnr, y_mvdr_Cyldircal_ASR]=evaluateScores(sig.x(:,cfg.ref),ymvdr_Cyldircal,cfg,sig);
% [y_mvdr_Diffuse_pesq, y_mvdr_Diffuse_fwsegsnr, y_mvdr_Diffuse_ASR]=evaluateScores(sig.x(:,cfg.ref),ymvdr_Diffuse,cfg,sig);


%[y_mvdrSuperDirective_pesq, y_mvdrSuperDirectiv_fwsegsnr, y_mvdrSuperDirectiv_ASR]=evaluateScores(sig.x(:,cfg.ref),ymvdr_superDirective,cfg,sig);
%[y_mwfFreeField_pesq, y_mwfFreeField_fwsegsnr, y_mwfFreeField_ASR]=evaluateScores(sig.x(:,cfg.ref),ymwf_freeField,cfg,sig);
%[y_mwf_pesq, y_mwf_fwsegsnr, y_mwf_ASR]=evaluateScores(sig.x(:,cfg.ref),ymwf,cfg,sig);
%[y_cdr_pesq, y_cdr_fwsegsnr, y_cdr_ASR]=evaluateScores(sig.x(:,cfg.ref),ycdr,cfg,sig);

%[pesqDSB, fwsegsnrDSB, asrDSB]=evaluateScores(sig.x(:,cfg.ref),ydsb,cfg,sig);
%[pesqInput, fwsegsnrInput, asrInput]=evaluateScores(sig.x(:,cfg.ref),yinput,cfg,sig);

fileID = fopen('resultsHRTF.txt','a');

fprintf(fileID,'\n');
if(cfg.nsrc==1)
    fprintf(fileID,'%6s, %g mics, source at %g, Interferer 0, Noise %g \n',cfg.RIRcond,cfg.nmic,cfg.look_azimuth,cfg.noise_type);
elseif(cfg.nsrc==2)
    fprintf(fileID,'%6s, %g mics, source at %g, Interferer %g, Noise %g \n',cfg.RIRcond,cfg.nmic,cfg.look_azimuth,5*cfg.position(2) - 5, cfg.noise_type);
else
    fprintf(fileID,'%6s, %g mics, source at %g, Interferer %g and %g, Noise %g \n',cfg.RIRcond,cfg.nmic,cfg.look_azimuth,5*cfg.position(2:end) - 5, cfg.noise_type);
end
fprintf(fileID,'                    |----PESQ scores----|-----FwSegSNR----|----ASR----------------\n');
fprintf(fileID,'           LSFIB_HRTF     |      %g        |      %.2f      |    %g             \n',y_rfsb_pesq,y_rfsb_fwsegsnr,y_rfsb_ASR);
fprintf(fileID,'           DSB_HRTF      |      %g        |      %.2f      |    %g             \n',pesqDSB,fwsegsnrDSB,asrDSB);
fprintf(fileID,'           MVDR_HRTF    |      %g        |      %.2f      |    %g             \n',y_mvdr_HRTF_pesq,y_mvdr_HRTF_fwsegsnr,y_mvdr_HRTF_ASR);
% fprintf(fileID,'MVDR_Diffuse |      %g        |      %.2f      |    %g             \n',y_mvdr_Diffuse_pesq,y_mvdr_Diffuse_fwsegsnr,y_mvdr_Diffuse_ASR);
% fprintf(fileID,'MVDR_Cylindrical |      %g        |      %.2f      |    %g             \n',y_mvdr_Cyldircal_pesq,y_mvdr_Cyldircal_fwsegsnr,y_mvdr_Cyldircal_ASR);
%fprintf(fileID,'Kuklasinski Freefield  |      %g        |      %.2f      |    %g             \n',y_mwfFreeField_pesq,y_mwfFreeField_fwsegsnr,y_mwfFreeField_ASR);
%fprintf(fileID,'   Kuklasinski      |      %g        |      %.2f      |    %g             \n',y_mwf_pesq,y_mwf_fwsegsnr,y_mwf_ASR);
%fprintf('   CDR Prop 2       |      %g        |      %.2f      |    %g             \n',y_cdr_pesq,y_cdr_fwsegsnr,y_cdr_ASR);
%fprintf(fileID,'   Unprocessed      |      %g        |      %.2f      |    %g             \n',pesqInput,fwsegsnrInput,asrInput);
fprintf(fileID,'----------------------------------------------------------------------------------\n');
fclose(fileID);
end