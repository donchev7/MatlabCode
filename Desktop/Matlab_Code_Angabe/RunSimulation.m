function RunSimulation(cfg,sig,flt)


[d_hrtf, d_rir] = GetConstants(cfg,flt);


for idx_mic=1:cfg.nmic
 X(:,:,idx_mic) = DFTAnaRealEntireSignal(sig.x(:,idx_mic),cfg.N,cfg.K,cfg.p);
end

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
Yinput = X(:,:,cfg.ref);
yinput = DFTSynRealEntireSignal(Yinput, cfg.N, cfg.K, cfg.p);

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



[y_rfsb_pesq, y_rfsb_fwsegsnr, y_rfsb_ASR]=evaluateScores(sig.x(:,cfg.ref),yrfsb,cfg,sig);
[pesqDSB, fwsegsnrDSB, asrDSB]=evaluateScores(sig.x(:,cfg.ref),ydsb,cfg,sig);
[y_mvdr_pesq, y_mvdr_fwsegsnr, y_mvdr_ASR]=evaluateScores(sig.x(:,cfg.ref),ymvdr,cfg,sig);
[y_mwf_pesq, y_mwf_fwsegsnr, y_mwf_ASR]=evaluateScores(sig.x(:,cfg.ref),ymwf,cfg,sig);
[pesqInput, fwsegsnrInput, asrInput]=evaluateScores(sig.x(:,cfg.ref),yinput,cfg,sig);

if strcmp(cfg.design, 'hrtf')
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
    fprintf(fileID,'           MVDR_HRTF    |      %g        |      %.2f      |    %g             \n',y_mvdr_pesq,y_mvdr_fwsegsnr,y_mvdr_ASR);
    fprintf(fileID,'   Kuklasinski      |      %g        |      %.2f      |    %g             \n',y_mwf_pesq,y_mwf_fwsegsnr,y_mwf_ASR);
    %fprintf('   CDR Prop 2       |      %g        |      %.2f      |    %g             \n',y_cdr_pesq,y_cdr_fwsegsnr,y_cdr_ASR);
    %fprintf(fileID,'   Unprocessed      |      %g        |      %.2f      |    %g             \n',pesqInput,fwsegsnrInput,asrInput);
    fprintf(fileID,'----------------------------------------------------------------------------------\n');
    fclose(fileID);
else
    fileID = fopen('resultsFreeField.txt','a');
    fprintf(fileID,'\n');
    if(cfg.nsrc==1)
        fprintf(fileID,'%6s, %g mics, source at %g, Interferer 0, Noise %g \n',cfg.RIRcond,cfg.nmic,cfg.look_azimuth,cfg.noise_type);
    elseif(cfg.nsrc==2)
        fprintf(fileID,'%6s, %g mics, source at %g, Interferer %g, Noise %g \n',cfg.RIRcond,cfg.nmic,cfg.look_azimuth,5*cfg.position(2) - 5, cfg.noise_type);
    else
        fprintf(fileID,'%6s, %g mics, source at %g, Interferer %g and %g, Noise %g \n',cfg.RIRcond,cfg.nmic,cfg.look_azimuth,5*cfg.position(2:end) - 5, cfg.noise_type);
    end
    fprintf(fileID,'                    |----PESQ scores----|-----FwSegSNR----|----ASR----------------\n');
    fprintf(fileID,'           LSFIB_FreeField     |      %g        |      %.2f      |    %g             \n',y_rfsb_pesq,y_rfsb_fwsegsnr,y_rfsb_ASR);
    fprintf(fileID,'           DSB_FreeField      |      %g        |      %.2f      |    %g             \n',pesqDSB,fwsegsnrDSB,asrDSB);
    fprintf(fileID,'MVDR_Diffuse |      %g        |      %.2f      |    %g             \n',y_mvdr_pesq,y_mvdr_fwsegsnr,y_mvdr_ASR);
    fprintf(fileID,'Kuklasinski_Freefield  |      %g        |      %.2f      |    %g             \n',y_mwf_pesq,y_mwf_fwsegsnr,y_mwf_ASR);
    fprintf(fileID,'   Unprocessed      |      %g        |      %.2f      |    %g             \n',pesqInput,fwsegsnrInput,asrInput);
    fprintf(fileID,'----------------------------------------------------------------------------------\n');
    fclose(fileID);
end
end