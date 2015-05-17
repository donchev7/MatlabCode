function [yrfsb,ymvdr,ymvdr_superDirective,ymwf_freeField,ymwf]= frequencyDomain(Xpadd,cfg,mue,d_hrtf,gammaIso,flt,tau,X)
for idx_freq=1:size(Xpadd,1)
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
    GammaCyldircal = Gamma_const_cylindrical;
    GammaDiffuse = Gamma_const_diffuse;

    d=d_hrtf(idx_freq,:,cfg.position(1)).';
    %EQ 2.25 Optimum Array Processing et. L. Van Trees
    k_vec = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tau;
    %v_k EQ 2.28 Optimum Array Processing et. L. Van Trees
    d0 = exp(-1i*k_vec);
    %Wmvdr_superDirective = filter for freefield 2D cyldirical noise and freefield steering vector d0
    Wmvdr_superDirective = myMVDR2(GammaDiffuse,d0,cfg.beamFormer);
    Wmvdr_freeField = myMVDR2(GammaCyldircal,d0,cfg.beamFormer);

    %WgammaIso = filter with evaluated cylindrical sound field from HRTF and
    %hrtf evaluated steering vector d
    WgammaIso = myMVDR2(gammaIso,d,cfg.beamFormer);
    for idx_frame=50:50:size(Xpadd,2)
    
        %X_aligned(idx_freq,:,:) = bsxfun(@times,d0(:,idx_freq)',squeeze(X(idx_freq,:,:)));
        %X_int_aligned(idx_freq,:,:)=bsxfun(@times,d0(:,idx_freq)',squeeze(X_int(idx_freq,:,:)));
        Pxx2 = cov(squeeze(Xpadd(idx_freq,idx_frame-49:idx_frame,:)));
        %Pxx = Thi_Y(squeeze(Thi_Ypsd(idx_freq,:,:)),squeeze(Thi_Ypsd(idx_freq,:,:)),cfg.nmic);
        %Thi_Y(:,:,idx_freq) = cov(squeeze(X(idx_freq,:,:)))./size(X,2);


        %Ymwf2(idx_freq,:) = squeeze(X(idx_freq,:,:))*WgammaIso;
        %S(idx_freq,:)=Ymwf2(idx_freq,:)*postFilter2(WgammaIso,gammaIso,d,Pxx,cfg.nmic);
        %Wmwf = the multichannel Wiener filter EQ 5b in Kuklasinski paper
        %but for freefield parameters
        Wmwf_freeField = postFilter(Wmvdr_freeField,GammaCyldircal,d0,Pxx2,cfg.nmic);

        %Wmwf as defined in Kuklasinski et. al EQ 5b 
        Wmwf = postFilter(WgammaIso,gammaIso,d,Pxx2,cfg.nmic);
        %Conventional MVDR beamformer with covariance matrix of microphone signals as input
        if cfg.nsrc > 1
            Wmvdr = myMVDR2(cov(squeeze(X_int(idx_freq,idx_frame-49:idx_frame,:))),d0,cfg.beamFormer);
        else

            Wmvdr = myMVDR2(Pxx2,d0,cfg.beamFormer);
        end
        Ymwf_freeField(idx_freq,idx_frame-49:idx_frame) = squeeze(X(idx_freq,idx_frame-49:idx_frame,:))*Wmwf_freeField;
        Ymwf(idx_freq,idx_frame-49:idx_frame) = squeeze(X(idx_freq,idx_frame-49:idx_frame,:))*Wmwf;
        Ymvdr(idx_freq,idx_frame-49:idx_frame) = squeeze(X(idx_freq,idx_frame-49:idx_frame,:))*Wmvdr;
        Ymvdr_superDirective(idx_freq,idx_frame-49:idx_frame) = squeeze(X(idx_freq,idx_frame-49:idx_frame,:))*Wmvdr_superDirective;
    end
    Yrfsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*flt.w.RFSB(:,idx_freq);
end
ymvdr = DFTSynRealEntireSignal(Ymvdr, cfg.N, cfg.K, cfg.p);
ymvdr_superDirective = DFTSynRealEntireSignal(Ymvdr_superDirective, cfg.N, cfg.K, cfg.p);
ymwf = DFTSynRealEntireSignal(Ymwf, cfg.N, cfg.K, cfg.p);
ymwf_freeField = DFTSynRealEntireSignal(Ymwf_freeField, cfg.N, cfg.K, cfg.p);
yrfsb = DFTSynRealEntireSignal(Yrfsb, cfg.N, cfg.K, cfg.p);

% [y_rfsb_pesq, y_rfsb_fwsegsnr, y_rfsb_ASR]=evaluateScores(sig,yrfsb,cfg);
% [y_mvdr_pesq, y_mvdr_fwsegsnr, y_mvdr_ASR]=evaluateScores(sig,ymvdr,cfg);
% [y_mvdrSuperDirective_pesq, y_mvdrSuperDirectiv_fwsegsnr, y_mvdrSuperDirectiv_ASR]=evaluateScores(sig,ymvdr_superDirective,cfg);
% [y_mwfFreeField_pesq, y_mwfFreeField_fwsegsnr, y_mwfFreeField_ASR]=evaluateScores(sig,ymwf_freeField,cfg);
% %[S_pesq, S_fwsegsnr, S_ASR]=evaluateScores(sig,Smwf,cfg);
% [y_mwf_pesq, y_mwf_fwsegsnr, y_mwf_ASR]=evaluateScores(sig,ymwf,cfg);
% 
% [pesqMAX, fwsegsnrMAX, asrMAX]=evaluateScores(sig,sig.s(:,1),cfg);
% 
% fprintf('\n');
% fprintf('\n');
% fprintf('                    |----PESQ scores----|-----FwSegSNR----|----ASR----------------\n');
% fprintf('           RFSB     |      %g        |      %.2f      |    %g             \n',y_rfsb_pesq,y_rfsb_fwsegsnr,y_rfsb_ASR);
% fprintf('           MVDR     |      %g        |      %.2f      |    %g             \n',y_mvdr_pesq,y_mvdr_fwsegsnr,y_mvdr_ASR);
% fprintf('MVDR_SuperDirective |      %g        |      %.2f      |    %g             \n',y_mvdrSuperDirective_pesq,y_mvdrSuperDirectiv_fwsegsnr,y_mvdrSuperDirectiv_ASR);
% fprintf('MVDR+MWF_FreeField  |      %g        |      %.2f      |    %g             \n',y_mwfFreeField_pesq,y_mwfFreeField_fwsegsnr,y_mwfFreeField_ASR);
% fprintf('   MVDR +  MWF      |      %g        |      %.2f      |    %g             \n',y_mwf_pesq,y_mwf_fwsegsnr,y_mwf_ASR);
% %fprintf('   S MWF            |      %g        |      %.2f      |    %g             \n',S_pesq,S_fwsegsnr,S_ASR);
% fprintf('   MAX optainable   |      %g        |      %.2f      |    %g             \n',pesqMAX,fwsegsnrMAX,asrMAX);
% fprintf('----------------------------------------------------------------------------------\n');
end