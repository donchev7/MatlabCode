function [Ydsb, Yrfsb, Ymvdr, Ymwf]= frequencyDomain(X,cfg,d_hrtf,d_rir,flt,tau)
mue = 10^(-cfg.wng_limit_db/10); %mue is used for regularization
for idx_freq=1:length(cfg.frange)
    if strcmp(cfg.design, 'freefield')
        beta = (2*pi*cfg.frange(idx_freq))/cfg.c;
        Gamma_tmp_Cylindrical = besselj(0,beta*cfg.micSpacing); %2D cylindrical noise model
        Gamma_tmp_diffuse = sinc(beta*cfg.micSpacing);
        if strcmp(cfg.beamFormer,'MVDR')
            Gamma_const_cylindrical = tril(Gamma_tmp_Cylindrical,-1)./(1+mue) + diag(diag(Gamma_tmp_Cylindrical)) + ...
            triu(Gamma_tmp_Cylindrical,1)./(1+mue);
            Gamma_const_diffuse = tril(Gamma_tmp_diffuse,-1)./(1+mue) + diag(diag(Gamma_tmp_diffuse)) + ...
            triu(Gamma_tmp_diffuse,1)./(1+mue);
        else
            Gamma_const_cylindrical = Gamma_tmp_Cylindrical;
            Gamma_const_diffuse = Gamma_tmp_diffuse;
        end
        GammaCylindrical = Gamma_const_cylindrical;
        GammaDiffuse = Gamma_const_diffuse;
        %EQ 2.25 Optimum Array Processing et. L. Van Trees
        k_vec = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tau;
        %v_k EQ 2.28 Optimum Array Processing et. L. Van Trees
        d0 = exp(-1i*k_vec);
    end
    
    d=squeeze(d_rir(idx_freq,cfg.position(1),:));
    d=d/d(cfg.ref);
    
    gammaIso = covarianceEstimate(squeeze(d_hrtf(idx_freq,:,:)));
    for idx=1:cfg.nmic
        gammaIso(idx,idx)=gammaIso(idx,idx)/max(gammaIso(idx,idx));
    end
    
    Pxx = covarianceEstimate(squeeze(X(idx_freq,:,:)));
    
    if strcmp(cfg.design, 'freefield')
        Wmvdr_Cylindrical = myMVDR2(GammaCylindrical,d0,cfg.beamFormer,cfg.nmic);
        Wmvdr_Diffuse = myMVDR2(GammaDiffuse,d0,cfg.beamFormer,cfg.nmic);
        
        %Wmwf = the multichannel Wiener filter EQ 5b in Kuklasinski paper
        %but for freefield parameters
        Wmwf_freeField = postFilter(Wmvdr_Cylindrical,GammaCylindrical,d0,Pxx,cfg.nmic);
        
        
        Ydsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d0/cfg.nmic);
%         Ymvdr_Cylindrical(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_Cylindrical;
%         Ymvdr_Diffuse(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_Diffuse;
        Ymvdr(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_Diffuse;
        Ymwf(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_Cylindrical*Wmwf_freeField;
    else
        Wmvdr_HRTF = myMVDR2(gammaIso,d,cfg.beamFormer,cfg.nmic);
        %Wmwf as defined in Kuklasinski et. al EQ 5b 
        Wmwf = postFilter(Wmvdr_HRTF,gammaIso,d,Pxx,cfg.nmic);
        Ymvdr(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_HRTF;
        Ymwf(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_HRTF*Wmwf;
        Ydsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d/cfg.nmic);
    end
    Yrfsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*flt.w.RFSB2(:,idx_freq);
end
end