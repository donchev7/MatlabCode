function [Ydsb,Yrfsb, Ymvdr,Ymvdr_Diffuse,Ymvdr_Cyldircal,Ymwf_freeField, Ymwf]= frequencyDomain(X,cfg,d_hrtf,d_rir,flt,tau)
Ydsb=zeros(size(X,1),size(X,2));
Yrfsb=zeros(size(X,1),size(X,2));
Ymvdr=zeros(size(X,1),size(X,2));
Ymvdr_Diffuse=zeros(size(X,1),size(X,2));
Ymvdr_Cyldircal=zeros(size(X,1),size(X,2));
Ymwf=zeros(size(X,1),size(X,2));
Ymwf_freeField=zeros(size(X,1),size(X,2));
% Ycdr=zeros(size(X,1),size(X,2));
 mue = 10^(-cfg.wng_limit_db/10); %mue is used for regularization
for idx_freq=1:length(cfg.frange)
    if strcmp(cfg.design, 'freefield')
        beta = (2*pi*cfg.frange(idx_freq))/cfg.c;
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
        %EQ 2.25 Optimum Array Processing et. L. Van Trees
        k_vec = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tau;
        %v_k EQ 2.28 Optimum Array Processing et. L. Van Trees
        d0 = exp(-1i*k_vec);
    end
    gammaIso = covarianceEstimate(squeeze(d_hrtf(idx_freq,:,:)));
    d=squeeze(d_rir(idx_freq,cfg.position(1),:));
    
    d=d/d(cfg.ref);
    for idx=1:cfg.nmic
        gammaIso(idx,idx)=gammaIso(idx,idx)/max(gammaIso(idx,idx));
    end
    Pxx = covarianceEstimate(squeeze(X(idx_freq,:,:)));
    
    Wmvdr_freeField = myMVDR2(GammaCyldircal,d0,cfg.beamFormer,cfg.nmic);
    Wmvdr_Diffuse = myMVDR2(GammaDiffuse,d,cfg.beamFormer,cfg.nmic);
    Wmvdr_GammaCyldrical = myMVDR2(GammaCyldircal,d,cfg.beamFormer,cfg.nmic);
    Wmvdr_superDirective = myMVDR2(GammaDiffuse,d0,cfg.beamFormer,cfg.nmic);
    %Wmvdr_superDirective = filter for freefield 2D cyldirical noise and freefield steering vector d0
    WgammaIso = myMVDR2(gammaIso,d,cfg.beamFormer,cfg.nmic);
    %Wmwf = the multichannel Wiener filter EQ 5b in Kuklasinski paper
    %but for freefield parameters
    Wmwf_freeField = postFilter(Wmvdr_freeField,GammaCyldircal,d0,Pxx,cfg.nmic);

    %Wmwf as defined in Kuklasinski et. al EQ 5b 
    Wmwf = postFilter(WgammaIso,gammaIso,d,Pxx,cfg.nmic);
    %Conventional MVDR beamformer with covariance matrix of microphone signals as input

    Wmvdr = myMVDR2(gammaIso,d,cfg.beamFormer,cfg.nmic);
    
    Ydsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d0/cfg.nmic);
    Yrfsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*flt.w.RFSB2(:,idx_freq);
    Ymvdr(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d0/cfg.nmic);
    Ymvdr_Diffuse(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d0/cfg.nmic);
    Ymvdr_Cyldircal(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d0/cfg.nmic);
    Ymwf_freeField(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf_freeField;
    Ymwf(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf;
end
end