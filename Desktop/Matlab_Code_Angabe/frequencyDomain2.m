function [yrfsb,ydsb, ymvdr,ymvdr_superDirective,ymwf_freeField,ymwf,ycdr]= frequencyDomain2(X,cfg,mue,d_hrtf,d_rir,flt,tau)
% c=1;
Ydsb=zeros(size(X,1),size(X,2));
Yrfsb=zeros(size(X,1),size(X,2));
Ymvdr=zeros(size(X,1),size(X,2));
Ymvdr_superDirective=zeros(size(X,1),size(X,2));
Ymwf=zeros(size(X,1),size(X,2));
Ymwf_freeField=zeros(size(X,1),size(X,2));
%Ycdr=zeros(size(X,1),size(X,2));
for idx_freq=1:length(cfg.frange)
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
    gammaIso = covarianceEstimate(squeeze(d_hrtf(idx_freq,:,:)));
    d=squeeze(d_rir(idx_freq,cfg.position(1),:));
    
    d=d/d(cfg.ref);
    %sd=max(gammaIso);
    for idx=1:cfg.nmic
        gammaIso(idx,idx)=gammaIso(idx,idx)/max(gammaIso(idx,idx));
    end
    %EQ 2.25 Optimum Array Processing et. L. Van Trees
    k_vec = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tau;
    %k_vecInt = (-2*pi*cfg.frange(idx_freq)/cfg.c)*tauInt;
    %v_k EQ 2.28 Optimum Array Processing et. L. Van Trees
    d0 = exp(-1i*k_vec);
   %dint = exp(-1i*k_vecInt);
    %A = [d0, dint];
    %e1=[1,0];
    %wNull = (e1*pinv(A)).';
    %X_aligned(idx_freq,:,:) = bsxfun(@times,d0(:,idx_freq)',squeeze(X(idx_freq,:,:)));
    %X_int_aligned(idx_freq,:,:)=bsxfun(@times,d0(:,idx_freq)',squeeze(X_int(idx_freq,:,:)));
    Pxx = covarianceEstimate(squeeze(X(idx_freq,:,:)));
    
    Wmvdr_freeField = myMVDR2(GammaCyldircal,d0,cfg.beamFormer,cfg.nmic);
    Wmvdr_superDirective = myMVDR2(GammaDiffuse,d0,cfg.beamFormer,cfg.nmic);
    %Wmvdr_superDirective = filter for freefield 2D cyldirical noise and freefield steering vector d0
    WgammaIso = myMVDR2(gammaIso,d,cfg.beamFormer,cfg.nmic);
    %Wmwf = the multichannel Wiener filter EQ 5b in Kuklasinski paper
    %but for freefield parameters
    Wmwf_freeField = postFilter(Wmvdr_freeField,GammaCyldircal,d0,Pxx,cfg.nmic);

    %Wmwf as defined in Kuklasinski et. al EQ 5b 
    Wmwf = postFilter(WgammaIso,gammaIso,d,Pxx,cfg.nmic);
    %Conventional MVDR beamformer with covariance matrix of microphone signals as input
    if cfg.nsrc > 1
        Wmvdr = myMVDR2(covarianceEstimate(squeeze(X_int(idx_freq,:,:))),d0,cfg.beamFormer,cfg.nmic);
    else

        Wmvdr = inv(Pxx)*d0/(d0'*inv(Pxx)*d0);%myMVDR2(Pxx,d0,cfg.beamFormer,cfg.nmic);
    end
    Ymwf_freeField(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf_freeField;
    Ydsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*(d0/cfg.nmic);
    Ymwf(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmwf;
    Ymvdr(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr;
    Ymvdr_superDirective(idx_freq,:) = squeeze(X(idx_freq,:,:))*Wmvdr_superDirective;
    %Ycdr(idx_freq,:)=Ydsb(idx_freq,:).*Wcdr(idx_freq,:);
    Yrfsb(idx_freq,:) = squeeze(X(idx_freq,:,:))*flt.w.RFSB(:,idx_freq);

    %Ynull(idx_freq,:) = squeeze(X(idx_freq,:,:))*wNull;
end
ymvdr = DFTSynRealEntireSignal(Ymvdr, cfg.N, cfg.K, cfg.p);
ymvdr_superDirective = DFTSynRealEntireSignal(Ymvdr_superDirective, cfg.N, cfg.K, cfg.p);
ymwf = DFTSynRealEntireSignal(Ymwf, cfg.N, cfg.K, cfg.p);
ymwf_freeField = DFTSynRealEntireSignal(Ymwf_freeField, cfg.N, cfg.K, cfg.p);
yrfsb = DFTSynRealEntireSignal(Yrfsb, cfg.N, cfg.K, cfg.p);
ydsb = DFTSynRealEntireSignal(Ydsb, cfg.N, cfg.K, cfg.p);
%ycdr = DFTSynRealEntireSignal(Ycdr, cfg.N, cfg.K, cfg.p);

end