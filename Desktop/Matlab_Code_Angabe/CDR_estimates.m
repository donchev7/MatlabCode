function [W_prop2, W_NoDOA] = CDR_estimates(X,d,Cnn,cfg,sig)
CDR_prop2 = zeros(size(X,1),size(X,2));
CDR_NoDOA = zeros(size(X,1),size(X,2));
CDR_prop2_Bf = zeros(size(X,1),size(X,2));
CDR_NoDOA_Bf = zeros(size(X,1),size(X,2));
Pxx = estimate_psd(X,cfg.alpha);
counter=0;
for m = 1:(cfg.nmic-1)
    for n = (m+1):cfg.nmic
         counter = counter + 1;
         TDOA = abs(finddelay(sig.x(:,m),sig.x(:,n)))/cfg.fs;
         Cxx = estimate_cpsd(X(:,:,m),X(:,:,n),cfg.alpha)./sqrt(Pxx(:,:,m).*Pxx(:,:,n));
         CDR_prop2 = CDR_prop2 + CDR_estimate_prop2(Cxx,Cnn(:,m,n),TDOA,cfg);
         CDR_NoDOA = CDR_NoDOA + CDR_estimate_NoDOA(Cxx,Cnn(:,m,n));
    end
end
CDR_prop2 = CDR_prop2/counter;
CDR_NoDOA = CDR_NoDOA/counter;

for idx_freq=1:length(cfg.frange)   
        CDR_prop2_Bf(idx_freq,:) = CDR_prop2(idx_freq,:) / ((d(:,idx_freq)'/(squeeze(Cnn(idx_freq,:,:))))*d(:,idx_freq));
        CDR_NoDOA_Bf(idx_freq,:) = CDR_NoDOA(idx_freq,:) / ((d(:,idx_freq)'/(squeeze(Cnn(idx_freq,:,:))))*d(:,idx_freq));
end

CDR_prop2_Bf = real(CDR_prop2_Bf);
CDR_NoDOA_Bf = real(CDR_NoDOA_Bf);
W_prop2 = CDR_prop2_Bf./(1+CDR_prop2_Bf);
W_NoDOA = CDR_NoDOA_Bf./(1+CDR_NoDOA_Bf);
W_prop2 = max(W_prop2,0.1);
W_NoDOA = max(W_NoDOA,0.1);
W_prop2 = min(W_prop2,1);
W_NoDOA = min(W_NoDOA,1);