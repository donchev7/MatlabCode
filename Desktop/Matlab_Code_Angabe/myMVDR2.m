function W = myMVDR2(Gamma,d0,beamFormer,nmic)

switch beamFormer
    case 'DSB'
        W = d0/nmic;
    case 'MVDR'
        W = inv(Gamma)*d0/(d0'*inv(Gamma)*d0);
end
% mue = 10^(-20/10); %regularization constant for diagonal loading;
% Gamma_const = tril(squeeze(Cnn(idx_nu,:,:)),-1)./(1+mue) + diag(diag(squeeze(Cnn(idx_nu,:,:)))) + ...
%     triu(squeeze(Cnn(idx_nu,:,:)),1)./(1+mue);
% W(:,idx_nu) = (inv(Gamma_const)*v_k) / ...
%     (v_k'*inv(Gamma_const)*v_k);
% W(:,idx_nu) = (inv(Pnn(:,:,idx_nu))*v_k) / ...
% (v_k'*inv(Pnn(:,:,idx_nu))*v_k);
% 
% V(:,idx_nu) = 1/(cfg.nmic-1) * trace((eye(cfg.nmic,cfg.nmic)-v_k*W(:,idx_nu)')*Thi_Y(:,:,idx_nu)*inv(Gamma_const));
% 
% S(:,idx_nu) = W(:,idx_nu)'*(Thi_Y(:,:,idx_nu)-V(:,idx_nu)*Gamma_const)*W(:,idx_nu);
% 
% Wmwf(:,idx_nu) = (S(:,idx_nu)/(S(:,idx_nu)+(V(:,idx_nu)*inv(v_k'*inv(Gamma_const)*v_k))))*...
%     W(:,idx_nu);


end