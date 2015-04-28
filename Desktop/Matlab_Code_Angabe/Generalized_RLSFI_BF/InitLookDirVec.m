%--------------------------------------------------------------------------
%                        Look Direction Vector                                  
%--------------------------------------------------------------------------
function [local] = InitLookDirVec(cfg,local)
    
% create vector d
%d includes all vectors d_{i}, i=0,...,I-1, with I being the
% number of protype look directions, and
% d_{i}=[D_{i}^{0}...D_{i}^{P}]^{T} includes all interpolation 
% values D (interpolation factors for each filter-ans-sum unit)
%d = zeros(cfg.nmbr_look_dir,cfg.P+1); 
p = 0:local.P;
for idx_look_dir=1:cfg.nmbr_look_dir
    local.d(idx_look_dir,:)=local.D_i(idx_look_dir).^p; % cfg.D_i: corresponding D for each look direction; D_i = (cfg.DesAng - 90)/90 for linear array and
    %cfg.d = cfg.D_i,cfg.nmbr_look_dir=1 cfg.P=0         % D = (cfg.DesAng - pi/cfg.N)/(pi/cfg.N) for circular array
end
%----------------------------------------------------------------------
% Apply kronecker product (see Eq. (4) in [Mabande et al,2010])
% Result D of dimension [I, N, N(P+1)] including all D_{i}, i=0,...,I-1
% of dimension [N, N(P+1)], with N = #Microphones   
for idx_look_dir = 1:cfg.nmbr_look_dir 
    local.D(idx_look_dir,:,:) = kron(eye(cfg.nmic),local.d(idx_look_dir,:));
end
%--------------------------------------------------------------------------        
end