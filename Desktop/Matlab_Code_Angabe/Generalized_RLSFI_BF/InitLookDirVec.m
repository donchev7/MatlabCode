%--------------------------------------------------------------------------
%                        Look Direction Vector                                  
%--------------------------------------------------------------------------
function [cfg] = InitLookDirVec(cfg)
    
    switch(cfg.int_choice)
        %------------------------------------------------------------------
        % polynomial interpolation according to equation (4) in [Mabande et al,2010]
        %------------------------------------------------------------------
        case 1
            % create vector d
            %d includes all vectors d_{i}, i=0,...,I-1, with I being the
            % number of protype look directions, and
            % d_{i}=[D_{i}^{0}...D_{i}^{P}]^{T} includes all interpolation 
            % values D (interpolation factors for each filter-ans-sum unit)
            %d = zeros(cfg.nmbr_look_dir,cfg.P+1); 
            p = 0:cfg.P;
            for idx_look_dir=1:cfg.nmbr_look_dir
                cfg.d(idx_look_dir,:)=cfg.D_i(idx_look_dir).^p; % cfg.D_i: corresponding D for each look direction; D_i = (cfg.DesAng - 90)/90 for linear array and
                %cfg.d = cfg.D_i,cfg.nmbr_look_dir=1 cfg.P=0         % D = (cfg.DesAng - pi/cfg.N)/(pi/cfg.N) for circular array
            end

        case 2 %not implemented yet
            warning('Warning: There is no other interpolation implemented but polynomial interpolation at the moment.');
            exit(1);
    end
    %----------------------------------------------------------------------
    % Apply kronecker product (see Eq. (4) in [Mabande et al,2010])
    % Result D of dimension [I, N, N(P+1)] including all D_{i}, i=0,...,I-1
    % of dimension [N, N(P+1)], with N = #Microphones   
    for idx_look_dir = 1:cfg.nmbr_look_dir 
        cfg.D(idx_look_dir,:,:) = kron(eye(cfg.N),cfg.d(idx_look_dir,:));
    end
    %--------------------------------------------------------------------------        
end