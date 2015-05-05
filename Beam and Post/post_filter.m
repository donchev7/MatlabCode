function [H_post,alpha,varargout] = post_filter(filter_nr,cn,aa,alpha,X,varargin)
% [H_post, alpha, new return vectors] = ...
% post_filter(filter_nr,cn,aa,alpha,X, input vectors, old return vectors)
% This function calculates the bins of the Post-Filter depending on
% which Filter-Realization is chosen.
% H_post
% alpha
% new return vectors
% Post-Filter bins
% new factor for Welch-estimation; if adaptive
% Welch-factor is used
% new calculated auto- and cross spectral density
% vectors or matrix for the current frame
% filter_nr
% cn
% aa
% alpha
% X
% input vectors
% old return vectors
% name of the chosen filter
% 'ZEL88' ...... Zelinski Filter based on Welch-estimated spectral
% density functions
% 'ZEL88p'...... Zelinski Filter based on Welch-estimation plus
% post-processing method (see Zelinski 1988)
% 'SIM92' ...... Simmer 1992
% 'APAB' ....... Adaptive Postfilter with Arbitrary Beamformer
% (see book "Micorphone Arrays" by Brandstein)
% 'APES' ....... Adaptive Postfilter Extension for
% Superdirective beamformers (see Bitzer et al. 1999)
% 'MCCC' ....... McCowan Filter (see McCowan 2003)
% 'GMCC' ....... McCowan Filter, correct solution for all angles (see
% diploma thesis)
% introduce Comfort Noise to reduce speech distortions
% (Minimum-Filter, see thesis section 4.6.1)
% 0 ..... no Comfort Noise (default)
% 1 ..... activate Comfort Noise
% use adaptive Welch-parameter
% (see thesis, section 4.6.2)
% 0 ..... use fixed Welch-parameter (default)
% 1 ..... use adaptive Welch-parameter
% previous factor for Welch-estimation
% time-aligned input signal bins
% other input vectors, depends on Post-Filter
% auto- and cross spectral density vectors or matrix
% of the previous frame
% hint:
% Please declare and initialize the return vectors in the master program

if nargin<2
help post_filter
return;
end
[N2,K] = size(X);
comb = 2/(K*(K-1));
% all possible sensor combinations
% Calc. minimum Postfilter function; to avoid zeros in the transfer
% function
f1 = ceil(512/16000*200);
f2 = ceil(512/16000*4500);
f = [f1:1:f2].';
H_min_func = [zeros(f1-1,1);0.7/(f2-f1).*(f-f1); 0.7.*ones(N2-f2,1)]; % see diploma thesis Boigner
% *************************************************************************
% Parameter, Minimum Postfilter function if choosen
% introduces Comfort Noise, but avoids speech degradation
if cn == 0
    H_min = 0.05;
else
    H_min = H_min_func;
end
%*************************************************************************
% Calculate Postfilters
switch filter_nr
    case 'SIM92'
    %**************************************************************
    % Estimation of postfiler transfer function =>
    % Simmer with Welch-estimation
    Y = varargin{1};
    % signal at beamformer output
    cross_dens = varargin{2};
    % Cross-PSD of time-aligned input signals
    Pyy= varargin{3};
    % Auto-PSD of beamformer output
    % Calc. Auto-PSD of beamformer output using Welch-formula
    Pyy = welch_est(Pyy,Y,Y,alpha);
    % Calc. Cross-PSD of time-aligned input signals using Welch-formula
    counter = 0;
    for m = 1:(K-1)
        for n = (m+1):K
            counter = counter + 1;
            cross_dens(:,counter) = welch_est(cross_dens(:,counter),X(:,m),X(:,n),alpha);
        end
    end
    % Calc. averaged Cross-PSD over all sensor comb.
    C = sum(real(cross_dens.')).';
    C = comb.*C;
    % Calc. Postfilter
    H_post = C./Pyy;
    if nargout == 4
        varargout{1} = cross_dens;
        varargout{2} = Pyy;
    else
        error('Wrong number of output arguments!!')
    end
case 'ZEL88'
    %**********************************************************************
    % Estimation of postfilter transfer function =>
    % Zelinski with Welch-estimation
    auto_dens = varargin{1};
    % Auto-PSD of time-aligned input signals
    cross_dens = varargin{2};
    % Cross-PSD of time-aligned input signals
    % Calc. Auto-PSD of time-aligned input signals using Welch-formula
    for m = 1:K
    auto_dens(:,m) = welch_est(auto_dens(:,m),X(:,m),X(:,m),alpha);
    end
    % Calc. Cross-PSD of time-aligned input signals using Welch-formula
    counter = 0;
    for m = 1:(K-1)
    for n = (m+1):K
    counter = counter + 1;
    cross_dens(:,counter) = welch_est(cross_dens(:,counter),X(:,m),X(:,n),alpha);
    end
    end
    % Calc. averaged Cross-PSD over all sensor comb.
    C = sum(real(cross_dens.')).';
    C = comb.*C;
    % Calc. averaged Auto-PSD over all sensor comb.
    A = sum(auto_dens.').';
    A = A./K;
    % Calc. Postfilter
    H_post = C./A;
    if nargout == 4
    varargout{1} = auto_dens;
    varargout{2} = cross_dens;
    else
    error('Wrong number of output arguments!!')
    end
case 'ZEL88p'
    %**********************************************************************
    % Estimation of postfilter transfer function =>
    % Zelinski with Welch-estimation incl. post-processing method
    auto_dens = varargin{1};
    % Auto-PSD of time-aligned input signals
    cross_dens = varargin{2};
    % Cross-PSD of time-aligned input signals
    % Calc. Auto-PSD of time-aligned input signals using Welch-formula
    for m = 1:K
    auto_dens(:,m) = welch_est(auto_dens(:,m),X(:,m),X(:,m),alpha);
    end
    % Calc. Cross-PSD of time-aligned input signals using Welch-formula
    counter = 0;
    for m = 1:(K-1)
    for n = (m+1):K
        counter = counter + 1;
        cross_dens(:,counter) = welch_est(cross_dens(:,counter),X(:,m),X(:,n),alpha);
    end
    end
    % Calc. averaged Cross-PSD over all sensor comb.
    C_tilde = sum(real(cross_dens.')).';
    C_tilde = comb.*C_tilde;
    % post-processing method (see Zelinski 1988)
    S_2 = max(C_tilde,0);
    S_2 = S_2.^2;
    S_2 = smooth(S_2);
    V = real(cross_dens).';
    V = min(V,0);
    for m = 1:N2
        M(m) = nnz(V(:,m));
    end
    M = max(M,1);
    V = sum(V.^2);
    V = (V./M).';
    V = smooth(V);
    alpha_k = S_2./(S_2 + V.*comb);
    C = alpha_k.*C_tilde;
    % Calc. averaged Auto-PSD over all sensor comb.
    A = sum(auto_dens.').';
    A = A./K;
    % Calc. Postfilter
    H_post = C./A;
    if nargout == 4
        varargout{1} = auto_dens;
        varargout{2} = cross_dens;
    else
        error('Wrong number of output arguments!!')
    end
case 'APAB'
    %**********************************************************************
    % Estimation of postfilter transfer function =>
    % APAB Postfilter
    Y = varargin{1};
    % signal at beamformer output
    Pxx = varargin{2};
    % Auto-PSD of time-aligned and averaged input signals
    Pyy = varargin{3};
    % Auto-PSD of beamformer output
    % average time-aligned input signals
    X = sum(X.').';
    X = X./K;
    % Calc. Auto-PSD of time-aligned and averaged input signals using
    % Welch-formula
    Pxx = welch_est(Pxx,X,X,alpha);
    % Calc. Auto-PSD of beamformer output using Welch-formula
    Pyy = welch_est(Pyy,Y,Y,alpha);
    % Calc. Postfilter
    H_post = Pyy./Pxx;
    if nargout == 4
        varargout{1} = Pxx;
        varargout{2} = Pyy;
    else
        error('Wrong number of output arguments!!')
    end
case 'APES'
    %**********************************************************************
    % Estimation of postfilter transfer function =>
    % APES Postfilter
    Y = varargin{1};
    % signal at DSB output
    Z = varargin{2};
    % signal at SDB output
    Pxx = varargin{3};
    % Auto-PSD of time-aligned and averaged input signals
    Pyy = varargin{4};
    Pzz = varargin{5};
    % Auto-PSD of DSB output
    % Auto-PSD of SDB output
    % Calc. Auto-PSD of DSB output using Welch-formula
    Pyy = welch_est(Pyy,Y,Y,alpha);
    % Calc. Auto-PSD of SDB output using Welch-formula
    Pzz = welch_est(Pzz,Z,Z,alpha);
    % Calc. Auto-PSD of time-aligned input signals using Welch-formula
    for m = 1:K
        Pxx(:,m) = welch_est(Pxx(:,m),X(:,m),X(:,m),alpha);
    end
    % Calc. nominator of APES Algorithm
    dum = sum(Pxx.').';
    dum = (1/K^2).*dum;
    nom = Pyy-dum;
    % average nominator over all sensor comb.
    nom = (K/(K-1)).*nom;
    % Calc. First postfilter
    H1 = nom./Pyy;
    % Calc. Second postfilter
    H2 = Pzz./Pyy;
    % Calc. Postfilter
    H_post = H1.*H2;
    if nargout == 5
        varargout{1} = Pxx;
        varargout{2} = Pyy;
        varargout{3} = Pzz;
    else
        error('Wrong number of output arguments!!')
    end
case 'MCC03'
    %**************************************************************
    % Estimation of the postfilter transfer function =>
    % McCowan 2003
    Gamma_comb = varargin{1};
    % Coherence matrix over all sensor comb.
    auto_dens = varargin{2};
    % Auto-PSD of time-aligned input signals
    cross_dens = varargin{3};
    % Cross-PSD of time-aligned input signals
    % Calc. Auto-PSD of time-aligned input signals using Welch-formula
    for m = 1:K
        auto_dens(:,m) = welch_est(auto_dens(:,m),X(:,m),X(:,m),alpha);
    end
    % Calc.
    counter =0;
    %for estimation of PSD of the speech signal (see McCowan 2003)
    for m=1:(K-1)
        for n = (m+1):K
            counter = counter + 1;
            % Calc. Cross-PSD of time-aligned input signals using
            % Welch-formula
            cross_dens(:,counter) = welch_est(cross_dens(:,counter),X(:,m),X(:,n),alpha);
            delta1(:,counter) = real(cross_dens(:,counter));
            % realpart of PSD
            delta2(:,counter) = auto_dens(:,m) + auto_dens(:,n); % average of Auto-PSDs m and n
        end
    end
    dum = 0.5.*Gamma_comb;
    delta2 = dum.*delta2;
    dum = delta1 - delta2;
    nom = 1 - Gamma_comb;
    % calc. denominator of speech PSD
    C_ss = dum./nom;
    % estimated PSD of speech signal
    % Calc. averaged speech-PSD over all sensor comb.
    C = sum(C_ss.').';
    C = comb.*C;
    % Calc. averaged Auto-PSD over all sensor comb.
    A = sum(auto_dens.').';
    A = A./K;
    % Calc. Postfilter
    H_post = C./A;
    if nargout == 4
        varargout{1} = auto_dens;
        varargout{2} = cross_dens;
    else
        error('Wrong number of output arguments!!')
    end
case 'GMCC'
    %**************************************************************
    % Estimation of the postfilter transfer function =>
    % McCowan03, correct solution,calc. in my diploma thesis
    Gamma_comb = varargin{1};
    % Coherence matrix over all sensor comb.
    d0 = varargin{2}.';
    % steering vector
    auto_dens = varargin{3};
    % Auto-PSD of time-aligned input signals
    cross_dens = varargin{4};
    % Cross-PSD of time-aligned input signals
    % Calc. Auto-PSD of time-aligned input signals using Welch-formula
    for m = 1:K
        auto_dens(:,m) = welch_est(auto_dens(:,m),X(:,m),X(:,m),alpha);
    end
    % Calc. estimation of PSD of the speech signal (see McCowan 2003)
    counter = 0;
    for m = 1:(K-1)
        for n = (m+1):K
            counter = counter + 1;
            % Calc. new Gama_comb with steering vector
            Gamma_comb(:,counter) = real(conj(d0(:,m)).*d0(:,n)).*Gamma_comb(:,counter);
            % Calc. Cross-PSD of time-aligned input signals using
            % Welch-formula
            cross_dens(:,counter) = welch_est(cross_dens(:,counter),X(:,m),X(:,n),alpha);
            delta1(:,counter) = real(cross_dens(:,counter));
            % realpart of PSD
            delta2(:,counter) = auto_dens(:,m) + auto_dens(:,n); % average of Auto-PSDs m and n
        end
    end
    dum = 0.5.*Gamma_comb;
    delta2 = dum.*delta2;
    dum = delta1 - delta2;
    nom = 1 - Gamma_comb;
    % calc. denominator of speech PSD
    C_ss = dum./nom;
    % estimated PSD of speech signal
    % Calc. averaged speech-PSD over all sensor comb.
    C = sum(C_ss.').';
    C = comb.*C;
    % Calc. averaged Auto-PSD over all sensor comb.
    A = sum(auto_dens.').';
    A = A./K;
    % Calc. Postfilter
    H_post = C./A;
    if nargout == 4
        varargout{1} = auto_dens;
        varargout{2} = cross_dens;
    else
        error('Wrong number of output arguments!!')
    end
    otherwise
        error('Wrong Filter-Name!!!')
end
% restrict Postfilter
H_post = max(H_post,H_min);
H_post = min(H_post,1);
% Calc. current frequency-dependent alpha
if aa == 1
    alpha_fact = 0.3;
    alpha = 0.98 - alpha_fact.*H_post;
end