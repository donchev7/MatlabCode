function E = PolyPhaseMat(h_p,K,option);
% E = PolyPhaseMat(h_p,K,option);
%
% PolyPhaseMat() builts a matrix to check the perfect reconstruction
% property for K-channel filterbanks, which are derived from a single 
% prototype filter h_p.
%
% Input parameters:
%      h_p      prototype filter
%      K        number of channels / polyphases
%      option   'uncon'      unconstraint
%               'uncon_corr' unconstraint corrected
%               'con_corr'   constraint corrected (default)
%
% Output parameters:
%      E        Convolutional matrix for polyphase components
%
% References:
%    [1] PP Vaidyanathan: Multirate Systems and Filter Banks. Prentice
%        Hall, 1993, pages 159-161.
%
% (c) St.Weiss, University of Strathclyde, 3.9.1997

h_p = h_p(:)';
Lp = size(h_p,2);
S = mod(Lp,K);
if S ~=0,
  disp('Warning: prototype filter length is not a multiple of the channel number K');
  h_p = [h_p  zeros(1,S)];
  Lp = length(h_p);
end;  
Ph = zeros(K,Lp/K);
Ph(:) = h_p; 

if ~exist('option');
  option = 'con_corr';
end;  
  
% polyphases of prototype filter

if strcmp(option,'uncon') == 1,
  % unconstraint convolutional matrix 
  E = zeros(2*Lp/K-1,Lp);
  for m = 1:Lp/K,
    E(m:m+Lp/K-1,m:Lp/K:Lp) = Ph'; 
  end;
elseif strcmp(option,'uncon_corr') == 1,
  % unconstraint convolutional matrix (corrected for multiplication with h_p)
  E = zeros(2*Lp/K-1,Lp);
  for m = 1:Lp/K,
    E(m:m+Lp/K-1,(m-1)*K+1:m*K) = fliplr(Ph');
  end;
elseif strcmp(option,'con_corr') == 1,
  % constraint convolutional matrix (corrected for multiplication with h_p)
  % the correction appears as an interleaving of the polyphase filters  
  E = zeros(2*Lp/K-1,Lp);
  for m = 1:Lp/K,
    E(m:m+Lp/K-1,m*K:-1:(m-1)*K+1) = Ph';
  end;
  % the following imposed a symmetry constraint on the filter  
  E = E(:,1:end/2) + E(:,end:-1:end/2+1); 
else 
  disp(error('requested option has not been implemented'));
end;  