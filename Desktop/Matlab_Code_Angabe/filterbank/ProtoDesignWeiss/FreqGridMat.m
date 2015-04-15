function T = FreqGridMat(L,Nw,w_s,W);
% T = FreqGridMat(L,Nw,w_s,W);
%
% Builts a matrix to calculate the real part of the frequency response
% estimated on a frequency grid using Nw equally spaced frequencies 
% between w_s and pi. A constraint for symmetry of the analyzed filter
% is included.
% Optional W imposes a weighting of the frequencies. 
%
% Input parameters:
%      L       length of filter to be analyzed
%      Nw      column dimension of T (number of frequency points)
%      w_s     stopband edge
%      W       (optional) weighting, dimension Nw.  
%
% Output parameters:
%      T       matrix for frequency response in stopband
%
% (c) St.Weiss, University of Strathclyde, 3.9.1997 

% create grid of discrete time and frequency grid
n = (0:1:L-1);
step = (pi-w_s)/Nw;
m = (w_s:step:pi-step);           % Nw equidistant frequency points

% built DCT matrix and impose symmetry constraint
T = cos(m'*n);
T = T(:,1:end/2) + T(:,end:-1:end/2+1);

% optional weighting of the frequency points
if exist('W'),
  if length(W) ~= Nw,
    disp(error(['weighting vector must conform with the number of' ...
	'frequencies analyzed']));
  end;
  T = diag(W)*T;
end;





