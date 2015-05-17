function xy = covarianceEstimate(x)
[m,n] = size(x);
%xc = bsxfun(@minus,x,sum(x,1)/m);  % Remove mean
xy = (x' * x) / (m-1);
%xy = (x' * x);
end
  