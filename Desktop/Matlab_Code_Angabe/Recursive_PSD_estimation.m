function Phi = Recursive_PSD_estimation(X1,X2,alpha)
if nargin == 2
    alpha = X2;
    X2=X1;
end

Phi = zeros(size(X1,1),size(X1,2));
for k=1:size(X1,1) %freq
    for l=1:size(X1,2)-1 %time
        Phi(k,l+1)=alpha*Phi(k,l)+(1-alpha)*X1(k,l+1)*X2(k,l+1)';
    end
end
% if nargin==2
%     Phi = real(Phi);
% end
%Phi = real(Phi);
end