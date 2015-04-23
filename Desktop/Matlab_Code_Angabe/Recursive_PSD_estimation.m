function Phi = Recursive_PSD_estimation(spectrum,alpha)

Phi = zeros(size(spectrum,1),size(spectrum,2),size(spectrum,3));

for idx_mic=1:size(spectrum,3) %mic
    for k=1:size(spectrum,1)-1 %freq
        for l=1:size(spectrum,2)-1 %time
            Phi(k,l+1,idx_mic)=alpha*Phi(k,l,idx_mic)+(1-alpha)*spectrum(k,l+1,idx_mic)*spectrum(k,l+1,idx_mic)';
        end
    end
end
    
end