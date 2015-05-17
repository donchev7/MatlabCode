function Thi_Y = Thi_Y(cpsd,psd,nmic)

Thi_Y = 0.5*(cov(cpsd)+cov(psd));
% Thi_Y = zeros(nmic,nmic);
% counter=0;
% for m=1:nmic
%     for n=1:nmic
%         if m==n
%             Thi_Y(m,n) = var(psd(:,m),1);
%         else
%             if(n>m)
%                 counter=counter+1;
%                 Thi_Y(m,n) = var(cpsd(:,counter),1);
%             else
%                 Thi_Y(m,n) = Thi_Y(n,m);
%             end
%         end
%     end
% end
end