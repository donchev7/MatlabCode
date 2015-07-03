function Wmwf = postFilter(W,Gamma,d0,Thi_Y,nmic)

V = 1/(nmic-1) * trace((eye(nmic,nmic)-d0*W')*Thi_Y*inv(Gamma));%EQ (7a)
%The variance of the interference has to be corrected by the beamformer
%suppression factor
Vo = V/(d0'*inv(Gamma)*d0);
S = W'*(Thi_Y-V*Gamma)*W;%EQ (7b)

Wmwf = (S/(S+Vo)); %EQ (5b)
end