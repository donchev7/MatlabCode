function Wmwf = postFilter2(W,Gamma,d0,Thi_Y,nmic)
V = 1/(nmic-1) * trace((eye(nmic,nmic)-d0*W')*Thi_Y*inv(Gamma));
S = W'*(Thi_Y-V*Gamma)*W;
Wmwf = S/(S+(V*inv(d0'*inv(Gamma)*d0)));
end