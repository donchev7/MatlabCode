function cfg = micSpacingCalculation(mics)
[K,Dim] = size(mics);
if (Dim < 2) || (K < 1)
    cfg = mics;
end
if Dim == 2
    rn = [mics zeros(K,1)];
else
    rn = mics;
end
% Define Matrix of the micorphone distances
xc = rn(:,1);
xc = xc(:,ones(K,1));
dxc = xc - xc.';
yc = rn(:,2);
yc = yc(:,ones(K,1));
dyc = yc - yc.';
if Dim == 2
    dR = sqrt(dxc.^2 + dyc.^2);
else
    zc = rn(:,3);
    zc = zc(:,ones(K,1));
    dzc = zc - zc.';
    dR = sqrt(dxc.^2 + dyc.^2 + dzc.^2);
end
cfg = dR;
end