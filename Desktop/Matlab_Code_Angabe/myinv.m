function X = myinv(A)
% myinv: returns inverse of A
% Started with bslashtx and edited to return entire inverse

%
 [freq time N] = size(A);
 I = eye(N,N);
 delta = 10e-4;
 B = delta*trace(A(:,:,1))/N;
 C = B*I;
 D = [C+A(:,:,1)];
 A=inv([(delta*trace(A(:,:,1))/N)*I+A(:,:,1)]);
 
 