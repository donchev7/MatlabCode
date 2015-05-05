function k = calc_parcor(x,P)
% k = calc_parcor(x,P)
% Calculate PARCOR coefficients
% (see diploma thesis, Appendix B.2)
% k vector with PARCOR coefficients for block x
% x
% P block of input signal vector x
% model order
% Calc. autocorrelation functions
N = length(x);
for i = 1:(P+1)
    sum_x = 0;
    for n = i:(N)
        sum_x = sum_x + x(n)*x(n-i+1);
    end
    r(i) = sum_x;
end
r = r(:);
% Initial conditions
e(1) = r(1);
a_new = 0;
% compute PARCOR coefficients and LP-parameters
for m = 1:P
    sum_a = 0;
    a = a_new;
    for i = 1:(m-1)
        sum_a = sum_a +a (i)*r(m+1-i);
    end
    k(m) = -(sum_a + r(m+1)) / e(m); % PARCOR coefficients
    for i = 1:(m-1)
        a_new(i) = a(i) + k(m)*a(m-i);
    end
    a_new(m) = k(m);
    e(m+1) = (1 - k(m)^2) * e(m);
end
k = k(:);