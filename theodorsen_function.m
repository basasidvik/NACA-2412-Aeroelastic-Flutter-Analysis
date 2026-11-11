function Ck = theodorsen_function(k)
% THEODORSEN_FUNCTION  Classical Theodorsen circulatory function C(k).
%
%   Ck = theodorsen_function(k)
%
% k is the positive reduced frequency, k = omega*c/(2*U). The zero-frequency
% limit is C(0)=1. The convention used here is
%
%   C(k) = H_1^(2)(2k) / (H_1^(2)(2k) + i H_0^(2)(2k)).

    if ~isscalar(k) || ~isreal(k) || k < 0
        error('k must be a nonnegative real scalar.');
    end

    if k < 1e-8
        Ck = 1.0;
        return;
    end

    H0 = besselh(0, 2, 2*k);
    H1 = besselh(1, 2, 2*k);
    Ck = H1 / (H1 + 1i*H0);
end
