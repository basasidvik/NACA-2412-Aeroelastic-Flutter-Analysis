function Sk = sears_function(k)
% SEARS_FUNCTION  Classical Sears gust-transfer approximation.
%
% The relation used here is
%   S(k) = C(k)[J0(k)-i J1(k)] + i J1(k),
% with C(k) evaluated using Theodorsen's function.

    if ~isscalar(k) || ~isreal(k) || k < 0
        error('k must be a nonnegative real scalar.');
    end
    if k < 1e-8
        Sk = 1.0;
    else
        Sk = theodorsen_function(k) * ...
             (besselj(0,k) - 1i*besselj(1,k)) + 1i*besselj(1,k);
    end
end
