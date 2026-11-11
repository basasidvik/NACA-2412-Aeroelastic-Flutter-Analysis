function Ck = theodorsen_rational_function(k)
% THEODORSEN_RATIONAL_FUNCTION  Low-order Jones rational approximation.
%
%   C(k) ~= 1 - 0.165/(1 - 0.0455*i/k)
%             - 0.335/(1 - 0.30*i/k).

    if ~isscalar(k) || ~isreal(k) || k < 0
        error('k must be a nonnegative real scalar.');
    end
    if k < 1e-8
        Ck = 1.0;
    else
        Ck = 1 - 0.165/(1 - 0.0455i/k) ...
               - 0.335/(1 - 0.30i/k);
    end
end
