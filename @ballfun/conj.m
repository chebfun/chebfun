function g = conj(f)
% CONJ Conjugate of a BALLFUN function
%   CONJ(f) is the conjugate of the BALLFUN function f
g = real(f)-1i*imag(f);
end
