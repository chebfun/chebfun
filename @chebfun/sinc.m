function F = sinc(F, pref)
%SINC   Sinc function of a CHEBFUN.
%   SINC(F) computes the sinc function of the CHEBFUN F, i.e., sin(pi*F)/(pi*F).
%
%   SINC(F, PREF) does the same but uses the CHEBPREF object PREF when computing
%   the composition.
%
% See also SIN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @(x) sin(pi*x)./(pi*x), pref);
end

end
