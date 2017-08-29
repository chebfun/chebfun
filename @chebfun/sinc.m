function F = sinc(F, pref)
%SINC   Sinc function of a CHEBFUN.
%   SINC(F) computes the sinc function of the CHEBFUN F, i.e., 
%       sinc(F) := sin(F)/(F).
%
%   SINC(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
%   Note that this definition of the SINC function differs from the MATLAB
%   implementation in the Signal Processing toolbox, which uses
%       sinc(F) := sin(pi*F)/(pi*F).
%
% See also SIN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfunpref();
end

% Call the compose method:
F = compose(F, @mysinc, pref);

end

function out = mysinc(x)
% Deal with the removable singularity at 0 explicitly.
out = sin(x)./x;
out(x == 0) = 1;
end
