function F = nchoosek(F, G, varargin)
%NCHOOSEK   Choose function for CHEBFUN objects.
%   NCHOOSEK(F, G) computes the continuous analogue of the CHOOSE function
%   for CHEBFUN objects F and G. The result is defined pointwise as
%
%       NCHOOSEK(F,G) = GAMMA(F+1) ./ (GAMMA(G+1) .* GAMMA(F-G+1))
%
%   NCHOOSEK(F, G, PREF) does the same but uses the CHEBFUNPREF object PREF
%   when computing the composition.
%
% See also GAMMA.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method on the continuous analogue of the choose function:
F = compose(F, @(f,g) gamma(f+1)./(gamma(g+1).*gamma(f-g+1)), G, varargin{:});

end
