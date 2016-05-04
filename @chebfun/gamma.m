function F = gamma(F, varargin)
%GAMMA   Gamma function of a CHEBFUN.
%   GAMMA(F) computes the composition of the gamma function with F. For
%   example, a chebfun representing the gamma function on the interval [0.1,3]
%   may be created with
%       x = chebfun('x', [0.1,3]);
%       f = gamma(x);
%
%   This method does not currently support introducing poles. To create a
%   chebfun with poles representing the gamma function, try for example
%       f = chebfun(@gamma,[-4 4.2],'splitting','on','blowup','on');
%
%   GAMMA(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @gamma, varargin{:});

end
