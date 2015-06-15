function  [ishappy, epslevel, cutoff] = bvpV5Check(f, varargin)
%BVPV5CHECK  
%   [ISHAPPY, EPSLEVEL, CUTOFF] = BVPV5CHECK(F, OP, VALUES) tests if the
%   This routine mimics the behavior of solvebvpLinear before the addition of
%   standardChop. 
%
% See also CLASSICCHECK, LOOSECHECK, STRICTCHECK, SAMPLETEST.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% call plateauCheck
[ishappy, epslevel, cutoff] = plateauCheck(f, varargin{:});

% create dummy chebfun and call simplify
g = f; simplify(g);

% set cutoff
cutoff = length(g);

end
