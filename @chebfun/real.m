function F = real(F)
%REAL   Complex real part of a CHEBFUN.
%   REAL(F) is the real part of F.
%
% See also IMAG, ISREAL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Handle the empty case:
if ( isempty(F) )
    return
end

for j = 1:numel(F)
    % Take real part of the pointValues:
    F(j).pointValues = real(F(j).pointValues);

    % Take real part of the FUNs:
    for k = 1:numel(F(j).funs)
        F(j).funs{k} = real(F(j).funs{k});
    end
end

end
