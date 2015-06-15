function F = imag(F)
%IMAG   Complex imaginary part of a CHEBFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Handle the empty case:
if ( isempty(F) )
    return
end

for j = 1:numel(F)
    % Take imaginary part of the pointValues:
    F(j).pointValues = imag(F(j).pointValues);

    % Take imaginary part of the FUNs:
    for k = 1:numel(F(j).funs)
        F(j).funs{k} = imag(F(j).funs{k});
    end
end

end
