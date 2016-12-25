function c = getCoeffs(source)
%GETCOEFFS    Get coefficients. Static, private method. 
%   C = GETCOEFFS(SOURCE) returns the Fourier/Chebyshev coefficients.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(source, 'chebmatrix') ) % Note: LINOP is a CHEBMATRIX.
    % Get the coefficients of each block of the CHEBMATRIX:
    c = cell(size(source));
    for k = 1:numel(c)
        try
            c{k} = toCoeff(source.blocks{k});
        catch
            c{k} = [];
        end
    end
else
    % If we didn't get a CHEBMATRIX passed in, still try to see whether we can
    % get coefficients.
    try
        c = {toCoeff(source)};
    catch
        c = {};
    end
end

end
