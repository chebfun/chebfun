function c = getCoeffs( source )
%GETCOEFFS     Get coefficients, private method. 
% 
% C = GETCOEFFS( SOURCE ) returns the Chebyshev T coefficients 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( isa(source, 'chebmatrix') )
    c = cell(size(source));
    for k = 1:numel(c)
        try
            c{k} = toCoeff(source.blocks{k});
        catch
            c{k} = [];
        end
    end
else
    try
        c = {toCoeff(source)};
    catch
        c = {};
    end
end
end
