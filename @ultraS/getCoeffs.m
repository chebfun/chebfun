function c = getCoeffs(source)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
if ( isa(source, 'chebmatrix') )
    c = cell(size(source));
    for k = 1:numel(c)
        try
            c{k} = blockCoeff(source.blocks{k}).coeffs;
        catch
            c{k} = [];
        end
    end
else
    try
        c = {blockCoeff(source).coeffs};
    catch
        c = {};
    end
end
end
