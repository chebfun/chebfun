function varargout = length(F)
%LENGTH   Length of a CHEBFUN3
%   The number of Chebyshev coefficients in each of the directions. If 
%   just one output has been asked, then it is the max number of 
%   coefficients in three dimensions.

if ( isempty(F) ) 
    % Output:
    if ( nargout <= 1 )
        varargout = {[]};
    else
        varargout = {[], [], []};  
    end
    return
end

if ( iszero(F) ) 
    if ( nargout <= 1 )
        varargout = {1};
    else 
        varargout = {1, 1, 1};
    end
    return
end

m = length(F.cols);
n = length(F.rows);
p = length(F.tubes);
% Output:
if ( nargout <= 1 )
    varargout = {max([m, n, p])};
else
    % ST decomposition
    varargout = {m, n, p};
end

end