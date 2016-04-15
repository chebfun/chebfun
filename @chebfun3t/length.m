function varargout = length( F )
%[out1, out2, out3] = length( F )
%LENGTH   Length of a Chebfun3t, i.e., the size of tensor of Chebyshev 
% coefficients. If just one output has been asked, then it is
% the max number of coefficients in three dimensions.

if ( isempty( F ) ) 
    % Output:
    if ( nargout <= 1 )
        varargout = { [] };
    else
        varargout = {[], [], []};  
    end
    return
end

coeffs = F.coeffs;

if numel(size(coeffs))<3                % If the function is univariate or
    [m,n] = size(coeffs); p = 1;        % bivariate, artificially put 1 as
                                        % the size in the 3rd dimension.
else
    [m,n,p] = size(coeffs);
end
    
% Output:
if ( nargout <= 1 )
    varargout = { max([m, n, p]) };
else
    varargout = {m, n, p};  
end

end