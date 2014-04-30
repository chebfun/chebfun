function vscale = getvscl( f )
% GETVSCL(F) returns an estimate of the vertical scale of a function. 
% 
% VSCALE = GETVSCL(F) estimates the vertical scale (also known as the dynamical
% range) of a function. This is required because a chebtech does not store its
% interpolation data. If F is an array-valued chebtech with K columns, then 
% VSCALE is a row vector of length K. 

% Check isempty: 
if ( isempty( f.coeffs ) ) 
    vscale = 0;
else
    % compute values
    values = f.coeffs2vals( f.coeffs );
    vscale = max(abs(values), [], 1);
end

end