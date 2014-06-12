function scl = scale( disc, uFun ) 
% SCALE   Estimate the vertical scale from a vector of values.
% 
%   SCL = SCALE( DISC, UFUN ) returns an estimate for the vertical scale 
%   of UFUN. COLLOC expects UFUN to be a cell array of sample data. 

% TODO: Should this be taking the absolute values of uFun? (This was
% originally copied from linop/expm.)
scl = max( cellfun(@max, uFun) );

end