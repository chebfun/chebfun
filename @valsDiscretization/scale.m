function scl = scale( disc, uFun ) 
% SCALE   Estimate the vertical scale from a vector of values.
% 
%   SCL = SCALE( DISC, UFUN ) returns an estimate for the vertical scale 
%   of UFUN. VALSDISCRETIZATION expects UFUN to be a cell array of sample data. 

scl = max( cellfun(@(v) max(abs(v)), uFun) );

end
