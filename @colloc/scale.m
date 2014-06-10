function scl = scale( disc, uFun ) 
% SCALE   Find the vertical scale from a 2nd kind grid.

scl = max( cellfun(@max, uFun) );

end