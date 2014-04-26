function dimAdjust = getDimAdjust(L)

% TODO: Document.

% This is used to determine the different discretisation sizes for
% different variables in a system of ODES. The approach taken is that
% thehighest appearing derivatives of each variable should be discretized at
% the same number of points.

if ( isa(L, 'linop') )
    
    % The input adjustment size of the (j,k) entry is max(diffOrder(:,k))
    dimAdjust = max(getDiffOrder(L), [], 1);
    dimAdjust = max(dimAdjust, 0);
    dimAdjust = repmat(dimAdjust, size(L, 1), 1);
    
else
    
    % Other type of course are not projected.
    dimAdjust = 0;
    
end

end