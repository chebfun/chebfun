function bcType = getBCType(N)

if ( isempty(N.bc) )
    bcType = 'bvp';
elseif ( strcmpi(N.bc, 'periodic') && isempty(N.lbc) && isempty(N.rbc) )
    bcType = 'periodic';
else
    bcType = 'general';
end

end
