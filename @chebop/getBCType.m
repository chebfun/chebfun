function bcType = getBCType(N)

if ( isempty(N.bc) )
    if ( isempty(N.lbc) && isempty(N.rbc) )
        bcType = 'bvp';
    elseif ( isempty(N.rbc) )
        bcType = 'ivp';
    elseif ( isempty(N.lbc) )
        bcType = 'fvp'; 
    else
        bcType = 'bvp';
    end
elseif ( strcmpi(N.bc, 'periodic') )
    bcType = 'periodic';
else
    bcType = 'general';
end

end
