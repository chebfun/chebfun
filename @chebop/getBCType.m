function bcType = getBCType(N)

if ( isempty(N.bc) )
    if ( isempty(N.lbc) && isempty(N.rbc) )
        bcType = 'nobc';
    elseif ( isempty(N.rbc) )
%         bcType = 'ivp'; % For now, just treat as BVP case.
        bcType = 'bvp';
    elseif ( isempty(N.lbc) )
%         bcType = 'fvp'; % For now, just treat as BVP case.
        bcType = 'bvp'; 
    else
        bcType = 'bvp';
    end
elseif ( strcmpi(N.bc, 'periodic') )
    bcType = 'periodic';
else
    bcType = 'general';
end

end