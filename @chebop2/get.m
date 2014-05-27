function val = get(f,propName)
% GET   Get chebop2 properties.
%

% Copyright 2012 by The University of Oxford and The Chebfun2 Developers.

val = [];
if numel(f) > 1
    val = cell(numel(f));
    for k = 1:numel(f)
        val{k} = get(f(k), propName);
    end
    return
end

switch propName
    case 'lbc'
        val = f.lbc;
    case 'rbc'
        val = f.rbc;
    case 'ubc'
        val = f.ubc;
    case 'dbc'
        val = f.dbc;
    case 'domain'
        val = f.domain;
    case 'op'
        val = f.op;
    case 'coeffs'
        val = f.coeffs;
    case 'coeffcell'
        val = f.coeffcell;
    case 'U'
        val = f.U; 
    case 'S' 
        val = f.S; 
    case 'V'
        val = f.V;
    otherwise
        error('CHEBFUN2:get:propnam',[propName,' is not a valid chebfun2 property.'])
end