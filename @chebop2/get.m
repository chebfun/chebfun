function val = get(f,propName)
% GET   GET method for the CHEBOP2 class
%
%   P = GET(N, PROP) returns the property P specified in the string PROP from
%   the CHEBOP2 N. Valid entries for the string PROP are:
%       'DOMAIN'         - The domain of defintion of N.
%       'OP'             - The partial differential operator of N.
%       'LBC'            - The left boundary constraints of N.
%       'RBC'            - The right boundary constraints of N.
%       'UBC'            - The top boundary constraints of N.
%       'DBC'            - The bottom boundary constraints of N.
% 
% The following are also supported: 
%       'COEFFS'         - The variable coefficients of N.
%       'U', 'S', 'V'    - The low rank structure of N.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
    case 'U'
        val = f.U; 
    case 'S' 
        val = f.S; 
    case 'V'
        val = f.V;
    otherwise
        error('CHEBFUN2:get:propnam',[propName,' is not a valid chebfun2 property.'])
end