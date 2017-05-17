function out = get(f, propName)
%GET   GET method for the CHEBOP2 class.
%   P = GET(N, PROP) returns the property P specified in the string PROP from
%   the CHEBOP2 N. Valid entries for the string PROP are:
%       'DOMAIN'         - The domain of definition of N.
%       'OP'             - The partial differential operator of N.
%       'LBC'            - The left boundary constraints of N.
%       'RBC'            - The right boundary constraints of N.
%       'UBC'            - The top boundary constraints of N.
%       'DBC'            - The bottom boundary constraints of N.
% 
% The following are also supported: 
%       'DIM'            - The size of the system.
%       'SCALE'          - The relative scale of the solution.
%       'COEFFS'         - The variable coefficients of N.
%       'XORDER'         - The order of differentiation in the x-direction.
%       'YORDER'         - The order of differentiation in the y-direction.
%       'U', 'S', 'V'    - The low rank structure of N.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    out = cell(numel(f));
    for k = 1:numel(f)
        out{k} = get(f(k), propName);
    end
    return
end

switch propName
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(propName);
    otherwise
        error('CHEBFUN:CHEBOP2:get:propName', ...
            [propName,' is not a valid chebop2 property.'])
end
