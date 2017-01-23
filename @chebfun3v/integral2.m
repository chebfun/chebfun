function v = integral2(F, S)
%INTEGRAL2   Flux integral of a Chebfun3v object through a 2D-surface.
%   INTEGRAL2(F, S) computes the flux integral of the Chebfun3v object F 
%   through the parametric surface S defined as a Chebfun2v object (with 3 
%   components):
%                           /
%       INTEGRAL2(F, S) =  |   < F, dS >
%                          /S
%                  //
%               = ||   < F(S(x,y)), cross(S_x(x,y), S_y(x,y)) > dxdy.
%                 //D
% 
% See also CHEBFUN3V/INTEGRAL, CHEBFUN3/INTEGRAL, CHEBFUN3/INTEGRAL2,
% CHEBFUN3/INTEGRAL3, CHEBFUN3/SUM, CHEBFUN3/SUM2 and CHEBFUN3/SUM3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check that F is a chebfun3v object with 3 components:
if ( F.nComponents ~= 3 )
    error('CHEBFUN:CHEBFUN3V:integral2:dimChebfun3v', ...
        'The chebfun3v object must have 3 components.')
end

% Check that S is a chebfun2v object with 3 components:
if ( isa(S, 'chebfun2v') )
    if ( S.nComponents ~= 3 )
        error('CHEBFUN:CHEBFUN3V:integral2:dimChebfun2v', ...
            'The parametrisation must have 3 components.')
    end
else
    error('CHEBFUN:CHEBFUN3V:integral2:notaChebfun2v', ...
        'The parametrisation must be a chebfun2v object.')
end

% Get components of F:
F1 = F.components{1};
F2 = F.components{2};
F3 = F.components{3};

% Get components of the surface S:
S1 = S(1);
S2 = S(2);
S3 = S(3);

% Build the composition F(S):
FS1_handle = @(x,y) feval(F1, feval(S1, x, y), feval(S2, x, y), ...
    feval(S3, x, y));

FS2_handle = @(x,y) feval(F2, feval(S1, x, y), feval(S2, x, y), ...
    feval(S3, x, y));

FS3_handle = @(x,y) feval(F3, feval(S1, x, y), feval(S2, x, y), ...
    feval(S3, x, y));

FS = chebfun2v(FS1_handle, FS2_handle, FS3_handle, S1.domain);

% Normal vector to the surface:
dS = normal(S);

% By definition:
v = sum2(dot(FS, dS));

end