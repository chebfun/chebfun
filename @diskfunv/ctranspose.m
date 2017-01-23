function F = ctranspose( F )
%CTRANSPOSE   Conjugate transpose of a DISKFUNV.
%    F = CTRANSPOSE(F) is the same as F'. The orientation of the DISKFUNV 
%    is transposed with complex conjugation. 
% 
%    Since all DISKFUNV objects are real-valued, this command is the same
%    as TRANSPOSE(F)
%
% See also DISKFUNV/CTRANSPOSE 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = transpose( F ); 
F = conj( F ); 

end
