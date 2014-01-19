function H = horzcat( F, G )
%HORZCAT Horizontal concatenation of chebfun2v objects.
% 
%  This is not allowed and returns an error.  This function exists so that
%  the error message is meaningful to a Chebfun2 user. 

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information.

error('CHEBFUN2V:HORZCAT','Horizontal concatenation of chebfun2v objects is not supported.')

end