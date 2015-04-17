function y = feval( f, theta, lambda ) 
% FEVAL    Evaluate a spherefun. 
% 
%  Y = FEVAL( F, THETA, LAMBDA )  evaluates a spherefun F at (THETA,
%  LAMBDA). 

% Shift the evaluation points, annoying: 
hn = pi/length(f.Cols);  % Have to adjust for the shift in y points. Ugly!
lambda = (lambda-hn)/pi-.5; 
theta = theta/pi;

% Check that sizes are correct:
if ( any( size( theta ) ~= size( lambda ) ) )
   error('SPHEREFUN:FEVAL:SIZES', ...
               'Matrix of evaluation points need to be of the same size. ') 
end

% Do the evaluation. 
y = feval(f.Cols, lambda(:)) * f.BlockDiag * feval(f.Rows, theta(:))';
y = diag( y );   % take just the diagonals. 
y = reshape( y, size(lambda,1), size(lambda,2)  ); 

end 