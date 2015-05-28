function y = feval( f, lambda, theta ) 
% FEVAL    Evaluate a spherefun. 
% 
%  Y = FEVAL( F, THETA, LAMBDA )  evaluates a spherefun F at (THETA,
%  LAMBDA). 

% % Shift the evaluation points, annoying: 
% hn = pi/length(f.Cols);  % Have to adjust for the shift in y points. Ugly!
% lambda = (lambda-hn)/pi-.5; 

% lambda = lambda/pi;
% theta = theta/pi - 0.5;

y = feval@separableApprox( f, lambda, theta );

% % Check that sizes are correct:
% if ( any( size( theta ) ~= size( lambda ) ) )
%    error('SPHEREFUN:FEVAL:SIZES', ...
%                'Matrix of evaluation points need to be of the same size. ') 
% end
% 
% % Get the low rank representation for f.
% [cols, D, rows] = cdr(f);
% 
% % Do the evaluation. 
% % y = feval(f.Cols, theta(:)) * f.BlockDiag * feval(f.Rows, lambda(:))';
% y = feval(f.Cols, theta(:)) * D * feval(f.Rows, lambda(:))';
% y = reshape( y, size(lambda,1), size(lambda,2)  ); 

end 