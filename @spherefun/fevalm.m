function y = feval_meshgrid( f, lambda, theta ) 
% FEVAL    Evaluate a spherefun at meshgrid of points. 
% 
%  Y = FEVAL( F, THETA, LAMBDA )  evaluates a spherefun F at 
%  a meshgrid. (More efficient for plotting if this structure 
% is exploited.) 
% 
% Private command. 

% TODO: Make this work for less standard meshgrids. 

% % Extract out row and column:
% lambda = lambda(1,:);
% theta = theta(:,1); 
% 
% % Shift the evaluation points, annoying: 
% hn = pi/length(f.Cols);  % Have to adjust for the shift in y points. Ugly!
% lambda = (lambda-hn)/pi-.5; 
% theta = theta/pi;

% % Do the evaluation. 
% y = feval(f.Cols, lambda(:)) * f.BlockDiag * feval(f.Rows, theta(:))';
% y = transpose( reshape( y, size(theta,1), size(lambda,2)  ) ); 

% lambda = lambda/pi;
% theta = theta/pi - 0.5;

y = fevalm@separableApprox( f, lambda, theta );

end 