function vals = feval(f, r, lambda, theta)
% FEVAL Evaluate a BALLFUNV
%   FEVAL(f, r, lambda, theta) is the array of values of the BALLFUNV f 
%   at the grid r x lambda x theta
F = f.comp;

% Get the size of the lists
Nr = length(r);
Nlam = length(lambda);
Nth = length(theta);

vals = zeros(Nr, Nlam, Nth, 3);

vals(:,:,:,1) = feval(F{1},r,lambda,theta);
vals(:,:,:,2) = feval(F{2},r,lambda,theta);
vals(:,:,:,3) = feval(F{3},r,lambda,theta);
end
