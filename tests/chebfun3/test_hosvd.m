function pass = test_hosvd(pref)
% Test HOSVD.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 10*pref.cheb3Prefs.chebfun3eps;

% Example from Page 233 of: 
% Hackbusch, Tensor Spaces and Numerical Tensor Calculus, Springer, 2012.
% Make also sure that singular values are represeted as columns, not rows.
f = chebfun3(@(x,y,z) x.*z + x.^2.*y, [0, 1, 0, 1, 0, 1]);
sv = hosvd(f);
exactsv1 = [sqrt((109/720)+(sqrt(46)/45)); sqrt((109/720)-(sqrt(46)/45))];
exactsv2 = [sqrt((109/720)+(sqrt(2899)/360)); sqrt((109/720)-(sqrt(2899)/360))];
exactsv3 = exactsv2;

pass(1) = norm(sv{1} - exactsv1) < tol;
pass(2) = norm(sv{2} - exactsv2) < tol;
pass(3) = norm(sv{3} - exactsv3) < tol;


% Check case of two outputs, by inspecting the _all orthogonality_ of the
% core tensor of second output:
f = chebfun3(@(x,y,z) x.*z + x.^2.*y);
[sv, g] = hosvd(f);
pass(4) = norm(squeeze(g.core(1,:,:) .* g.core(2,:,:))) < tol;
pass(5) = norm(squeeze(g.core(:,1,:) .* g.core(:,2,:))) < tol;
pass(6) = norm(g.core(:,:,1) .* g.core(:,:,2)) < tol;

% Check case of five outputs:
[sv, Score, Scols, Srows, Stubes] = hosvd(f);
pass(7) = norm(squeeze(Score(1,:,:) .* Score(2,:,:))) < tol;
pass(8) = norm(squeeze(Score(:,1,:) .* Score(:,2,:))) < tol;
pass(9) = norm(Score(:,:,1) .* Score(:,:,2)) < tol;

end