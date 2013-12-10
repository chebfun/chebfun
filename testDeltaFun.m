deltafun

%%
% shouldn't work
deltafun(1)

%%
deltafun(1,2)
%%
% shouldn't work
deltafun( rand(3,1), rand(1,4) )
%%
% shouldn't work
deltafun( rand(3,1), rand(3,1) )
%%
% should work
deltafun( rand(3,3), rand(3,1) )
deltafun( rand(3,3), rand(1,3) )
%%
% shouldn't work
deltafun( rand(3,1), rand(1,3) )
%%
% should work
deltafun( rand(1,3), rand(3,1) )
%%
% Test simplify
mag = [rand(3,3); zeros(3,3)]; loc = rand(3,1); 
%%mag(end,end) = 1e-12;
d = deltafun(mag, loc, chebfun(0)); d.delta.magnitude
d = simplify(d); d.delta.magnitude
%%
% Dirac delta and derivatives
d = dirac(deltafun);
d = diff(d,5);
d.delta.magnitude
%%
% The dirac delta function and inner products
loc = 0;
mag = 1;
f = chebfun(0);
d = deltafun(mag, loc, f)
x = chebfun('x')
ip(d, x)
f = sin(x);
ip(diff(d), f )

%%
mag = rand(1,5);
loc = rand(1,5);
d = deltafun( mag, loc, chebfun(0));
ip( d, sin(2*pi*x))
%%
% Some more inner products
loc = rand(1, 5);
mag = rand(1, 5);
d = deltafun( mag, loc, chebfun(0))
ip(d, 1) - sum(mag)
%%
loc = [-.5 .5];
mag = [1 2];
x = chebfun('x');
d = deltafun( mag, loc, 0*x)
ip(d, x) - sum(mag.*loc)

% %%



%%
% tol = 1e-10;
% 
% va = rand(1, 5);
% vb = va;
% va = [ 1 2 3 4 ];
% vb = [ 2 3 4 5];
% 
% Au = rand(2,4);
% Bu = rand(2,4);
% 
% [va, idx] = sort(va);
% A = Au(:, idx);
% 
% [vb, idx] = sort(vb);
% B = Bu(:, idx);
% 
% 
% m = 1; n = 1; i = 1;
% lenA = length(va);
% lenB = length(vb);
% 
% vc = zeros(1, lenA + lenB);
% C  = zeros(size(A,1), length(vc)); % or zeros(size(B))
% 
% while ( m <= lenA || n <= lenB )
%     % If one of the arrays is exausted, copy the second one into output, no
%     % duplication is assumed in va or vb, so:
%     if ( m > lenA )
%         while ( n <= lenB )
%             vc(i) = vb(n);
%             C(:, i) = B(:, n);
%             n = n + 1;
%             i = i + 1;
%         end
%         break;
%     end
%     
%     if ( n > lenB )
%         while ( m <= lenA )
%             vc(i) = va(m);
%             C(:, i) = A(:, n);
%             m = m + 1;
%             i = i + 1;
%         end
%         break;
%     end
%     
%     % None of the arrays is exausted;
%     
%     % Duplication scenarios
%     p = ( va(m) == 0 && vb(n) == 0 );
%     vcMax = max(abs([va(m), vb(n)]));
%     p = p | ( abs((va(m)-vb(n)))/vcMax < tol );
%     % If there is a duplicate
%     if ( p )
%         C(:,i) = A(:, m) + B(:, n);
%         vc(i) = va(m);
%         m = m + 1;
%         n = n + 1;% or vb(n);
%     elseif ( va(m) < vb(n) )
%         vc(i) = va(m);
%         C(:, i) = A(:, m);
%         m = m + 1;
%     else
%         vc(i) = vb(n);
%         C(:, i) = B(:, n);
%         n = n + 1;
%     end
%     i = i + 1;
% end
% C = C(:, 1:i-1)
% vc = vc(1:i-1)
% 
% 
% 
% 
