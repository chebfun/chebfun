function f = laplacian( f ) 
% LAPLACIAN     Scalar laplacian of a spherefun 
% 
% F = LAPLACIAN( F )

realf = isreal(f);

% TODO: Handle the case where the user selects latitude instead of
% co-latittude.

[cols,D,rows] = cdr(f);

% We are going to work at the tech level to make things faster.
coltechs = cols.funs{1}.onefun;
rowtechs = rows.funs{1}.onefun;

% Compute the derivatives
d_coltechs = diff(coltechs)/pi;
dd_coltechs = diff(d_coltechs)/pi;
dd_rowtechs = diff(rowtechs,2)/pi^2;

% We will do everyting in value space. Sorry Alex.

% Evaluate at the half grid points in theta, so the poles are not included.
m = length(cols);
h = 2*pi/m;
shift = h/2/pi;

coltechs = circshift(coltechs,-shift);
coltech_vals = coltechs.values;
d_coltechs = circshift(d_coltechs,-shift);
d_coltech_vals = d_coltechs.values;
dd_coltechs = circshift(dd_coltechs,-shift);
dd_coltech_vals = dd_coltechs.values;

% Evaluate sin(theta) and cos(theta) terms at the half grid points.
th = trigpts(m,[-1,1]) + h/2/pi;
sinth = sin(pi*th);
costh = cos(pi*th);

% Evalute rows and the derivatives at the grid
rowtech_vals = rowtechs.values.';
dd_rowtech_vals = dd_rowtechs.values.';

f = dd_coltech_vals*D*rowtech_vals + bsxfun(@rdivide,...
    bsxfun(@times,d_coltech_vals*D*rowtech_vals,costh) + ...
    bsxfun(@rdivide,coltech_vals*D*dd_rowtech_vals,sinth),sinth);

% Without factoring 1/sin(th) out
% f = dd_coltech_vals*D*rowtech_vals + ...
%     bsxfun(@times,d_coltech_vals*D*rowtech_vals,costh./sinth) + ...
%     bsxfun(@rdivide,coltech_vals*D*dd_rowtech_vals,sinth.^2);
% 

% Factoring out sin(th)^2 term.
% f = bsxfun(@rdivide,bsxfun(@times,dd_coltech_vals*D*rowtech_vals,sinth.^2) + ...
%     bsxfun(@times,d_coltech_vals*D*rowtech_vals,costh.*sinth) + ...
%     coltech_vals*D*dd_rowtech_vals,sinth.^2);


%
% Shift back to regular grid points in theta.
%

% Just do the works ourselves rather than call trigtech as the code will be
% faster and we know how to easily shift back by h/2.

n = length(rows);
idcol1 = m:-1:m/2+1;
idcol2 = m/2:-1:1;
% Enforce exact symmetry
f = 0.5*(f + [f(idcol1,n/2+1:n) f(idcol1,1:n/2) ; f(idcol2,n/2+1:n) f(idcol2,1:n/2)]);

% Wave numbers ordering accoriding to MATLAB's FFT
m1 =  floor((m-1)/2);
m2 = (m/2)*ones(rem(m+1,2));
waveNum = [(0:m1)  m2 (-m1:-1)]';

f = ifft(bsxfun(@times,fft(f),exp(1i*shift*pi).^waveNum));

if realf
    f = real(f);
end

f(m/2,:) = mean(f(m/2,:));
f(m,:) = mean(f(m,:));

f = spherefun(f(m/2:m,:));

end 