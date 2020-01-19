% fourier analysis and SIF comparison

function [Kv, erv] = mapping(va, v1, filt)

% removing the zero values and making identical analy. field
v0 = ones( size( va ) );
v0( va == 0 ) = 0;
v1(va == 0) = 0;

% considering mean here changes the process a lot.
% mean could be mean( va(:)) = a number
% or mean( va, 1 ) = an array. array is sturdy and stable.
siz_v = size(va);
vm = ones( siz_v(1), 1 ) * mean( va, 1 );
va = va - vm;

% calculation of the weighted average, weight by the abs val
ord = 0.5;
vw = sum( va .* abs(va).^ord, 1 );
vws = sum( (abs(va)).^ord, 1 );
vws( vws == 0 ) = 1;
vw = vw ./ vws;
va = va - ones( siz_v(1), 1 ) * vw;

% fourier transformation
fv = fft2(va);
fv1 = fft2(v1);
fvs = sort( abs(fv1(:)) );

opts = optimset('Display','off');

[i, j] = find( abs(fv1) > fvs( ceil( numel(fv)*filt ) ) );
fvm = zeros( size(fv) );
fvm( i, j ) = fv( i, j );
fun = @(c, x) c(1)*x(:,1);

% least square error evaluation during the mapping
% mapping includes complex numbers (frequency)
x = fv1(:);
[Kvt, ~] = lsqcurvefit( fun, 0, x, fvm(:),[],[], opts);
erv = sum( sum( ( abs( real(Kvt(1)*v1) - va ) ).^0.5 ) );
Kv = abs( Kvt );

end