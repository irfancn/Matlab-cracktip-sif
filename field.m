% to compare and find the K

function [Kv, erv] = field(va,x,y,filt,cent,m,c,E,nu,rho)
K1 = 1;             % for mode I
K2 = 0;             % for mode II

% material properties and wave velocities
mu = E / (2*(1 + nu));
k = (3 - nu)/(1 + nu);    % plane stress
cs = sqrt(mu/rho);
cl = sqrt( (k+1)/(k-1) * mu/rho );

% coef. for the displacement field
b1 = sqrt( 1 - (c/cl)^2 ); b2 = sqrt( 1 - (c/cs)^2 );
D = 4*b1*b2 - (1 + b2*b2)^2;
B1 = (1 + b1*b1)/D; B2 = 2*b2/D;

[t1, r1] = cart2pol(x, b1*y);
[t2, r2] = cart2pol(x, b2*y);

h1 = 2*b1*b2 / (1 + b2*b2);
h2 = (1 + b2*b2)/2;

% analytical field
v = K1*B1 * ( -b1*sqrt(r1) .* sin(t1/2) + h1/b2*sqrt(r2) .* sin(t2/2) ) +...
    K2*B2 * ( b1*sqrt(r1) .* cos(t1/2) - h2/b2*sqrt(r2).*cos(t2/2) );

% removing the center stip from the displacment to reduce error
vmin = min( abs( v( 1 : cent(2)-m, : ) ) );
v = [v(1 : cent(2) - m, :); v(cent(2) + m:end, :)];
va = [va(1 : cent(2) - m, :); va(cent(2) + m:end, :)];

% calling the mapping - overlapping in the Fourier space
[Kv, erv] = mapping(va, v, filt);
Kv = Kv * mu/sqrt(2/pi);
