% To find the crack tip & K for a defined path

clc, clear
load 'data'  % contains 'crackpath (path) and displ. (Uy)

filt = 0.5;                 % filter in terms of %
m = 10;                     % cut out of strip of width 2m-1
rat = 7.6e-5;               % pixel to meter ration
len = 100;                  % length of the domain of inter.
wid = 60;                   % width of the domain
dom = 0.2;                  % horizontal offset of the domain
c = 453.5;                  % crack prop. velocity
E = 6e9;                    % Young's modulus
nu = 0.33;                  % Poisson's ratio
rho = 1700;                 % density
crtip = [288, 457];         % previous cracktip/notch (pixel)

Ki = zeros(length(path), 1);         % SIF declaration
err = zeros(length(path), 1);         % error distribution

for i = 1 : length(path)
% selection of the domain of interest
rect = [crtip(1,1)-dom*len, crtip(1,2)-wid/2, len, wid];
v = uy( rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3) );
domsiz = [rect(4)+1, rect(3)+1];
% creating the analytical domain
yspan = linspace(0, domsiz(1)*rat, domsiz(1) );
xspan = linspace(0, domsiz(2)*rat, domsiz(2) );
[xd, yd] = meshgrid( xspan, yspan );
xdn = xd - ( path(i, 1) - rect(1) ) * rat;
ydn = yd - ( path(i, 2) - rect(2) ) * rat;
cent = round(path(i, :)) + [1 1] - rect(1:2);
% calling the function to create analytic. field and map
[Kv, er] = field(v, xdn, ydn, filt, cent, m, c, E, nu, rho);
Ki(i) = abs( Kv(1) );
err(i) = er;
end

figure(1), plot(err)
xlabel('crack path (pixel)'), ylabel('Error distribution')

[~, k1] = min(err);         % point corresponding to min. err.
K = Ki(k1);                 % SIF at the instance
cr_tip = path(k1, :);       % crack tip for this instance

% display of the crack path and tip
disp(['Crack tip is ([x, y] in pixels) = ', num2str(cr_tip)])
disp(['SIF = ' num2str(K)]);

figure(2), imshow('image.jpg'), hold on
plot(path(:, 1), path(:, 2), 'g.')
plot(cr_tip(1), cr_tip(2), 'r*')
disp(['The crack tip is shown in the image. The slight ',...
    'offset in the path is due to the rigid body motion']);