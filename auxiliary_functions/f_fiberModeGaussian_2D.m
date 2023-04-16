function [E, H] = f_fiberModeGaussian_2D( w0, lambda0, xvec, zvec, theta, d0, nclad )
% somewhat adapted from Cale's code
%
% Generate Gaussian-beam mode profile at a plane through y = d0 at angle of
% theta (in degrees)
%
% the coordinate system is a bit confusing - z = longitudinal direction of grating =
% direction the fiber is tilted in
% x = the other transverse dimension of the fiber
%
% NOTE that the angle theta is with respect to the positive y axis. When
% overlapping downwards, a "positive theta" in our typical definition (from
% the downwards normal) is actually a negative theta with respect to the
% positive y axis.
%
% currently solving for TE mode only (Ex polarized)
%
% (using H.A. Haus, Waves & Fields, Chapter 5)
% CMG November 21st, 2014
%
%
% Inputs:   
%   w0  
%       type: double, scalar
%       desc: 1/e beam RADIUS at waist, in units 'units' (MFD/2)         
%   lambda0
%       type: double, scalar
%       desc: free space wavelength, in units 'units'
%   xvec
%       type: double, array
%       desc: x coordinates, n by 1
%   zvec
%       type: double, array
%       desc: z coordinates, n by 1    
%   theta 
%       type: double, scalar
%       desc: angle from normal in degrees (in the yz plane)
%   d0  
%       type: double, scalar
%       desc: distance from beam waist to slice
%   nclad 
%       type: double, scalar
%       desc: cladding index
%
% Outputs: 
%   E - struct with x, y, z components
%   H - struct with x, y, z components


% Constants, in units of meters
c       = 3e8;                                  % m/s
mu0     = 4*pi * 1e-7;                          % H/m
% omega0  = 2*pi*c/( lambda0 * unit_scale);       % rad/s.
lambda  = lambda0 / nclad;         % wavelength in cladding
k       = 2*pi / lambda;                          % 1/m
k0      = 2*pi / lambda0;
% w0      = w0;                      % radius
% d0      = d0;                      % offset


% Convert to radians
theta = (pi/180)*theta;

% Scale coordinates
% xvec = xvec;                                              % units m
% zvec = zvec;                                              % units m

% coordinates in fiber frame
zprime = zvec.*cos(theta) - d0*sin(theta);
yprime = zvec.*sin(theta) + d0*cos(theta);  % normal direction
xprime = xvec;
[ Zprime, Xprime ] = meshgrid( zprime, xprime );    % dimensions x' vs z'
if isrow( yprime )
    Yprime = repmat( yprime, [ length(xprime), 1 ] );
else
    Yprime = repmat( yprime.', [ length(xprime), 1 ] );
end

% b (confocal parameters) is used instead of z0 so that z0 = -1j.*b removes the singularity of the solution on the real z axis (see Haus pg 109)
b = k*w0^2/2;                                                                                   

% Equation (5.2) in Haus [1/meters]
u00 =   1j .* sqrt(k*b/pi) .* ( 1./(Yprime + 1j.*b) ).*...
    exp( -1j.*k.*( Xprime.^2 + Zprime.^2 )./( 2*(Yprime + 1j.*b) ) );     

% normalize the mode to intensity, makes things nicer
% dx      = xvec(2) - xvec(1);
% u00     = u00/sqrt( dx * sum( abs( u00 ).^2 ) );

% calculate fields, the suffix is the polarization
E_xprime = -1i * u00 .* exp( -1i * k * Yprime );
E_yprime = - ( Xprime ./ ( Yprime + 1i * b ) ) .* E_xprime;
% E_yprime = 1i * ( Xprime ./ ( Yprime + 1i * b ) ) * u00 .* exp( -1i * k * yprime );
H_zprime = ( nclad ./ ( mu0 * c ) ) .* E_xprime;
H_yprime = - ( -Zprime./( Yprime + 1i*b ) ) .* H_zprime;
% H_xprime = ( -1i * k/mu0 ) * u00 .* exp( -1i * k * yprime );
% H_zprime = ( -Zprime./( Yprime + 1i*b ) ) .* H_xprime; 

% project fields back into original coordinates
Ex = E_xprime;
Ey = E_yprime * cos(theta);
Ez = E_yprime * sin(theta);
Hx = zeros( size( Ex ) );
Hy = H_yprime * cos(theta) + H_zprime * sin(theta);
Hz = H_yprime * sin(theta) + H_zprime * cos(theta);

% return as structs
E = struct( 'x', Ex, 'y', Ey, 'z', Ez );
H = struct( 'x', Hx, 'y', Hy, 'z', Hz );


end     % end fiberModeGaussian()












































