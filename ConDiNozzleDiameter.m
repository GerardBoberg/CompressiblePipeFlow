function [ diameter ] = ConDiNozzleDiameter( x, L, D_t )
%ConDiNozzleDiameter The diameter of the converging-diverging nozzle at x.
%   x   --- A vector of the x locations of the nozzle diameter.
%   D_t --- The diameter of the nozzle at the throat.
%             Defaults to 0.05 meters
%   L   --- Half of the length of the nozzle. 
%             Defaults to half the maximum value of x.


% error check inputs
if( nargin <= 1 )
    L = 0.5 * max( x );  % default Length
end
if( nargin == 2 )
    D_t = 0.05;  % default throat diameter
end

% actual equation for the nozzle geometry
%   varies between 1 and 3 times throat diameter.
diameter = 2 * D_t .* ( 1 + ( 1/2 * cos( (pi/L) .* x) ) );

end

