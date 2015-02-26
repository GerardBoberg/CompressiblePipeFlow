function [ diameter ] = ConDiNozzleDiameter( x, diameter_throat, length )
%ConDiNozzleDiameter The diameter of the converging-diverging nozzle at x.
%   Detailed explanation goes here

diameter = 2 * diameter_throat * ( 1 + ( 1/2 * cos( (pi/length) * x) ) );

end

