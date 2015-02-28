function [ ] = PlotNozzle( x, D )
%PlotNozzle Creates a pretty graph of the nozzle
%
%   x --- the vector of x locations of the nozzle's geometry
%   D --- the Diameter at each x point

r = 0.5 .* D; % radius is half the diameter

hold on;
area( x,  r ); % plot the upper and lower bounds
area( x, -r );
ylabel( 'Radius, meters' );
xlabel( 'Length, meters' );
title( 'Nozzle shape' );
hold off;

end

