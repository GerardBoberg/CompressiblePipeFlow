%% CompressiblePipeFlow 
% Written by Gerard Boberg and Zane Paterson
%   for the class Aero 303, taught by Dr. Marshal at Cal Poly, SLO
%

clc;
clear all;
close all;

%% Project description:
% The purpose of this project is to model a compressible gas flow thru a
%   pipe with a converging-diverging nozzle.
%
% Friction and Isentropic Area change are modeled.
% Heat exchange and Shockwaves are NOT modeled.
%
% In addtion, pure isentropic area change and pure fanno flow methods are
%   also provided.

%% Step 1: find the choke point
% At this point we ASSUME that the flow is choked. 
% If it turns out that it's not, we'll deal with that later.

f = 0.005;   % given friction factor
gamma = 1.4; % gamma of air

Dt = 0.05; % meters, Diameter of the Throat 
L  = 2;    % meters, Length of nozzle
n = 10001;  % number of steps to numerically solve with
x = linspace( 0, L, n );

% Calculate the diameter variation over our geometry
D = ConDiNozzleDiameter( x, Dt, L / 2.0 ); % The equation for the nozzle expects L/2

% Find the choke point
[ choke_location, choke_index, ~ ] = FindChokePoint( x, D, f, gamma );

figure();
PlotNozzle( x, D );
hold on;
plot( [choke_location, choke_location], [ 1.5, -1.5], 'r' );
axis( [0, 2, -1.5, 1.5] );
display( ['The choke location is at: ', num2str( choke_location ), ' meters'] );

%% Step 2: Compare to 1D Fanno Flow
%
% Fanno flow only has friction, with no area change.

% Start with a constant diameter pipe
D_f = ConstantDiameter( x, Dt * 3 );

%% Step 3: Run the analysis 3 times for 3 different inlet conditions


L     = 1;         % meters, 1/2 of span
span  = [0, 2*L];  % meters
init  = 0.05;      % initial mach number
Dfun  = @ConDiNozzleDiameter;
Dt    = 0.05;  % meters
f     = 0.005; % friction factor
gamma = 1.4;   % ratio of specific heats for air
ss_ac = false;

% n       = 10001;
% x_array = linspace( span(1), span(end), n ); 
% D_array = Dfun( x_array, Dt, L );
% [ x_choke, ~, ~ ] = FindChokePoint( x, D, f, gamma );
x_choke = choke_location;

[ M_plus, M_minus ] = SolveMachSingularRoots( x_choke, Dfun, Dt, L, f, gamma )

% 1 - (B * dx) = t;
% 1 - t        = B * dx
%(1 - t) / B   = dx
target_step = .995;
delta_x = ( 1 - target_step ) / M_plus;
invinit = target_step;

invspan = [ x_choke - delta_x , 0 ];
[ x_out, Mach ] = PropogateMachNumber( invspan, invinit, Dfun, Dt, L, f, gamma,...
                                            ss_ac, x_choke, M_plus, M_minus );

figure();
plot( x_out, Mach );
ts = strcat( 'M_i = ', num2str( init ), ' // Backwards Propogation ' );
title( ts );


%[ x_out, Mach ] = PropogateMachNumber( span, init, Dfun, Dt, L, f, gamma,...
%                                            ss_ac, x_choke, M_plus, M_minus );

figure();
plot( x_out, Mach );
ts = strcat( 'M_i = ', num2str( init ), ' //  Forwards Propogation ' );
title( ts );



