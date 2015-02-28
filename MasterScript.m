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

L  = 2;    % meters, Length of nozzle
Dt = 0.05; % meters, Diameter of the Throat 
n = 1001; % number of steps to numerically solve with
x = linspace( 0, L, n );

% Calculate the diameter variation over our geometry
D = ConDiNozzleDiameter( x, L / 2.0, Dt ); % The equation for the nozzle expects L/2

% Find the choke point
[ choke_location, choke_index, ~ ] = FindChokePoint( x, D, f, gamma );

figure();
PlotNozzle( x, D );
hold on;
plot( [choke_location, choke_location], [ 1.5*Dt, -1.5*Dt], 'r--' );
