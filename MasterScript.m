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


%% Verification Cases
% Table of Mach# and Thermodynamic properties for:
%   1) Isentropic thru nozzle vs. Method of Beans
%   2) Fanno flow thru a pipe vs. Method of Beans


%% Example cases
%
%
% 1) purely subsonic flow  vs. isentropic same conditions
% 2) Choked flow, goes subsonic afterwards vs. isentropic choked
% 3) Choked flow, goes supersonic afterwards vs. isentropic choked