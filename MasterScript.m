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

% General parameters that are consistent for all cases
L    = 1;          % meters
span = [ 0, 2*L ]; % meters
T_i  = 300;        % Kelvin
P_i  = 101.325;    % kPa
Dfun = @ConDiNozzleDiameter; % cos-based converging-diverging nozzle
Dt   = 0.05;       % meters
f    = 0.01;      % friction factor
gamma= 1.4;        % ratio of specific heats for air

% -- Case 1: purely subsnoic flow
M_i   = 0.062;
ss_ac = false;
[ x_out_1, Mach_out_1, T_1, T0_1, P_1, P0_1 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);

figure();
plot( x_out_1, Mach_out_1, x_out_1, T_1, x_out_1, T0_1, x_out_1,...
      P_1, x_out_1, P0_1 );
title( 'Purely subsonic flow' );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i' );
                                            
% -- Case 2: choked flow, subsonic after choke
M_i   = 0.064;
ss_ac = false;
[ x_out_2, Mach_out_2, T_2, T0_2, P_2, P0_2 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);

figure();
plot( x_out_2, Mach_out_2, x_out_2, T_2, x_out_2, T0_2, x_out_2,...
      P_2, x_out_2, P0_2 );
title( 'Choked flow, with subsonic exit' );

figure();
dMdx_2 = diff( Mach_out_2 ) ./ diff( x_out_2 );
plot( x_out_2(2:end), dMdx_2 );
hold on;
plot( [ 1.072, 1.072 ], [ max( dMdx_2), min( dMdx_2 ) ], 'r' );
hold off;
title( 'dMdx, choked subsonic' );
%legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i' );
                                            
% -- Case 3: choked flow, supersonic after choke
M_i   = 0.064;
ss_ac = true;
[ x_out_3, Mach_out_3, T_3, T0_3, P_3, P0_3 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);

figure();
plot( x_out_3, Mach_out_3, x_out_3, T_3, x_out_3, T0_3, x_out_3,...
      P_3, x_out_3, P0_3 );
title( 'Choked flow, with supersonic exit' )
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i' )