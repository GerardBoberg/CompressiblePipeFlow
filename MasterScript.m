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

%% -- Case 1: purely subsnoic flow
% With friction
M_i   = 0.062;
ss_ac = false;
[ x_out_1, Mach_out_1, T_1, T0_1, P_1, P0_1 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);

figure();
plot( x_out_1, Mach_out_1, x_out_1, T_1, x_out_1, T0_1, x_out_1,...
      P_1, x_out_1, P0_1 );
title( 'Purely subsonic flow with friction' );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South');

%--- Without friction
M_i   = 0.062;
ss_ac = false;
[ x_out_1f, Mach_out_1f, T_1f, T0_1f, P_1f, P0_1f ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                0, gamma, ss_ac);

figure();
plot( x_out_1f, Mach_out_1f, x_out_1f, T_1f, x_out_1f, T0_1f, x_out_1f,...
      P_1f, x_out_1f, P0_1f );
axis( [span(1), span(end), 0, 1] );
title( 'Purely subsonic flow, frictionless' );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );


display( '----- Case 1 -----' )
display( 'Exit condition thermo properties = X * isentropic values' );
display( ['M_friction  = ', num2str( Mach_out_1(end)/Mach_out_1f(end) ), ' * M_isen' ] );
display( ['P_friction  = ', num2str( P_1(end)/P_1f(end) ), ' * P_isen' ] );
display( ['P0_friction = ', num2str( P0_1(end)/P0_1f(end) ), ' * P0_isen' ] );
display( ['T_friction  = ', num2str( T_1(end)/T_1f(end) ), ' * T_isen' ] );
display( ['T0_friction = ', num2str( T0_1(end)/T0_1f(end) ), ' * T0_isen' ] );
                                            
%% -- Case 2: choked flow, subsonic after choke
% With friction
M_i   = 0.064;
ss_ac = false;
[ x_out_2, Mach_out_2, T_2, T0_2, P_2, P0_2 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);

figure();
plot( x_out_2, Mach_out_2, x_out_2, T_2, x_out_2, T0_2, x_out_2,...
      P_2, x_out_2, P0_2 );
title( 'Choked flow, with subsonic exit, with friction' );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );

figure();
dMdx_2 = diff( Mach_out_2 ) ./ diff( x_out_2 );
plot( x_out_2(2:end), dMdx_2 );
title( 'dMdx, choked subsonic' );


%--- Without friction
M_i   = 0.064;
ss_ac = false;
[ x_out_2f, Mach_out_2f, T_2f, T0_2f, P_2f, P0_2f ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                0, gamma, ss_ac);

figure();
plot( x_out_2f, Mach_out_2f, x_out_2f, T_2f, x_out_2f, T0_2f, x_out_2f,...
      P_2f, x_out_2f, P0_2f );
title( 'Choked flow, with subsonic exit, frictionless' );
axis( [ span(1), span(end), 0, 1] );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );


display( '----- Case 2 -----' )
display( 'Exit condition thermo properties = X * isentropic values' );
display( ['M_friction  = ', num2str( Mach_out_2(end)/Mach_out_2f(end) ), ' * M_isen' ] );
display( ['P_friction  = ', num2str( P_2(end)/P_2f(end) ), ' * P_isen' ] );
display( ['P0_friction = ', num2str( P0_2(end)/P0_2f(end) ), ' * P0_isen' ] );
display( ['T_friction  = ', num2str( T_2(end)/T_2f(end) ), ' * T_isen' ] );
display( ['T0_friction = ', num2str( T0_2(end)/T0_2f(end) ), ' * T0_isen' ] );
                                            
%% -- Case 3: choked flow, supersonic after choke
% With friction
M_i   = 0.064;
ss_ac = true;
[ x_out_3, Mach_out_3, T_3, T0_3, P_3, P0_3 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);

figure();
plot( x_out_3, Mach_out_3, x_out_3, T_3, x_out_3, T0_3, x_out_3,...
      P_3, x_out_3, P0_3 );
title( 'Choked flow, with supersonic exit, with friction' )
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'NorthWest' )

%--- Without friction
M_i   = 0.064;
ss_ac = true;
[ x_out_3f, Mach_out_3f, T_3f, T0_3f, P_3f, P0_3f ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                0, gamma, ss_ac);

figure();
plot( x_out_3f, Mach_out_3f, x_out_3f, T_3f, x_out_3f, T0_3f, x_out_3f,...
      P_3f, x_out_3f, P0_3f );
title( 'Choked flow, with supersonic exit, frictionless' )
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'NorthWest' )


display( '----- Case 3 -----' )
display( 'Exit condition thermo properties = X * isentropic values' );
display( ['M_friction  = ', num2str( Mach_out_3(end)/Mach_out_3f(end) ), ' * M_isen' ] );
display( ['P_friction  = ', num2str( P_3(end)/P_3f(end) ), ' * P_isen' ] );
display( ['P0_friction = ', num2str( P0_3(end)/P0_3f(end) ), ' * P0_isen' ] );
display( ['T_friction  = ', num2str( T_3(end)/T_3f(end) ), ' * T_isen' ] );
display( ['T0_friction = ', num2str( T0_3(end)/T0_3f(end) ), ' * T0_isen' ] );

