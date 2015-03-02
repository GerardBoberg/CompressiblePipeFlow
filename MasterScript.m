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

%% Isentropic Area change verification case

%% Fanno Flow verification case
% This is an example, done in class week 4, extended out to the sonic point
%
% There is a pipe of constant D = 0.15 meters, with f = 0.005
%   with inlet conditions of M1 = 0.3, P1 = 1atm, T1 = 273Kelvin
% In the class example, the pipe is 30 meters long. 
% Here, I will extend it out to L*, 33.97 meters, and a bit extra to 40 meters

% Parameters for Fanno verification
L    = 45;         % meters
span = [ 0, L ];   % meters
T_i  = 273;        % Kelvin
P_i  = 1;          % atm
Dfun = @ConstantDiameter; % constant diameter pipe
Dt   = 0.15;       % meters
f    = 0.005;       % friction factor
gamma= 1.4;        % ratio of specific heats for air

M_i = 0.3;
% Analytic Fanno Flow
[ x_out_ff, Mach_out_ff, T_ff, T0_ff, P_ff, P0_ff ] = AnalyticFannoFlow( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma );

% Method of beans to verify
[ x_out_fb, Mach_out_fb, T_fb, T0_fb, P_fb, P0_fb ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, true);

figure();
subplot( 1, 2, 1 ), plot( ...
        x_out_ff, Mach_out_ff, x_out_ff, T_ff,...
        x_out_ff, T0_ff, x_out_ff, P_ff, x_out_ff, P0_ff );
axis( [span(1), span(end), 0, 1 ] );
title( 'Pure Fanno Flow' )
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );
subplot( 1, 2, 2 ), plot( ...
        x_out_fb, Mach_out_fb, x_out_fb, T_fb, x_out_fb, T0_fb,...
        x_out_fb, P_fb, x_out_fb, P0_fb );
axis( [span(1), span(end), 0, 1 ] );
title( 'Method of Beans, no area change' )
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );
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
f    = 0.005;      % friction factor
gamma= 1.4;        % ratio of specific heats for air

% Graph the nozzle and choke point
% The choke point function is discrete. needs discrete inputs.
n           = 100001; % increase for more percision. Needs to be odd.
x_arr_choke = linspace( span(1), span(end), n );
D_arr_choke = Dfun( x_arr_choke, Dt, L );

[ x_choke, ~, ~ ] = FindChokePoint( x_arr_choke, D_arr_choke, f, gamma );
figure();
PlotNozzle( x_arr_choke, D_arr_choke );
hold on;
plot( [ x_choke, x_choke ], [-1.5, 1.5], 'r' );
axis( [0,2,-1.5,1.5] );

%% -- Case 1: purely subsnoic flow
% With friction
M_i   = 0.060;
ss_ac = false;
[ x_out_1, Mach_out_1, T_1, T0_1, P_1, P0_1 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);


%--- Without friction
M_i   = 0.060;
ss_ac = false;
[ x_out_1f, Mach_out_1f, T_1f, T0_1f, P_1f, P0_1f ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                0, gamma, ss_ac);
% Plot them both on a subplot
figure();
subplot( 1, 2, 1 ), plot(...
      x_out_1, Mach_out_1, x_out_1, T_1,...
      x_out_1, T0_1, x_out_1, P_1, x_out_1, P0_1 );
axis( [ span(1), span(end), 0, 1] );
title( 'Purely subsonic flow with friction' );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South');
xlabel( 'x, meters' );
ylabel( 'normalized value' );

subplot( 1, 2, 2 ), plot(...
      x_out_1f, Mach_out_1f, x_out_1f, T_1f,...
      x_out_1f, T0_1f, x_out_1f, P_1f, x_out_1f, P0_1f );
axis( [span(1), span(end), 0, 1] );
title( 'Purely subsonic flow, frictionless' );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );
xlabel( 'x, meters' );
ylabel( 'normalized value' );


display( '----- Case 1 -----' )
display( 'Exit condition thermo properties = X * isentropic values' );
display( ['M_friction  = ', num2str( (1 - Mach_out_1(end)/Mach_out_1f(end)) ), ' * M_isen' ] );
display( ['P_friction  = ', num2str( (1 - P_1(end)/P_1f(end) )), ' * P_isen' ] );
display( ['P0_friction = ', num2str( (1 - P0_1(end)/P0_1f(end) )), ' * P0_isen' ] );
display( ['T_friction  = ', num2str( (1 - T_1(end)/T_1f(end) )), ' * T_isen' ] );
display( ['T0_friction = ', num2str( (1 - T0_1(end)/T0_1f(end) )), ' * T0_isen' ] );
                                            
%% -- Case 2: choked flow, subsonic after choke
% With friction
M_i   = 0.062;
ss_ac = false;
[ x_out_2, Mach_out_2, T_2, T0_2, P_2, P0_2 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);

%figure();
%dMdx_2 = diff( Mach_out_2 ) ./ diff( x_out_2 );
%plot( x_out_2(2:end), dMdx_2 );
%title( 'dMdx, choked subsonic' );


%--- Without friction
M_i   = 0.0644;
ss_ac = false;
[ x_out_2f, Mach_out_2f, T_2f, T0_2f, P_2f, P0_2f ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                0, gamma, ss_ac);

figure();
subplot( 1, 2, 1 ), plot(...
      x_out_2, Mach_out_2, x_out_2, T_2, x_out_2, T0_2,...
      x_out_2, P_2, x_out_2, P0_2 );
axis( [ span(1), span(end), 0, 1] );
title( 'Choked subsonic exit, with friction' );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );
xlabel( 'x, meters' );
ylabel( 'normalized value' );
subplot( 1, 2, 2 ), plot(...
      x_out_2f, Mach_out_2f, x_out_2f, T_2f, x_out_2f, T0_2f,...
      x_out_2f, P_2f, x_out_2f, P0_2f );
title( 'Choked subsonic exit, frictionless' );
axis( [ span(1), span(end), 0, 1] );
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'South' );
xlabel( 'x, meters' );
ylabel( 'normalized value' );


display( '----- Case 2 -----' )
display( 'Exit condition thermo properties = X * isentropic values' );
display( ['M_friction  = ', num2str( 1 - Mach_out_2(end)/Mach_out_2f(end) ), ' * M_isen' ] );
display( ['P_friction  = ', num2str( 1 - P_2(end)/P_2f(end) ), ' * P_isen' ] );
display( ['P0_friction = ', num2str( 1 - P0_2(end)/P0_2f(end) ), ' * P0_isen' ] );
display( ['T_friction  = ', num2str( 1 - T_2(end)/T_2f(end) ), ' * T_isen' ] );
display( ['T0_friction = ', num2str( 1 - T0_2(end)/T0_2f(end) ), ' * T0_isen' ] );
                                            
%% -- Case 3: choked flow, supersonic after choke
% With friction
M_i   = 0.062;
ss_ac = true;
[ x_out_3, Mach_out_3, T_3, T0_3, P_3, P0_3 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac);


%--- Without friction
M_i   = 0.06447;
ss_ac = true;
[ x_out_3f, Mach_out_3f, T_3f, T0_3f, P_3f, P0_3f ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                0, gamma, ss_ac);

figure();
subplot( 1, 2, 1 ), plot(...
      x_out_3, Mach_out_3, x_out_3, T_3, x_out_3, T0_3,...
      x_out_3, P_3, x_out_3, P0_3 );
title( 'Choked supersonic exit, with friction' )
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'NorthWest' )
xlabel( 'x, meters' );
ylabel( 'normalized value' );
subplot( 1, 2, 2 ), plot(...
      x_out_3f, Mach_out_3f, x_out_3f, T_3f, x_out_3f, T0_3f,...
      x_out_3f, P_3f, x_out_3f, P0_3f );
title( 'Choked supersonic exit, frictionless' )
legend( 'Mach number', 'T / T_i', 'T_0 / T_0_,_i', 'P / P_i', 'P / P_0_,_i',...
        'Location', 'NorthWest' )
xlabel( 'x, meters' );
ylabel( 'normalized value' );


display( '----- Case 3 -----' )
display( 'Exit condition thermo properties = X * isentropic values' );
display( ['M_friction  = ', num2str( 1 - Mach_out_3(end)/Mach_out_3f(end) ), ' * M_isen' ] );
display( ['P_friction  = ', num2str( 1 - P_3(end)/P_3f(end) ), ' * P_isen' ] );
display( ['P0_friction = ', num2str( 1 - P0_3(end)/P0_3f(end) ), ' * P0_isen' ] );
display( ['T_friction  = ', num2str( 1 - T_3(end)/T_3f(end) ), ' * T_isen' ] );
display( ['T0_friction = ', num2str( 1 - T0_3(end)/T0_3f(end) ), ' * T0_isen' ] );

