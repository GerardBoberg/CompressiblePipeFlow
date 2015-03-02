function [ x_out, Mach_out, T, T0, P, P0 ] = MethodOfBeans( span,...
                                                M_i, T_i, P_i,...
                                                Dfun, Dt, L,...
                                                f, gamma, ss_ac)
%MethodOfBeans Numerically solves a compressible 1-D flow situation. 
%   
%   Assumes the flow is choked, iterates backwards from the choke point to
%     the inlet condition to see if the assumption was correct.
%   If it is choked, then iterate forwards from the choke point to the exit
%     conditions, for both super sonic and subsonic exits.
%   If it isn't choked, then iterate forwards from the inlet to the exit.
%
%   span --- a 2-long array of the area to iterate over. 
%               [0, 2L] is expected.
%   M_i   --- The Mach number at the inlet, the initial condition for Mach
%   T_i   --- The Temperature at the inlet, the initial condition for Temp
%   P_i   --- The Pressure at the inlet, the initial condition for Pressure
%   Dfun  --- A function handle of the form F( x[], Dt, L ) that returns the
%                diameter of the pipe at any location x.
%   Dt    --- The throat diameter. Used for Dfun.
%   L     --- The characteristic length. Used for Dfun.
%   f     --- Friction factor. 0 for frictionless, 0.005 small, 0.4 large
%   gamma --- The ratio of specific heats for the gas. 1.4 for air.
%   ss_ac --- Is the flow super sonic after the choke point? true / false


%% Step 1: assume choked, Find the choke point and Singular roots
% Start by finding the choke point
% The choke point function is discrete. needs discrete inputs.
n           = 100001; % increase for more percision. Needs to be odd.
x_arr_choke = linspace( span(1), span(end), n );
D_arr_choke = Dfun( x_arr_choke, Dt, L );

% The nozzle has constant diameter. Give a false x_choke around 0+delta
%  so that it will determine that it is not choked and propogate forwards.
if( max( abs( diff( D_arr_choke ) ) ) <= 1e-8 )
    delta_x = ( span(end) - span(1) ) / 1000;
    x_choke = span(1) + delta_x;
else
    [ x_choke, ~, ~ ] = FindChokePoint( x_arr_choke, D_arr_choke, f, gamma );
end

[ M_plus, M_minus ] = SolveMachSingularRoots( x_choke, Dfun, Dt,...
                                              L, f, gamma );


%% Step 2: assume choked, interate backwards to check our assumption
% Now iterate backwards from the choke point to the inlet.

% The differential Mach equation is singular at the choke point, so the
%   first step we pre-calculate using a linear assumption.
% Start by finding the step size required to get past Mach .995
% 1 - (B * dx) = t;
% 1 - t        = B * dx
%(1 - t) / B   = dx
target_step = .995;
delta_x  = ( 1 - target_step ) / M_plus;
delta_x  = abs( delta_x );
c2i_span = [ x_choke - delta_x , 0];  % One step before M=1
c2i_M    = target_step;

% Propogate 
[ x_back, Mach_back ] = PropogateMachNumber( c2i_span, c2i_M,...
                                                  Dfun, Dt, L,...
                                                  f, gamma,...
                                                  ss_ac, x_choke,...
                                                  M_plus, M_minus );
                                              
%% Step 3: Check our choked assumption
% Case A: Flow is not choked. Propogate upstream from the inlet.
% Case B: Flow is choked. Run super or subsonic upstream analysis.
% Case C: Non-physical inlet condition. Shockwave will occur, that is not
%           modeled.
M_check = Mach_back( end );
tol_chocked_assumption = 1e-3;
%is_choked = false;
% Case B -- perfectly choked.
if( abs( M_check - M_i ) < tol_chocked_assumption )
    is_choked = true;


% Case A -- Flow is never choked.
elseif( M_check > M_i )
    is_choked = false;

% Case C -- Error out. We're not handling supersonic shockwaves.
else%if( M_check < M_i )
    e_message = [ 'Inlet Condition is non-physical! Maximum physical',...
                    'inlet condition is Mach_i = ', num2str(M_check),...
                    '. M_i given = ', num2str( M_i ) ];
    error( e_message );
end

%% Step 4: Acquire full Mach-number variation.
%
% Case A -- Flow is NOT choked. Simply propogate forward from inlet.
if( is_choked == false )
    
    % Propogate from inlet conditions
    [ x_out, Mach_out ] = PropogateMachNumber( span, M_i,...
                                                  Dfun, Dt, L,...
                                                  f, gamma,...
                                                  ss_ac, x_choke,...
                                                  M_plus, M_minus );
    
% Case B -- Flow IS choked. Use the propogation data from step 2, and
%               append the post-sonic-point propogation.
else % if( is_choked == true )
    
    % The differential Mach equation is singular at the choke point, so the
    %   first step we pre-calculate using a linear assumption.
    % Start by finding the step size required to get past Mach 1.005
    if( ss_ac == true ) % super sonic_after choke
        target_step = 1.005;
        delta_x = ( 1 - target_step ) / M_plus;
    else
        target_step = 0.995;
        delta_x = ( 1 - target_step ) / M_minus;
    end
    
    delta_x  = abs ( delta_x );
    c2i_span = [ x_choke + delta_x , span(2)];  % one step beyond M=1
    c2i_M    = target_step;

    % Propogate 
    [ x_fwd, Mach_fwd ] = PropogateMachNumber( c2i_span, c2i_M,...
                                                      Dfun, Dt, L,...
                                                      f, gamma,...
                                                      ss_ac, x_choke,...
                                                      M_plus, M_minus );

    % Take our propogation and place it into the output variables
    x_out = x_back(end:-1:1); % x_back is backwards. flip it around.
    x_out = [ x_out; x_choke; x_fwd ]; % appened back, M=1, fwd
    
    Mach_out = Mach_back(end:-1:1); % Mach_back is also backwards. flip.
    Mach_out = [ Mach_out; 1; Mach_fwd ]; % append back, M=1, fwd
   
end

%% Step 5: Calculate Thermodynamic properties
% We now have full Mach number distribution. Thermo properties are
%   functions of Mach number, and the area at that point.

Afun = @(x)( (pi/4) * Dfun(x, Dt, L).^2 ); % Lambda function for the area
[ T, P, T0, P0 ] = MethodOfBeansThermo( x_out, Mach_out, Afun,...
                                              T_i, P_i, gamma );

% Normalize outputs of thermo functions by the inlet conditions
T  = T  ./ T(1);
P  = P  ./ P(1);
T0 = T0 ./ T0(1);
P0 = P0 ./ P0(1);

% That's everything.
end

