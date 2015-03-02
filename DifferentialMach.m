function [ dMdx ] = DifferentialMach( x, M, Dfun, Dt, L, f, gamma, ss_ac,...
                                                x_choke, M_plus, M_minus )
%DifferentialMach Calculates the differential Mach equation. Use with ode45
%
%   Specifically, the Mach equation for adiabatic flow is:
%
%   dM     2M + (gamma - 1)M^3    {  gamma M^2 f       dA/dx    }
%   --  = --------------------- * { ------------- - ----------- }
%   dx       2 ( 1 - M^2 )        {       2 D            A      }
%
%   Which is singular as M approaches 1. error out when (1-M^2) <= 0.0001
%
%   This function propogates that equation over the interval x.
%
%   x     --- The current x-location. In ode45, this is called 't'
%   M     --- The current Mach number. In ode45, this is called 'y'
%   Dfun  --- A function handle to the func used to generate pipe geometry.
%              Should return the diameter of the pipe at any x value.
%   Dt    --- Diameter input for Dfun.
%   L     --- Length input for Dfun.
%   f     --- Friction factor. Defaults to 0: isentropic flow.
%   gamma --- Ratio of specific heats of the gas. Defaults to 1.4
%   ss_ac --- Super Sonic_After Choke: 

%% Check inputs

if( nargin < 8 )
    ss_ac = true; % Assume super sonic_after choke in unspecified
end

if( nargin < 7 )
    gamma = 1.4; % gamma for air
end

if( nargin < 6 )
    f = 0; % friction factor for isentropic flow
end


% Check for the singular condition, dividing by zero at M = 1
% We need to do some funky math to handle it.

delta = 0.0098; % .9951 mach or 1.0051 mach
if( abs(1 - M^2 ) <= delta )
    % error( 'DifferentialMach is singular as M approaches 1!' );
        
    dMdx = M_plus;
    % less than Mach AND approaching choke
    if( (M < 1) && ( x < x_choke )  )
        dMdx = M_plus;
        
    % greater than Mach AND approaching choke
    elseif( ( M >= 1 ) && ( x < x_choke ) ) 
        dMdx = M_minus;
    
    % super sonic after choke, and after choke
    elseif( (ss_ac == true ) && ( x > x_choke) )
        dMdx = M_plus;
    
    % Not ss after choke, and after choke
    elseif( (ss_ac == false) && ( x >= x_choke ) )
        dMdx = M_minus;
    end
    
    %display( [ 'dMdx = ', num2str( dMdx ) ] );
    return; % This is a singular edge case. Div zero if don't return now.
end

%% Calculate Diameter and Area values from @Dfun

D = Dfun( x, Dt, L ); % Diameter
A = ( pi/4 ) * D.^2;  % Area, assuming circular cross section

% Find dA/dx for this point using a finite difference.
delta = 0.0001; % meters
range  = [ x-delta , x + delta ];
A_range = ( pi/4 ) * Dfun( range, Dt, L ).^2 ; 
dAdx = ( A_range(2) - A_range(1) ) /  ( 2*delta );

%% Calculate the Differential Mach Equation
mach_term =              2 * M + ( gamma - 1 ) * M^3;      % numerator over
mach_term = mach_term / (2 * ( 1 - M^2 ) );                % denominator

friction_term = gamma * M^2 * f / ( 2*D ) * 4;

area_term = dAdx / A;

dMdx = mach_term * ( friction_term - area_term );

% Catch singularity slipping thru
%if( dMdx > M_plus )
%    dMdx = M_plus;
%end
%
%if( dMdx < M_minus )
%    dMdx = M_minus;
%end
%display( [ ' // dMdx = ', num2str( dMdx ), ' // x = ', num2str(x) ] );

end

