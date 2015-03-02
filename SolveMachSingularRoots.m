function [ M_plus, M_minus ] = SolveMachSingularRoots( x_choke, Dfun, Dt,...
                                                         L, f, gamma )
%SolveMachSingularRoots Finds the two differential roots at the soinc point
%
%   Specifically, the Mach equation for adiabatic flow is singular at the 
%     point where M = 1. Via LoHop's rule, the derivative at M = 1 is the
%     following Quadratic equation:
%
%               4                gamma f        d( ( dAdx / A ) M->1 )
%   0   = ------------- X^2 + ------------- X + ---------------------- +...
%          (gamma + 1)         Dfun(M->1)                 dx
%
%           gamma f     d( 1 / D(M->1) )
%         ------------  ----------------
%              2               dx
%
%   Where X is the limit of  (dM/dx) M->1.
%   This produces two roots, a positive M_plus and a negative M_minus.
%

%% Check inputs
% Make assumptions about inputs we don't know.

if( nargin < 4 )
    gamma = 1.4; % gamma for air
end

if( nargin < 3 )
    f = 0; % isentropic, frictionless flow
end

%% Solve for quadratic coefficients
% The equation takes the form of 
%        a x^2 + b x + c = 0
% Setup the numbers into that form

% a and b are easy
a = 4 / ( gamma + 1 ); 

b = gamma * f / Dfun( x_choke, Dt, L );

% c is hard. We need to do finite-differences for each derivative.
delta = 1e-4; % .0001

% First, find d/dx ( 1 / D) 
range = [ x_choke - delta, x_choke + delta ];
dDdx  = diff( 1 ./ Dfun( range, Dt, L ) ) ./ ( 2 * delta );

% when f == 0, dDdx will also be zero at the choke point.
if( abs( dDdx ) <= 1e-6 )
    dDdx = 1; 
    f    = 0;
end


% We need to find the:  derivative of { (1/A) * derivative of A  }
% Where A is a Lambda funcction for Area = pi/4 D^2
A = @( x, Dt, L )( (pi/4) .* ( Dfun( x, Dt, L ).^2 ) );

% So we need a Lambda function for the derivative of Area, normalized by A:
dAdx = @( x, Dt, L)...
    ( diff( A( [x-delta, x+delta], Dt, L ) ) ./...
                                      ( 2 * delta * A( x, Dt, L) ) );
  
% And take a finite difference of that derivative lambda function. 
d2dx2 = ( dAdx( x_choke + delta, Dt, L )...
        - dAdx( x_choke - delta, Dt, L ) )...
        ./ ( 2*delta );
% dAdx     = diff( A( range, Dt, L ) ) ./ ( 2 * delta );
% dAdx     = dAdx / A( choke, Dt, L );
%ddAdxAdx = dAdx;

% That mess is done. Put it into c
area_term     = -1 * d2dx2;
friction_term = gamma * f * 0.5 / dDdx;

c = area_term + friction_term;


%% Determine which roots is M_plus and M_minus 
% Should produce two roots, one above and one below zero.
% solves for x in ( a x^2 + b x + c = 0 )

r = roots( [ a, b, c] );

if( r(1) > r(2) )
    M_plus  = r(1);
    M_minus = r(2);
else
    M_plus = r(2);
    M_minus = r(1);
end

