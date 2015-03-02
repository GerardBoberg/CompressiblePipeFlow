function [ x_out, Mach ] = PropogateMachNumber( span, init, Dfun, Dt,...
                                                  L, f, gamma,...
                                                  ss_ac, x_choke,...
                                                  M_plus, M_minus )
%PropogateMachNumber Uses ODE45 to propogate the differential Mach equation
%
%   Specifically, the Mach equation for adiabatic flow is:
%
%   dM     2M + (gamma - 1)M^3    {  gamma M^2 f       dA/dx    }
%   --  = --------------------- * { ------------- - ----------- }
%   dx       2 ( 1 - M^2 )        {       2 D            A      }
%
%   Which is singular as M approaches 1
%
%   This function propogates that equation over the interval x.

%% Run the numeric approximation -- ode45
dM      = @DifferentialMach;
x_span  = span;
M_0     = init;
options = odeset( 'RelTol', 1e-4, 'AbsTol', 1e-5, 'Stats', 1 );
%Dfun    = Dfun;
%Dt      = Dt;
%L       = L;
%f       = f;
%gamma   = gamma;
% DifferentialMach( x, M, Dfun, Dt, L, f, gamma )
[ x_out, Mach ] = ode45( dM, x_span, M_0, options, Dfun, Dt, L, f, gamma,...
                                            ss_ac, x_choke, M_plus, M_minus );


end

