function [ T, P, T0, P0 ] = MethodOfBeansThermo( x, Mach, Afun,...
                                                 T_i, P_i, gamma )
%MethodOfBeansThermo Propogates Thermodynamic properties given Mach array.
%   
%   x    --- The array of x-locations that corrospond to the Mach array.
%   Mach --- The array of Mach numbers to propogate thermo properties over.
%   Afun --- A function handle to the area of the pipe at a given x. F(x).
%   T_i  --- The static temperature at x(1), Mach(1)
%   P_i  --- The static pressure at x(1), Mach(1)
%  gamma --- Ratio of specific heats of the gas. Use 1.4 for air.
%
% Outputs the absolute temperature and pressure, static and stagnation at
%   every ( x, M ) point. Recomend normalizing outputs for graphs.

%% Check inputs
if( nargin < 5 )
    gamma = 1.4; % gamma for air
end


%% Step 1: Calculate Temperature variation
%
% T2 / T1 = (T01/T02) * (1 + (g-1/2)M_1 ^2 ) / (1 + (g-1/2)M_2^2 )
% Let T1 = const = T_i, M_1 = M(1)
% For no heat transfer, T01 / T02 = 1 ==> T0 is costant
%
% T0 / T = 1 + ( (gamma - 1) / 2 ) M^2
% T0 = T * ( T0 / T )
T0_i =  T_i * ( 1 + ( (gamma - 1)/2 )*Mach(1)^2 );
T0   = T0_i * ones( size( Mach ) );

g = (gamma - 1) / 2;
T = T_i * 1 * ( 1 + g * Mach.^2 ) ./ ( 1 + g * Mach(1)^2 );

%% Step 2: Calculate Pressure Variation
%
% P2 / P1 = (A1 / A2) * ( M1 / M2) * sqrt( T02 / T01 ) *...
%                              sqrt( (1+gM1^2)/(1+gM2^2) )
% For no heat transfer, T01 / T02 = 1

P = ( Afun( x(1) ) ./ Afun( x ) );
P = P .* ( Mach( 1 ) ./ Mach );
P = P .* sqrt( ( 1 + g * Mach( 1 )^2 ) ./ ( 1 + g .* Mach.^2 ) );
P = P .* P_i;

% Now to find out P stagnation:
%
% T0 / T = 1 + ( (gamma - 1) / 2 ) M^2
% T0 = T * ( T0 / T )
T0overT = T0 ./ T;

% P / P0 = ( T / T0 )^(gamma/(gamma-1))
% P0 = P * 1/(P/P0)
p_exp      = gamma / ( gamma - 1 );
PoverP0    = ( 1 ./ T0overT ).^( p_exp );
P0         = P .* ( 1 ./ PoverP0 );


end

