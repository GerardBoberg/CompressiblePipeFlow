function [ x_out, Mach_out, T, T0, P, P0 ] = AnalyticFannoFlow( span,...
                                                M_i, T_i, P_i,...
                                                ~, Dt, ~,...
                                                f, gamma )
%AnalyticFannoFlow Propogates a fanno flow
%   Detailed explanation goes here

%% Check inputs
L    = abs( span(end) - span(1) );

if( f < 1e-10 )
    error( 'Fanno flow needs a positive friction factor!' )
end

%% Solve for Mach number and locations
% Solve for L* from our initial condition. Use that to propogate distance
%   into Mach number.
% 4f L* / D = (1-M^2)/gammaM^2 + (gamma+1/2*gamma) *
%              ln( (gamma+1/2)M^2 / (1 + gamma-1/2 M^2 ) )
% L* = D/4f * (...)
fL4D      = @( M )( ( (1 - M.^2) ./ (gamma*M.^2) ) +...
                (gamma+1) / (2*gamma) .*...
                log( ((gamma+1)/2) .* M.^2 ./ (1+((gamma-1)/2) .* M.^2) ) );
            
fL4D_star = fL4D( M_i );
L_star_i = fL4D_star * ( Dt / ( 4 * f ) );

%% Fanno flow mach propogations
% Fanno flow is choked. Find L values for locations up to sonic point, and
%   10 points for the flat Mach 10 past it.
if( L_star_i <= L )
   Mach_out = linspace( M_i, 1, 100 );
   x_out    = span(1) + L_star_i - ( (Dt/(4*f)) .* fL4D( Mach_out ) );
   
   x_out    = [ x_out(1:end-1), linspace( x_out(end), span(end), 10 ) ];
   Mach_out = [ Mach_out(1:end-1), linspace( Mach_out(end), 1, 10 ) ];
else
    error( 'TODO: implement non-choked fanno' );
end

%% Thermodynamic properties

% Fanno pressure relations
% P2 / P1 = (M1/M2) sqrt( (2+(gamma-1)M1^2) / (2+(gamma-1)M262) )
P = P_i * (M_i ./ Mach_out) .*...
                sqrt( (2+(gamma-1)*M_i^2) ./ (2+(gamma-1)*Mach_out.^2) );

% Fanno Temperature relations
% T2 / T1 = (1 + (g-1/2)M_1 ^2 ) / (1 + (g-1/2)M_2^2 )
g = (gamma - 1) / 2;
T = T_i * 1 * ( 1 + g * M_i^2 ) ./ ( 1 + g * Mach_out.^2 );

% For no heat transfer, T0_1 = T0_n
% Using isentropic relations,
% T0 / T = 1 + ( (gamma - 1) / 2 ) M^2
% T0 = T * ( T0 / T )
T0_i =  T_i * ( 1 + ( (gamma - 1) / 2 )* M_i^2 );
T0   = T0_i .* ones( size( T ) );

% P / P0 = ( T / T0 )^(gamma/(gamma-1)) 
% P0 = P * 1/(P/P0)
p_exp = gamma / ( gamma - 1 );
P0    = P .* ( 1 ./ ( ( T ./ T0 ).^( p_exp ) ) );

% Normalize Thermo properties by inlet conditions
T  = T  ./  T(1);
P  = P  ./  P(1);
T0 = T0 ./ T0(1);
P0 = P0 ./ P0(1);
end

