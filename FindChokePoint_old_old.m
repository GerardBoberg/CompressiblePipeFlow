function [ x_choke, index_choke, error ] = FindChokePoint( x, D, f, gamma )
%FindChokePoint Finds the choke point of a converging-diverging nozzle
%   Uses MATLAB's built-in derivative estimator, diff, as the input to the
%     governing equation for the sonic point:
%
%      ( (1/A)(dA / dx) )   = 2 * gamma * ( f / D )
%                    M -> 1                     M -> 1
%
%   x     --- Array of locations where the pipe's geometry is calculated
%   D     --- Array of diameter values at every x location
%   f     --- Friction factor of fanno flow.
%               Defaults to 0, an edge case for no friction.
%   gamma --- The gamma value of the gas, ratio of specific heats.
%               Defaults to 1.4, the gamma for dry air.


%% Check inputs to the function

% User didn't give us a gamma value for the gas.
if( nargin <= 3 )
    gamma = 1.4; % gamma for air = 1.4. Use as default value
end

% Edge case! Frictionless isentropic flow
if( (nargin <= 2) || ( f == 0 ) )
    % the choke will always be at the throat -- the min diameter
    [ ~, i ] = min( abs( D ) );  
    
    x_choke     = x( i );
    index_choke =    i;
    return; % no further calcs needed for this edge case.
end


%% Calculate the sonic point

A    = (pi/4) * D.^2;     % Area, m^2
dAdx = diff( A ) * 100;   % diff returns percentage derivative
dAdx = [ dAdx(1), dAdx ]; %   that is 1 spot shorter b/c can't estimate
                          %   derivative of A(1) due to no A(0) info.
                              

                          
% The actual equation being calculated, solved for zero.
error = ( dAdx ./ A ) - ( 2 * gamma * f ./ D );


plot( x, dAdx ./ A, x, (2*f*gamma ./ D), x, error );
legend( 'dAdx/A', 'f/D', 'err' );

% Now solve find the zero location.
%    There are actually two (or more?) zeros, and we only want the first
%    one that cross from negative values to positive ones.

index_choke = 1; % pre-allocate, so survive's loop's scope
se = sign( error ); % array of [ -1, -1, -1, 0, 1, 1, 1]

for i = 2:length(error)       % search the sign array
    if( se( i-1 ) < se( i ) ) % for the spot where -1 turns to +1
        
        % We've found the sonic point. It's between the + and - roots.
        %   figure out which one is closer to zero:
        difference = abs( error( i-1 ) ) - abs( error( i ) );
        if( difference >=0 )
            index_choke = i; % (i-1) has a larger error. Use i.
        else %if( difference < 0 )
            index_choke = i-1; % i has a larger error. Use (i-1)
        end
    end
end

% And lookup where that x value actually is
x_choke     = x( index_choke );

end

