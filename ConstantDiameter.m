function [ diameter ] = ConstantDiameter( x, D, ~ )
%ConDiNozzleDiameter The diameter of a pipe. Used for Fanno Flow.
%
%   This extremely simple method is used with pure Fanno Flow where there
%     is no area change along the pipe. 
%
%   x   --- A vector of the x locations of the pipe diameter.
%   D   --- The diameter of the pipe.
%             Defaults to 0.15 m


% error check inputs
if( nargin < 2 )
    D = 3 * 0.05;  % default throat diameter
end

% It's a constant Diameter Pipe. What did you expect?
diameter = D * ones( size( x ) );

end

