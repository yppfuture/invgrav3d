function func = meshfunc(curve, weight, varagin) %#ok<INUSD>
%MESHFUNC Returns function of set parmaters for building graded meshes.
%   Meshfunc takes a set of parameters and an optional interactive mode
%   flag and returns a function for building graded meshes. The returned
%   function takes a regularly spaced vector x and returns a two vectors.
%
%   gmesh = meshfunc(curve, wieght, interactive)
%   [x, dx] = gmesh( linspace(0, 1, n) )
%
%   x samples values on the graded mesh and dx is equal to diff(x)
%
%   The Inputs to meshfunc are curve - which adjusts the steepness in the
%   change of dx spacing. Weight, which adjusts the ratio between the
%   larger dx spacing at the boundaries and the smallest dx spacing in the
%   centre dx(centre) / dx(boundary).
%   Meshfunc also takes the optional parameter interactive, which should
%   set to either true, or false or left out. Default is false.
%   IF set to true, it allows user to choose a c and w paramter
%   dynamically. When the user presses the Accept button a return func
%   will be created using the parameters from the last iteration, the
%   parameters shown on the plot and as the default values in the dialogue
%   box.


%   Ben Postlethwaite 2012
%   benpostlethwaite.ca

c = curve;
w = weight;
n = 40;
interactive = false;

if nargin == 3
    interactive = varagin;
end

if interactive
    out = meshgui();
    c = out(1);
    w = out(2);
end

% Gmesh builds a fat parabola with minimum dictated by wieght
% The w/(1-w) is necessary to that the scaling is preserved.
gmesh = @(x) ((x - (x(end) - x(1))/2 - x(1)).^c + (w / (1-w) ) * ((x(end) - x(1))/2).^c);
% getdx turns gmesh into an x and dx
getdx = @(x) (x(end) - x(1)) * gmesh(x) / sum(gmesh(x));
%getdx = @(x) gmesh(x);
getx = @(x) cumsum([x(1), getdx(x)]);
% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);
% Get dims performs a few function calls to make things easier.
func = @(x) deal( getx(centre(x)), getdx(centre(x)) );

end



