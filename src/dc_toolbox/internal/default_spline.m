function paramValuesSpline = default_spline(paramValues,knots)
assert(length(knots)>1,'need at least 4 splines for default spline function')

% we add 3 knots (because cubic, 4th order splines) in the
% beginning and the end 
knots = [repmat(knots(1),1,3) knots repmat(knots(end),1,3)];


% This functino always removes either first or last spline. We therefore
% need to recover it by running it twice and concatenating
a = Bernstein(paramValues,knots,[],4,[],0);
b = Bernstein(paramValues,knots,[],4,[],1);

paramValuesSpline = a;
paramValuesSpline(b(:)==1) = 1;
