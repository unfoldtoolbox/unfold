function dmt = cyclical_spline(x, knots,bounds)
%     """Builds an unconstrained cubic regression spline design matrix.
%     Returns design matrix with dimensions ``len(x) x n``
%     for a cubic regression spline smoother
%     where
%      - ``n = len(knots)`` for natural CRS
%      - ``n = len(knots) - 1`` for cyclic CRS
%     .. note:: See 'Generalized Additive Models', Simon N. Wood, 2006, p. 145
%     :param x: The 1-d array values.
%     :param knots: The 1-d array knots used for cubic spline parametrization,
%      must be sorted in ascending order.
%     :param cyclic: Indicates whether used cubic regression splines should
%      be cyclic or not. Default is ``False``.
%     :return: The (2-d array) design matrix.
%     """

cyclic = 1;
n = length(knots);
if nargin == 2 || isempty(bounds)
    bounds = [min(knots),max(knots)];
end
if cyclic
    x = map_cyclic(x, bounds(1), bounds(2));
    n = n- 1;
end

[ajm, ajp, cjm, cjp, j] = compute_base_functions(x, knots);
j1 = j+1;
if cyclic
    j1(j1 == n+1) = 1;
end

i = eye(n);
if cyclic
    f = get_cyclic_f(knots);
else
    f = get_natural_f(knots);
end

% This works only in R2016b, but is easier to read
% dmt = ajm .* i(j, :)' + ajp .* i(j1, :)'+ cjm .* f(j, :)'+ cjp .* f(j1, :)';

% reimplemented for 2014b
mult = @(x,y)bsxfun(@times,double(x),double(y));
dmt = mult(ajm,i(j,:)') + mult(ajp,i(j1,:)') + mult(cjm,f(j,:)') + mult(cjp,f(j1,:)');
dmt = dmt';




end

function [ajm, ajp, cjm, cjp, j] = compute_base_functions(x, knots)
%     """Computes base functions used for building cubic splines basis.
%     .. note:: See 'Generalized Additive Models', Simon N. Wood, 2006, p. 146
%       and for the special treatment of ``x`` values outside ``knots`` range
%       see 'mgcv' source code, file 'mgcv.c', function 'crspl()', l.249
%     :param x: The 1-d array values for which base functions should be computed.
%     :param knots: The 1-d array knots used for cubic spline parametrization,
%      must be sorted in ascending order.
%     :return: 4 arrays corresponding to the 4 base functions ajm, ajp, cjm, cjp
%      + the 1-d array of knots lower bounds indices corresponding to
%      the given ``x`` values.
%     """
j = find_knots_lower_bounds(x, knots);

h = knots(2:end) - knots(1:end-1);
hj = h(j);
xj1_x = knots(j+1) - x;
x_xj = x - knots(j);

ajm = xj1_x ./ hj;
ajp = x_xj ./ hj;

cjm_3 = xj1_x .* xj1_x .* xj1_x ./ (6. .* hj);
cjm_3(x > max(knots)) = 0.;
cjm_1 = hj .* xj1_x ./ 6;
cjm = cjm_3 - cjm_1;

cjp_3 = x_xj .* x_xj .* x_xj ./ (6. .* hj);
cjp_3(x < min(knots)) = 0;
cjp_1 = hj .* x_xj ./ 6;
cjp = cjp_3 - cjp_1;

end


function lb = find_knots_lower_bounds(x, knots)
%     """Finds knots lower bounds for given values.
%     Returns an array of indices ``I`` such that
%     ``0 <= I[i] <= knots.size - 2`` for all ``i``
%     and
%     ``knots[I[i]] < x[i] <= knots[I[i] + 1]`` if
%     ``np.min(knots) < x[i] <= np.max(knots)``,
%     ``I[i] = 0`` if ``x[i] <= np.min(knots)``
%     ``I[i] = knots.size - 2`` if ``np.max(knots) < x[i]``
%
%     :param x: The 1-d array values whose knots lower bounds are to be found.
%     :param knots: The 1-d array knots used for cubic spline parametrization,
%      must be sorted in ascending order.
%     :return: An array of knots lower bounds indices.
%     """
lb = [];
for xi = 1:length(x)
    tmp  =find(knots>x(xi),1)- 1;
    if isempty(tmp)
        tmp = length(knots);
    end
    lb(xi) = tmp; 
end

% I[i] = 0 for x[i] <= np.min(knots)
lb(lb == 0) = 1;

% I[i] = knots.size - 2 for x[i] > np.max(knots)
lb(lb == length(knots)) = length(knots) - 1;

end


function map = get_cyclic_f(knots)
%     """Returns mapping of cyclic cubic spline values to 2nd derivatives.
%     .. note:: See 'Generalized Additive Models', Simon N. Wood, 2006, pp 146-147
%     :param knots: The 1-d array knots used for cubic spline parametrization,
%      must be sorted in ascending order.
%     :return: A 2-d array mapping cyclic cubic spline values at
%      knots to second derivatives.
%     """
h = knots(2:end) - knots(1:end-1);
n = length(knots) - 1;
b = zeros(n);
d = zeros(n);

b(1, 1) = (h(n - 1) + h(1)) ./ 3.;
b(1, n) = h(n - 1) ./ 6.;
b(n, 1) = h(n - 1) ./ 6.;

d(1, 1) = -1. ./ h(1) - 1. ./ h(n - 1);
d(1, n) = 1. ./ h(n - 1);
d(n, 1) = 1. ./ h(n - 1);

for i = 2:n
    b(i, i) = (h(i - 1) + h(i)) ./ 3.;
    b(i, i - 1) = h(i - 1) ./ 6.;
    b(i - 1, i) = h(i - 1) ./ 6.;
    
    d(i, i) = -1. ./ h(i - 1) - 1. ./ h(i);
    d(i, i - 1) = 1. ./ h(i - 1);
    d(i - 1, i) = 1. ./ h(i - 1);
    
end
map = (b\d);

end


function x = map_cyclic(x, lbound, ubound)
%     """Maps values into the interval [lbound, ubound] in a cyclic fashion.
%     :param x: The 1-d array values to be mapped.
%     :param lbound: The lower bound of the interval.
%     :param ubound: The upper bound of the interval.
%     :return: A new 1-d array containing mapped x values.
%     :raise ValueError: if lbound >= ubound.
%     """
if lbound >= ubound
    error()
end
%         raise ValueError("Invalid argument: lbound (%r) should be "
%                          "less than ubound (%r)."
% (lbound, ubound))


x(x > ubound) = lbound + mod((x(x > ubound) - ubound), (ubound - lbound));
x(x < lbound) = ubound - mod((lbound - x(x < lbound)),(ubound - lbound));
end

