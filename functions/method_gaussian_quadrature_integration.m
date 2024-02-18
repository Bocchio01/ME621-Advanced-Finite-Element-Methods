% roots = vpasolve(legendreP(7,x) == 0)

function result = method_gaussian_quadrature_integration(integrand_func, fem)
% METHOD_GAUSSIAN_QUADRATURE_INTEGRATION
%   This function computes the integral of a given function using the
%   Gaussian quadrature method.
%
%   result = METHOD_GAUSSIAN_QUADRATURE_INTEGRATION(integrand_func, fem)
%
%   integrand_func: function handle to the integrand
%   fem: struct containing the number of Gauss points required
%
%   result: value of the integral

switch fem.N_Gauss_points
    case 1
        weights = [2];
        roots = [0];
    case 2
        weights = [1, 1];
        roots = [-1/sqrt(3), 1/sqrt(3)];
    case 3
        weights = [5/9, 8/9, 5/9];
        roots = [-sqrt(3/5), 0, sqrt(3/5)];
    otherwise
        error('Order of Gaussian quadrature not implemented (too high)')
end

result = 0;
for root_idx = 1:length(roots)
    
    result = result + integrand_func(roots(root_idx)) * weights(root_idx);
    
end

end