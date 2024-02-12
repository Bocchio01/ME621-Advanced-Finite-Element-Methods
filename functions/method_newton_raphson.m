function U = method_newton_raphson(U, R, Kt_inv)

if (isscalar(U))
    assert(isscalar(R) && isscalar(Kt_inv), ...
        'U, delta_P, and Kt_inv must be scalars or vector.');
elseif (isvector(U))
    assert(isvector(R) && ismatrix(Kt_inv) && ...
        size(U, 1) == size(Kt_inv, 1) && ...
        size(U, 2) == size(R, 2) && ...
        size(Kt_inv, 2) == size(R, 1), ...
        'Invalid dimensions for matrix multiplication.');
else
    assert(false, 'Dimensions of input not supported yet')
end

U = U - Kt_inv * R;

end

