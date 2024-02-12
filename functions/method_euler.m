function U = method_euler(U, delta_P, Kt_inv)

if (isscalar(U))
    assert(isscalar(delta_P) && isscalar(Kt_inv), ...
        'U, delta_P, and Kt_inv must be scalars or vector.');
elseif (isvector(U))
    assert(isvector(delta_P) && ismatrix(Kt_inv) && ...
        size(U, 1) == size(Kt_inv, 1) && ...
        size(U, 2) == size(delta_P, 2) && ...
        size(Kt_inv, 2) == size(delta_P, 1), ...
        'Invalid dimensions for matrix multiplication.');
else
    assert(false, 'Dimensions of input not supported yet')
end

U = U + Kt_inv * delta_P;

end
