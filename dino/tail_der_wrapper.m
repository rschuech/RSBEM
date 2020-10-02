function integrand = tail_der_wrapper(t, time, geom)

if size(t,2) > 1
    t = t';
end

[~,~,der] = tail_parameterized(t, time, geom);

 integrand = sqrt(sum( der .^2 , 2) );

 