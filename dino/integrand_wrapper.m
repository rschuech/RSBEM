function integrand = integrand_wrapper(t, time, geom)

integrand = NaN(size(t));

for i = 1:numel(t)
    
[~,~,der] = tail_parameterized(t(i), time, geom);

 integrand(i) = sqrt(sum( der .^2 , 2) );

end
 
 