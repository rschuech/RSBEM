
% finding LCM of transverse and tail frequencies
  % UNFINISHED......

% When x and y are rational numbers, here's one way to do it:
% 
% 1. Write x and y with a common denominator; say x = a/d and y = b/d
% 2. Find the least common multiple of a and b by using the method for integers; let us say c = lcm(a, b)
% 3. The least common multiple of x and y is c/d. 


transverse_omega = 45.4;  transverse_period = roundn(2 * pi / transverse_omega, -2)
tail_omega = 46;          tail_period = roundn(2 * pi / tail_omega, -2)



%%


d = 1E6 / pi;  %just make sure you get integers after multiplying both omegas by d (so make sure d is big enough)

a = transverse_period * d
b = tail_period * d

c = lcm(a,b);

lcm_omega = c / d

