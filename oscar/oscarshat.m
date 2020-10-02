amp = 0.2;
lambda = 2.34;
length = 6;

% lambda*t/(2*math.pi) = length  % tail 6 micron long along centerline

%tmax =  length / lambda * 2 * pi;

nlambda = length / lambda;


temp_input.performance.verbose = true;
temp_input.tail.amp = amp;
temp_input.tail.lambda = lambda;
temp_input.tail.nlambda = nlambda;

[mesh_succeed] = gen_mesh_wrapper(temp_input)



%% Oscar
new_sweep.AR1 = (1:12)';
new_sweep.AR2 = zeros(size(new_sweep.AR1));
new_sweep.amp = repmat(amp,size(new_sweep.AR1));
new_sweep.lambda = repmat(lambda,size(new_sweep.AR1));
new_sweep.nlambda = repmat(nlambda,size(new_sweep.AR1));

save('E:\Hull\Oscar_sweep.mat','new_sweep');

%%

width = 0.71;
length = [0.8051	0.9234	1.0592	1.2148	1.3934	1.5981	1.833	2.1024	2.4114	2.7658	3.1723	3.6385	4.1733	4.7866	5.4901	6.297	7.2224	8.2839	9.5014	10.8978	12.4994	14.3365	16.4435	18.8602	21.632	24.8113	28.4578	32.6402	37.4373	40];
cyl_height = length - width;  %total length - 2*radius

AR1 = length / width;
V = pi * (width/2)^2 * cyl_height + 4/3*pi*(width/2)^3;
