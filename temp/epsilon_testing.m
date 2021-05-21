

epsilons = sort([ 10.^(-9:1:0)  5*10.^(-9:1:0)]);

[Eps_body,Eps_tail] = ndgrid(epsilons,epsilons);

I = randperm(size(Eps_body,1));
J = randperm(size(Eps_body,2));
conds = NaN(size(Eps_body));

for l = randperm(numel(conds))
        
        eps_body = Eps_body(l);
        eps_tail = Eps_tail(l);
        [eps_body eps_tail];
    
    
    settings_inputs_blob_sweep;
    main_test2;
    
    conds(l) = condA;
    %%
    figure(3453);
   p =  pcolor(conds); shading flat; h = colorbar; drawnow;
   set(gca,'ColorScale','log');
    %%
end