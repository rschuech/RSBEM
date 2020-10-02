
clear guesses

for i = 1:size(un,1)
    
   [~,ind] = ismember(un(i,:),X_Y_unq,'rows');
   
   if ind ~=0
       guesses(i,:) = [best_amps(ind) best_lambdas(ind) best_nlambdas(ind)];
       continue
   end
   
   [~,ind] = ismember(un(i,:),bodies,'rows');
   if ind ~= 0
       guesses(i,:) = tails(ind,:);
       continue
   end
   
   disp(['do not have ',num2str(un(i,:))]);
   i
   
end