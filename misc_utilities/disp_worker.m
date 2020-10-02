function [] = disp_worker()

     task = getCurrentTask;
    
    ID = task.ID;
    disp(['Worker = ',num2str(ID),'      ','col_i = ',num2str(col_i)]);