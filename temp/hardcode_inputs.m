
Inputs.body.AR = [5.75 0.65]';
Inputs.tail.amp = 0.634;  Inputs.tail.lambda = 4.32;  Inputs.tail.nlambda = 1.345;


i = 1;

Inputs(i).paths.namebase.body = [Inputs(i).body.shape,'_AR1_',num2str(Inputs(i).body.AR(1)),'_AR2_',num2str(Inputs(i).body.AR(2)),Inputs(i).body.suffix];

Inputs(i).paths.namebase.tail = ['tail_','radius_',num2str(Inputs(i).tail.radius,digits),'_amp_',num2str(Inputs(i).tail.amp,digits),'_lambda_',num2str(Inputs(i).tail.lambda,digits),'_nlambda_',num2str(Inputs(i).tail.nlambda,digits),Inputs(i).tail.suffix];
rest = Inputs(i).paths.suffix;

Inputs(i).paths.namebase.full = [Inputs(i).paths.namebase.body,'_',Inputs(i).paths.namebase.tail,'_motorBC_',Inputs(i).tail.motorBC, rest];
