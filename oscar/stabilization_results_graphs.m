% Expr 8

width = 0.87;
lengths = 1:0.5:10;
cyl_height = lengths - width;  %total length - 2*radius

AR1 = lengths / width;
AR2 = zeros(size(AR1));

load( 'C:\Hull\Results\Oscar\Oscar experiment 8.mat' );

speed = [Results.Avg_Speed];
tau = [Results.tau_a];

figure(8); clf;  set(gcf,'Position',[  -1390         164         857         705]);
yyaxis left
h1 = plot(lengths,speed,'o-','Color','r','MarkerFaceColor','r');  hold on;  ax1 = gca;  ax1.YColor = 'r';
xlabel('cell length (\mum)');  ylabel('swimming speed (\mum/s)');
yyaxis right
h2 = plot(lengths,tau,'o-','Color','b','MarkerFaceColor','b');  hold off;  ax2 = gca;  ax2.YColor = 'b';
ylabel('\tau for loss of orientation (s)');
grid on
ax1_pos = ax1.Position; % position of first axes
ax3 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','YColor','none');
ax3.XLim = [AR1(1) AR1(end)];
ax3.XTick = round( linspace(AR1(1),AR1(end),10) , 2 );
xlabel2 = annotation('textbox','Position',[0.454947368421051 0.965 0.157082706766917 0.030379746835443],...
    'String','cell aspect ratio','LineStyle','none','FontSize',11);

