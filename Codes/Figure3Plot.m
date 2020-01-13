xNum  = 200;
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 1; % xMin = 0
tMax  = 1;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData);
u = 0.3*xMesh-0.75*tMesh+1;
figure
for tn = floor(linspace(1,tNum,10))
    plot(xData,u(tn,:),'LineWidth',1.5,'Color',[0 0.4470 0.7410,tn/tNum])
    hold on 
end

xlabel('$x$','Interpreter','Latex')
ylabel('$u$','Interpreter','Latex')
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)
%%
xNum  = 200;
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 1; % xMin = 0
tMax  = 1;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData);
u = 1+exp(-0.95*tMesh-xMesh);
figure
for tn = floor(linspace(1,tNum,10))
    plot(xData,u(tn,:),'LineWidth',1.5,'Color',[0 0.4470 0.7410,tn/tNum])
    hold on 
end

xlabel('$x$','Interpreter','Latex')
ylabel('$u$','Interpreter','Latex')
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)
%%
xNum  = 200;
% xNum  = 301;
tNum  = floor(xNum^(7/8));
xMax  = 10; % xMin = 0
tMax  = 1;
dx    = xMax/(xNum-1);
dt    = tMax/(tNum-1);
xData = 0:dx:xMax;
tData = 0:dt:tMax;
[xMesh,tMesh] = meshgrid(xData,tData);
u = 0.3*xMesh./(1-0.3*tMesh)+1-0.3*tMesh;
figure
for tn = floor(linspace(1,tNum,10))
    plot(xData,u(tn,:),'LineWidth',1.5,'Color',[0 0.4470 0.7410,tn/tNum])
    hold on 
end

xlabel('$x$','Interpreter','Latex')
ylabel('$u$','Interpreter','Latex')
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)