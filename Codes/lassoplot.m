[LL,fitinfo] = lasso(noiseFmat,dUtVec,'Lambda',linspace(0,2.5,200));
figure
flags = zeros(1,size(LL,2));
for i = 1:size(LL,1)
    for j = 1:size(LL,2)
        if (i==4 || i==7)
            rectangle('Position',[fitinfo.Lambda(j),i,1,1],'FaceColor',(1-[0,1,1])*abs(LL(i,j))/(max(abs(LL(i,:)))+eps),'EdgeColor','none')
            hold on
        else
            rectangle('Position',[fitinfo.Lambda(j),i,1,1],'FaceColor',(1-[1,1,0])*abs(LL(i,j))/(max(abs(LL(i,:)))+eps),'EdgeColor','none')
        end
        
%         if ((abs(LL(i,j))/max(max(abs(LL)))<eps) && (flags(j)==0))
%             rectangle('Position',[fitinfo.Lambda(j),i,0.1,1],'FaceColor',[0,0,0],'EdgeColor','none')
%             flags(j) =1;
%         end
    end
end
xlim([min(fitinfo.Lambda),max(fitinfo.Lambda)])
xlabel('$\lambda$','Interpreter','Latex')
% xt = get(gca,'XTickLabel');
% xtt = flipud(xt);
% set(gca, 'XTickLabel', xtt)

yticks((1:10)+0.5)
yticklabels({'const','$u$','$u^2$','$u_x$','$u_x^2$','$uu_x$','$u_{xx}$','$u_{xx}^2$','$u_xu_{xx}$','$uu_{xx}$'})
set(gca,'TickLabelInterpreter','latex')
title('$\sigma=0.5, N=100$','Interpreter','Latex')
ylim([1,10])
set(gca,'fontsize',20)

set(gcf, 'Position',  [100, 100, 580, 300])

%%

figure
for i = 1:size(LL,1)
    if i == 4
        h1 = semilogx(fitinfo.Lambda,LL(i,:),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
    elseif i == 7
        h2 = semilogx(fitinfo.Lambda,LL(i,:),'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
    else    
        semilogx(fitinfo.Lambda,LL(i,:),'-.','Color',[0.5,0.5,0.5],'LineWidth',1.5)
    end
    hold on
end
xlim([0,2.5])
xlabel('$\lambda$','Interpreter','Latex')
hleg1 = legend(h1,{'$u_x$'},'Interpreter','Latex');
set(gca,'TickLabelInterpreter','latex')
title('$\sigma=0.03,M=1000$','Interpreter','Latex')

set(gca,'fontsize',20)