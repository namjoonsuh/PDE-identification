N = logspace(1,8,1000);
K = 10;
Pmax = 2;
unorm = 1;
% sigma = 0.3;
figure
h(1) = line([10,10^8],[1,1],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',1.5);
hold on
pltNo = 2;
for sigma = [0.1,0.3,0.5]
    P = (8*K+2)*N.^((13+Pmax)/7).*K.*exp(-(N.^(1/7)-unorm).^2/(2*sigma^2));
    h(pltNo) = plot(N,P,'LineWidth',1.5);
    pltNo = pltNo + 1;
    hold on
end
legend([h(2),h(3),h(4)],{'$\sigma=0.1$','$\sigma=0.3$','$\sigma=0.5$'},'Interpreter','Latex')
xticks([10,10^4,10^8])
xticklabels({'$10$','$10^4$','$10^8$'})
yticks([10^-20,1,10^20])
yticklabels({'$10^{-20}$','$1$','$10^{20}$'})
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('$N$','Interpreter','Latex')
ylim([10^-20,10^20])
xlim([10,10^8])
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)
%%

N = logspace(1,8,1000);
K = 10;
Pmax = 2;
% unorm = 1;
sigma = 0.3;
figure
h(1) = line([10,10^8],[1,1],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',1.5);
hold on
pltNo = 2;
for unorm = [1,3,5]
    P = (8*K+2)*N.^((13+Pmax)/7).*K.*exp(-(N.^(1/7)-unorm).^2/(2*sigma^2));
    h(pltNo) = plot(N,P,'LineWidth',1.5);
    pltNo = pltNo + 1;
    hold on
end
legend([h(2),h(3),h(4)],{'$||u||_{L^{\infty}(\Omega)}=1$',...
    '$||u||_{L^{\infty}(\Omega)}=3$',...
    '$||u||_{L^{\infty}(\Omega)}=5$'},'Interpreter','Latex')
xticks([10,10^4,10^8])
xticklabels({'$10$','$10^4$','$10^8$'})
yticks([10^-20,1,10^20])
yticklabels({'$10^{-20}$','$1$','$10^{20}$'})
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('$N$','Interpreter','Latex')
ylim([10^-20,10^20])
xlim([10,10^8])
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)

%%
N = logspace(1,8,1000);
% K = 10;
Pmax = 2;
unorm = 1;
sigma = 0.3;
figure
h(1) = line([10,10^8],[1,1],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',1.5);
hold on
pltNo = 2;
for K = [10,50,100]
    P = (8*K+2)*N.^((13+Pmax)/7).*K.*exp(-(N.^(1/7)-unorm).^2/(2*sigma^2));
    h(pltNo) = plot(N,P,'LineWidth',1.5);
    pltNo = pltNo + 1;
    hold on
end
legend([h(2),h(3),h(4)],{'$K=10$',...
    '$K=50$',...
    '$K=100$'},'Interpreter','Latex')
xticks([10,10^4,10^8])
xticklabels({'$10$','$10^4$','$10^8$'})
yticks([10^-20,1,10^20])
yticklabels({'$10^{-20}$','$1$','$10^{20}$'})
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('$N$','Interpreter','Latex')
ylim([10^-20,10^20])
xlim([10,10^8])
set(gca,'TickLabelInterpreter', 'Latex');
set(gca,'fontsize',30)