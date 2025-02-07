n = 2;
p = 2;
x_0 = [0 0]';
rng(876)
M = randn([2,2]);
y = randn([2,1]);
L = 8;
T = 10;
hvec = 2.^-[1,2,3,4];
tau = 0.01;
f = @(x) 1/(2*n)*norm(M*x - y)^2;
gradf = @(x) 1/n * M' * (M*x - y);

j = 1;
e = 1;

figure;
for h=hvec
    
    [xxGD,tt] = GD(x_0,gradf,T,h);
    
    yyGD=zeros(size(tt));
    
    for i=1:length(tt)
        yyGD(i)=f(xxGD(:,i));
    end
    subplot(4,2,j)
    yyNGD = zeros(1,length(tt));
    yyEM = zeros(1,length(tt));
    hold on 
    for l=1:L
        [xxNGD,tt] = NGD(x_0,gradf,T,h,tau);
        plot(xxNGD(1,:),xxNGD(2,:),'r-o')
        for i=1:length(tt)
            yyNGD(i)=yyNGD(i)+f(xxNGD(:,i));
        end
        [xxEM,ttEM] = EM2(x_0,M,y,T,2^-8,tau);
        plot(xxEM(1,1:2^(10-e):length(ttEM)),xxEM(2,1:2^(10-e):length(ttEM)),'g-o')
        k = 1;
        for i=1:2^(10-e):length(ttEM)
            yyEM(k)=yyEM(k)+f(xxEM(:,i));
            k=k+1;
        end
    end
    yyNGD = yyNGD/L;
    yyEM = yyEM/L;
    plot(xxGD(1,:),xxGD(2,:),'b-o')
    xlabel('$x_1$', 'Interpreter','latex','FontSize',14)
    ylabel('$x_2$', 'Interpreter','latex','FontSize',14)
    xmin = M\y;
    x1 = linspace(xmin(1)-0.5,xmin(1)+0.5,100);
    x2 = linspace(xmin(2)-0.8,xmin(2)+0.5,100);
    [X1,X2] = meshgrid(x1,x2);
    Z = arrayfun(@(x1,x2) f([x1;x2]),X1,X2);
    levels = linspace(0,0.5,20);
    contour(X1,X2,Z,levels)
    hold off
    j = j + 1;
    subplot(4,2,j)
    hold on 
    plot(tt,yyNGD,'r-*')
    plot(tt,yyEM,'g-*')
    plot(tt,yyGD,'b-*')
    hold off
    title(['$h= 2^{-',num2str(e),'}$'], 'Interpreter','latex','FontSize',15,'Position',[4.5,0.06,0], 'HorizontalAlignment','center')
    xlabel('$t$', 'Interpreter','latex','FontSize',14)
    ylabel('$f(x)$', 'Interpreter','latex','FontSize',14)
    legend('NGD','SDE','GD',location='northeast')
    j = j + 1;
    e = e + 1;
    k = k + 1;
end

set(gcf,'Units','centimeters','Position',[2,2,19,22]);
exportgraphics(gcf,'q84.png','Resolution',600)