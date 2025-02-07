warning off;

n = 5;
p = 2;
x_0 = [0 0]';
rng(439)
M = randn([n,p]);
y = null(M');
y = y(:,2);
L = 8;
T = 10;
hvec = 2.^-[1,2,3,4];
K = floor(T./hvec);

f = @(x) 1/(2*n)*norm(M*x - y)^2;
gradf = @(x) 1/n * M' * (M*x - y);

j = 1;
e = 1;
figure;
for h=hvec
    disp(h)
    [xxGD,tt] = GD(x_0,gradf,T,h);
    
    yyGD=zeros(size(tt));
    
    for i=1:length(tt)
        yyGD(i)=f(xxGD(:,i));
    end
    subplot(4,3,j)
    yySGD = zeros(1,length(tt));
    yyEM = zeros(1,length(tt));
    hold on 
    for l=1:L
        [xxSGD,tt] = SGD(x_0,M,y,T,h);
        plot(xxSGD(1,:),xxSGD(2,:),'r-o')
        for i=1:length(tt)
            yySGD(i)=yySGD(i)+f(xxSGD(:,i));
        end
        [xxEM,ttEM] = EM(x_0,M,y,T,h);
        plot(xxEM(1,1:2^(10-e):length(ttEM)),xxEM(2,1:2^(10-e):length(ttEM)),'g-o')
        k = 1;
        for i=1:2^(10-e):length(ttEM)
            yyEM(k)=yyEM(k)+f(xxEM(:,i));
            k=k+1;
        end
    end
    yySGD = yySGD/L;
    yyEM = yyEM/L;
    plot(xxGD(1,:),xxGD(2,:),'b-o')
    xlabel('$x_1$', 'Interpreter','latex','FontSize',14)
    ylabel('$x_2$', 'Interpreter','latex','FontSize',14)
    xmin = (M'*M)\(M'*y);
    if e==1
        x1 = linspace(xmin(1)-1,xmin(1)+0.6,100);
        x2 = linspace(xmin(2)-1.2,xmin(2)+0.8,100);
    else
        x1 = linspace(xmin(1)-0.6,xmin(1)+0.5,100);
        x2 = linspace(xmin(2)-0.8,xmin(2)+0.5,100);
    end
    [X1,X2] = meshgrid(x1,x2);
    Z = arrayfun(@(x1,x2) f([x1;x2]),X1,X2);
    levels = linspace(0,0.5,20);
    contour(X1,X2,Z,levels)
    hold off
    j = j + 1;
    subplot(4,3,j)
    hold on 
    plot(tt,yySGD,'r-*')
    plot(tt,yyEM,'g-*')
    plot(tt,yyGD,'b-*')
    ylim([0.1,0.2])
    hold off
    title(['$h= 2^{-',num2str(e),'}$'], 'Interpreter','latex','FontSize',15, 'HorizontalAlignment','center')
    xlabel('$t$', 'Interpreter','latex','FontSize',14)
    ylabel('$f(x)$', 'Interpreter','latex','FontSize',14)
    if j~=2
        lgd = legend('SGD','SDE','GD','Location','northeast','FontSize', 6);
    else
        lgd = legend('SGD','SDE','GD','Location','southeast','FontSize', 6);
    end
    lgd.ItemTokenSize = [10,10];
    j = j + 1;
    N=10000;
    X=zeros(2,N);
    for k=1:N
        [x,tt] = EM(x_0,M,y,T,h);
        X(:,k) = x(:,end);
    end
    subplot(4,3,j)
    
    
    nbins = [20,20]; 
    [counts, centers] = hist3(X', 'Nbins', nbins);
    
    imagesc(centers{1}, centers{2}, counts'); 
    xlabel('$x_1$','Interpreter','latex','FontSize',14);
    ylabel('$x_2$','Interpreter','latex','FontSize',14);
    
    sigma2 = h/2 * f([0, 0]');
    [xGrid, yGrid] = meshgrid(linspace(min(X(1,:)), max(X(1,:)), 100), linspace(min(X(2,:)), max(X(2,:)), 100));

    gaussian_density = exp(-0.5 * ((xGrid.^2 + yGrid.^2) / sigma2)) / (2 * pi * sigma2);
    
    hold on;
    contour(xGrid, yGrid, gaussian_density, [0.1,1.5,10,30,47], 'LineColor', 'r'); 
    hold off;
    
    if e==1
        xlim([-0.88,0.88])
        ylim([-1.1,1.1])
    end 

    if e==2
        xlim([-0.58,0.58])
        ylim([-0.7,0.7])
    end 

    if e==3
        xlim([-0.5,0.5])
        ylim([-0.4,0.4])
    end 

    if e==4
        xlim([-0.26,0.26])
        ylim([-0.3,0.3])
    end 

    j = j + 1;
    e = e + 1;
end

set(gcf,'Units','centimeters','Position',[2,2,19,22]);
exportgraphics(gcf,'q62.png','Resolution',600)