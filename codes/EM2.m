function [xh,tt] = EM2(x_0, M, y, T, h, tau)
    n = length(y);
    tt = 0:h:T;
    tt = tt';
    xh = zeros(2,length(tt));
    xh(:,1) = x_0;
    for i = 2:length(tt)
        xh(:,i) = xh(:,i-1) - h * 1/n *M' * (M*xh(:,i-1)-y) + sqrt(2*tau*h) * randn(size(x_0));
    end
end