function [xh,tt] = NGD(x_0, gradf, T, h, tau)
    tt = 0:h:T;
    tt = tt';
    xh = zeros(2,length(tt));
    xh(:,1) = x_0;
    for i = 2:length(tt)
        xh(:,i) = xh(:,i-1) - h * gradf(xh(:,i-1)) + sqrt(2*h*tau)*randn(size(x_0));
    end
end