function [xh,tt] = SGD(x_0, M, y, T, h)
    tt = 0:h:T;
    tt = tt';
    xh = zeros(2,length(tt));
    xh(:,1) = x_0;
    for i = 2:length(tt)
        R = M*xh(:,i-1)-y;
        j = randi([1,size(M,1)]);
        xh(:,i) = xh(:,i-1) - h * R(j) * M(j,:)';
    end
end