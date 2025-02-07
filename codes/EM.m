function [xh,tt] = EM(x_0, M, y, T, h)
    n = length(y);
    dt = 2^-10;
    tt = 0:dt:T;
    tt = tt';
    xh = zeros(2,length(tt));
    xh(:,1) = x_0;
    for i = 2:length(tt)
        R = M*xh(:,i-1) - y;
        sigma = 1/n * M' * (diag(R.^2) - 1/n * (R * R')) * M;
        [U,D,V] = svd(sigma);
        s12 = U * sqrt(D) * V';
        xh(:,i) = xh(:,i-1) - dt * 1/n*M'*(M*xh(:,i-1)-y) + sqrt(dt) * sqrt(h) * s12 * randn([2,1]);
    end
end