n = 2; 
p = 2; 
T = 10; 
dt = 2^(-8);
x_0 = [0; 0];
rng(351)

M = randn([n, p]);
y = M * [1; 0];
tt = 0:dt:T;
hvec = 2.^-[1, 2, 3, 4];

R = @(x) M * x - y;
f = @(x) 1/(2*n) * norm(R(x))^2;
gradf = @(x) 1/n * M' * R(x);
hessf = @(x) 1/n * (M' * M);
sigma = @(x) 1/n * M' * (diag(R(x).^2) - (1/n) * R(x) * R(x)') * M;

err = zeros(1, length(hvec));

for k = 1:length(hvec)
    h = hvec(k);
    th = 0:h:T;
    phi1 = zeros(size(th));
    phi2 = zeros(size(th));

    for trials = 1:1000
        disp(trials)

        W = sqrt(dt) * randn(p, length(tt));
        W(:,1) = 0;
        W = cumsum(W, 2); 

        X = zeros(p, length(tt));
        X(:,1) = x_0;

        for i = 2:length(tt)
            grad = gradf(X(:,i-1));
            hess_grad = hessf(X(:,i-1)) * grad;
            s = sigma(X(:,i-1));
            [U, D, V] = svd(s);
            s12 = U * sqrt(D) * V';

            X(:,i) = X(:,i-1) ...
                     - dt * grad ...
                     - (h/2) * dt * hess_grad ...
                     + sqrt(h) * s12 * (W(:,i) - W(:,i-1));
        end

        phi1 = phi1 + 3 * X(1, 1:h/dt:end) + 2 * X(2, 1:h/dt:end);

        [xh, ~] = SGD(x_0, M, y, T, h);
        phi2 = phi2 + 3 * xh(1,:) + 2 * xh(2,:);
    end

    phi1 = phi1 / 1000;
    phi2 = phi2 / 1000;

    err(k) = max(abs(phi1 - phi2));
end

figure;
loglog(hvec, err, '-o');
hold on;
loglog(hvec, hvec, '--', 'DisplayName', 'h');
loglog(hvec, hvec.^2, '--', 'DisplayName', 'h^2');
legend("Result", "h", "h^2",location='southeast',fontsize=9);
hold off;

set(gcf,'Units','centimeters','Position',[2,2,10,6]);
exportgraphics(gcf,'q7.png','Resolution',600)