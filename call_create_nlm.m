ds = A_Ist_detail;
modelfun = @(b,x) x(:,1).*exp((-1*x(:,2)/b(1)).^2);
beta0 = [1100];

% ds = f_max_detail;
% modelfun = @(b,x) b(1) + b(2)./x(:,1) + b(3)./x(:,2) + b(4)./(x(:,1)).^2 + b(5)./(x(:,2)).^2 + b(6)./(x(:,1).*x(:,2));
% beta0 = [1000 1000 1000 1000 1000 1000];

model = create_nlm(ds, modelfun, beta0);