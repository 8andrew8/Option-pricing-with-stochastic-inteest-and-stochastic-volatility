function SSE = SI_FFT(x, data, alpha, largeNumber, pricing)

    kappa  = x(1);
    theta  = x(2);
    sigma  = x(3);
    v0      = x(4);
    eta    = 0;
    T            = data(:,6)/250;
    marketoption = data(:,3);
    K            = data(:,2);    % 向量
    S            = data(1,4);    % 純量
    r            = data(1,7);    % 純量

    % SI model
    model_price = SI_Options( S, v0, K, r, T, kappa, theta ...
        , sigma, alpha, largeNumber);

    if pricing == 1
        SSE = model_price;
        return
    end

%     temp = [marketoption, model_price, marketoption - model_price]
%     pause
    SSE = (marketoption - model_price)./marketoption ;    
    SSE = sum(SSE.^2);    
    % 使用fmincon要回傳一個值, 若使用lsqnonlin要回傳一列向量
end
%% function FFT_call = SI_Options( S, v, K, r, tau, kappa_r, theta_r, eta, sigma_r, alpha, largeNumber )
function FFT_call = SI_Options( S, v, K, r, tau, kappa_r, theta_r , sigma_r, alpha, largeNumber)
    k = log(K);
    f = @(x) exp(-1i*x.*k).*(exp(-r*tau).*cf_SI(x-1i*(alpha+1), S, v, r, tau, kappa_r, theta_r , sigma_r ))./(alpha^2+alpha-x.^2+1i*(2*alpha+1)*x);
    FFT_call = exp(-alpha*k)/pi.*real(integral(f, 0, largeNumber, 'ArrayValued', true));
end

%% Characteristic Function for SI model
function value = cf_SI(phi ,S, v, r, tau, kappa_r, theta_r, sigma_r)
    s = log(S);
    A = @(phi) -1i*phi*kappa_r*theta_r*tau./(sigma_r^2/2-kappa_r)+...
        ((1i*phi*kappa_r*theta_r)*(exp((sigma_r^2/2-kappa_r)*tau)-1)./(theta_r^2/2-kappa_r)^2);
    B = @(phi) -1i*phi*(1-exp(sigma_r^2*tau./2-kappa_r*tau))./(sigma_r^2/2-kappa_r);

    value = exp(1i*phi*s + 1i*phi*(1i*phi-1)*v*tau./2 + A(phi) + B(phi)*r);
end