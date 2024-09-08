function SSE = SV_FFT(x, data, alpha, largeNumber, pricing)

    kappa  = x(1);
    theta  = x(2);
    sigma  = x(3);
    rho    = x(4);
    v0      = x(5);
    eta    = 0;
    T            = data(:,6)/250;
    marketoption = data(:,3);
    K            = data(:,2);    % 向量
    S            = data(1,4);    % 純量
    r            = data(1,7);    % 純量

    % SV model
    model_price = SV_Options( S, v0, K, r, T, kappa, theta, eta, ...
            sigma, rho, alpha, largeNumber);

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

%% function FFT_call = SV_Options( S, v, K, r, tau, kappa_v, theta_v, eta, sigma_v, rho, alpha, largeNumber )
function FFT_call = SV_Options( S, v, K, r, tau, kappa_v, theta_v, eta, sigma_v, rho, alpha, largeNumber)
    k = log(K);
    f = @(x) exp(-1i*x.*k).*(exp(-r*tau).*cf_SV(x-1i*(alpha+1), S, v, r, tau, kappa_v, theta_v, eta, sigma_v, rho))./(alpha^2+alpha-x.^2+1i*(2*alpha+1)*x);
    FFT_call = exp(-alpha*k)/pi.*real(integral(f, 0, largeNumber, 'ArrayValued',true));
end

%% Characteristic Function for SV model
function value = cf_SV(phi ,S, v, r, tau, kappa_v, theta_v, eta, sigma_v, rho)
    s = log(S);
    d = @(phi) sqrt((rho*sigma_v*1i*phi-kappa_v-eta).^2-sigma_v^2*(-1i*phi-phi.^2)) ;
    g = @(phi) (kappa_v+eta-rho*sigma_v*1i*phi+d(phi)) ./ (kappa_v+eta-rho*sigma_v*1i*phi-d(phi));
    A = @(phi) r*1i*phi*tau + kappa_v*theta_v/sigma_v^2*( (kappa_v+eta-rho*sigma_v*1i*phi+d(phi))*tau...
               -2*log( (1-g(phi).*exp(d(phi)*tau))./(1-g(phi))) );
    B = @(phi) (kappa_v+eta-rho*sigma_v*1i*phi+d(phi))./sigma_v^2.*(1-exp(d(phi)*tau))./(1-g(phi).*exp(d(phi)*tau));

    value = exp(1i*phi*s + A(phi) + B(phi)*v);
end
