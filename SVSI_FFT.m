function SSE = SVSI_FFT(x, data, alpha, largeNumber, pricing)

    kappa_v    = x(1);
    theta_v    = x(2);
    sigma_v    = x(3);
    kappa_r    = x(4);
    theta_r    = x(5);
    sigma_r    = x(6);
    rho        = x(7);
    v0         = x(8); 
    eta      = 0;

    T            = data(:,6)/250;
    marketoption = data(:,3);
    K            = data(:,2);    % 向量
    S            = data(1,4);    % 純量
    r            = data(1,7);    % 純量

    % SVSI model
    model_price = SVSI_Options( S, v0, log(K), r, T, kappa_r, theta_r, ...
        kappa_v, theta_v, eta, sigma_r, sigma_v, rho, alpha, largeNumber);

    if pricing == 1
        SSE = model_price;
        return
    end        
    
    SSE = (marketoption - model_price)./marketoption ;    
    SSE = sum(SSE.^2);    
    % 使用fmincon要回傳一個值, 若使用lsqnonlin要回傳一列向量

end
%% function FFT_call =SVSI_Options( S, v, log(K), r, tau, kappa_r, theta_r, kappa_v, theta_v, eta, sigma_r, sigma_v, rho, alpha, largeNumber )
function FFT_call = SVSI_Options( S, v, k, r, tau, kappa_r, theta_r, kappa_v, theta_v, eta, sigma_r, sigma_v, rho, alpha, largeNumber)
    f = @(x) exp(-1i*x.*k).*(exp(-r*tau).*cf_SVSI(x-1i*(alpha+1),S, v, r, tau, kappa_r, theta_r, kappa_v, theta_v, eta, sigma_r, sigma_v, rho))./(alpha^2+alpha-x.^2+1i*(2*alpha+1)*x);
    FFT_call = exp(-alpha*k)/pi.*real(integral(f, 0, largeNumber, 'ArrayValued',true));
end

%% Characteristic Function for SVSI model
function value = cf_SVSI(phi ,S, v, r, tau, kappa_r, theta_r, kappa_v, theta_v, eta, sigma_r, sigma_v, rho)
%     s = log(S);
%     d_r = @(phi) sqrt(2*sigma_r^2*(1i*phi-1)-kappa_r^2);
%     d_v = @(phi) sqrt((rho*sigma_v*1i*phi-kappa_v-eta).^2-sigma_v^2*(-1i*phi-phi.^2)) ;
%     g_r = @(phi) (kappa_r+eta+d_r(phi)) ./ (kappa_r+eta-d_r(phi));
%     g_v = @(phi) (kappa_v+eta-rho*sigma_v*1i*phi+d_v(phi)) ./ (kappa_v+eta-rho*sigma_v*1i*phi-d_v(phi));
%     A = @(phi) r*tau + kappa_v*theta_v/sigma_v^2*( (kappa_v+eta-rho*sigma_v*1i*phi+d_v(phi))*tau...
%                -2*log( (1-g_v(phi).*exp(d_v(phi)*tau))./(1-g_v(phi))) ) + kappa_r*theta_r/sigma_r^2 ...
%                *( (kappa_r+eta+d_r(phi))*tau - 2*log( (1-g_r(phi).*exp(d_r(phi)*tau))./(1-g_r(phi))) ) ;
%     B = @(phi) (kappa_v+eta-rho*sigma_v*1i*phi+d_v(phi))./sigma_v^2.*(1-exp(d_v(phi)*tau))./(1-g_v(phi).*exp(d_v(phi)*tau));
%     C = @(phi) (1i*phi-1)*(1-exp(d_r(phi)*tau)).*(1-g_r(phi))./(-d_r(phi))./(1-g_r(phi)*exp(d_r(phi)*tau));
%     
%     value = exp(1i*phi*s + A(phi) + B(phi)*v + C(phi)*r);
    s = log(S);
    d = @(phi) sqrt((rho*sigma_v*1i*phi-kappa_v-eta).^2 + sigma_v^2*(1i*phi + phi.^2)) ;
    g = @(phi) (kappa_v+eta-rho*sigma_v*1i*phi+d(phi)) ./ (kappa_v+eta-rho*sigma_v*1i*phi-d(phi));
    C = @(phi) kappa_v*theta_v/sigma_v^2*((kappa_v + eta - rho*sigma_v*1i*phi + d(phi))*tau ... 
        -2*log( (1-g(phi).*exp(d(phi)*tau))./(1 - g(phi))) ) - kappa_r*theta_r*1i*phi*tau/(1/2*sigma_r^2-kappa_r)...
        + kappa_r*theta_r*1i*phi/(1/2*sigma_r^2-kappa_r)^2*(exp(1/2*sigma_r^2-kappa_r)*tau-1);
    D = @(phi) (kappa_v + eta - rho*sigma_v*1i*phi + d(phi))./sigma_v^2.*(1 - exp(d(phi)*tau))./(1 - g(phi).*exp(d(phi)*tau));
    E = @(phi) (-1i*phi)*(1 - exp(1/2*sigma_r^2*tau-kappa_r*tau))/(1/2*sigma_r^2-kappa_r);
    
    value = exp(1i*phi*s + C(phi) + D(phi)*v + E(phi)*r);
end
