function SSE = SVJ_FFT(x, data, alpha, largeNumber, pricing)

    kappa    = x(1);
    theta    = x(2);
    sigma    = x(3);
    rho      = x(4);
    lambda_S = x(5);
    JS_mu    = x(6);
    JS_sigma = x(7);
    v0        = x(8); 
    eta      = 0;

    T            = data(:,6)/250;
    marketoption = data(:,3);
    K            = data(:,2);    % 向量
    S            = data(1,4);    % 純量
    r            = data(1,7);    % 純量


    % SVJ model
    model_price = SVJ_Options( log(S), v0, log(K), r, T, kappa,...
        theta, eta,sigma, rho, lambda_S,...
        JS_mu, JS_sigma, alpha, largeNumber);

    if pricing == 1
        SSE = model_price;
        return
    end        
    
    SSE = (marketoption - model_price)./marketoption ;    
    SSE = sum(SSE.^2);    
    % 使用fmincon要回傳一個值, 若使用lsqnonlin要回傳一列向量

end
%% function FFT_call =SVJ_Options( log(S), v, log(K), r, tau, kappa, theta, eta, sigma, rho, lambda_S, JS_mu, JS_sigma, alpha, largeNumber )
function FFT_call = SVJ_Options( s, v, k, r, tau, kappa, theta, eta, sigma, rho, lambda_S, JS_mu, JS_sigma, alpha, largeNumber)
    J_bar = exp(JS_mu+0.5*JS_sigma^2)-1;
    f = @(x) exp(-1i*x.*k).*(exp(-r*tau).*cf_SVJ(x-1i*(alpha+1),s, v, r, tau, kappa, theta, eta,...
        sigma, rho, lambda_S, JS_mu, JS_sigma, J_bar))./(alpha^2+alpha-x.^2+1i*(2*alpha+1)*x);
    FFT_call = exp(-alpha*k)/pi.*real(integral(f, 0, largeNumber, 'ArrayValued',true));
end

%% Characteristic Function for SVJ model
function value = cf_SVJ(phi ,s, v, r, tau, kappa, theta, eta, sigma, rho, lambda_S, JS_mu, JS_sigma, J_bar)
    u = 1i*phi;
    d = sqrt((rho*sigma*u-kappa-eta).^2-sigma^2*(-u-phi.^2)) ;
    y = exp(d*tau);
    x = kappa+eta-rho*sigma*u+d;
    sigma_2 = sigma^2;
    
    g = x ./ (x-d-d);
    A = ( (r-lambda_S*J_bar)*u + lambda_S*(exp(u*JS_mu+0.5*(JS_sigma*u).^2)-1) )*tau + ...
             kappa*theta/sigma_2*( x*tau -2*log( (1-g.*y)./(1-g)) );
    B = x./sigma_2.*(1-y)./(1-g.*y);
    
    value = exp(u*s + A + B*v);
end
