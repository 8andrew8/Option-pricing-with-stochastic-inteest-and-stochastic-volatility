function [estimated_parameters, SSE] = calibration(data, parameters, model)
warning('off','all')
alpha = 0.05

% default options
% 1. MaxFunEvals = 100*numberOfVariables, Maximum number of function evaluations allowed.
% 2. MaxIter = 400, Maximum number of iterations allowed.
% 3. TolFun = 1e-6, Termination tolerance on the function value.
% 4. TolX = 1e-6, Termination tolerance on x.
% 5. Algorithm: 
%    fmincon:'interior-point'(default), 'trust-region-reflective', 'sqp', 'active-set'.
% Examples:
% options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-10);
% options = optimset('Display','iter','MaxIter',40,'Algorithm','active-set');
% SV: options = optimset('MaxIter',100,'TolFun',1e-5,'TolX',1e-5);
% SVJ: options = optimset('MaxIter',120,'TolFun',1e-5,'TolX',1e-5);

switch model
    case 0 
        %BS
        BS_implied_vol = blsimpv(data(:,4), data(:,2), ...
                                 data(:,7), data(:,6)/250, data(:,3));
                        % blsimpv(S, K, r, T, option)
        x(1,1) = mean(BS_implied_vol);
        SSE = sum( (BS_implied_vol-x(1,1)).^2 );

    case 1
        %SV        
        x0 = [1, 0.1, 0.1,  0,  0.1];
        %       kappa  theta sigma  rho  v0
        ub =   [ 150,    2,    2,    0.9,   2  ];
        lb =   [ eps,   eps,  eps,  -0.9,   eps  ];
        options = [];
        [x, SSE] = fmincon(@SV_FFT, x0, [], [], [], [], lb, ub, [], options, data, alpha, 200, 0);
%         [x, SSE] = lsqnonlin(@SV_FFT, x0, lb, ub, options, data, alpha, 100, 0);
        
    case 2
        %SI
        
        x0 = [1, 0.1, 0.1, 0.1];
        %       kappa  theta sigma    v0
        ub =   [ 150,    2,    2,     2  ];
        lb =   [ eps,   eps,  eps,   eps  ];
        options = [];
        [x, SSE] = fmincon(@SI_FFT, x0, [], [], [], [], lb, ub, [], options, data, alpha, 200, 0);
%         [x, SSE] = lsqnonlin(@SI_FFT, x0, lb, ub, options, data, alpha, 100, 0);
    
    case 3
        %SVSI

        x0 = [1, 0.1, 0.1, 1, 0.01, 0.01, 0,  0.1];
        %       kappa_v, theta_v, sigma_v, kappa_r. theta_r, sigma_r, rho,  v0
        ub =   [ 150,      2,      2,       150,        1,      1,     0.9,  2  ];
        lb =   [ eps,    eps,     eps,       0,       eps,     eps,   -0.9, eps  ];
        options = [];
        [x, SSE] = fmincon(@SVSI_FFT, x0, [], [], [], [], lb, ub, [], options, data, alpha, 200, 0);
%         [x, SSE] = lsqnonlin(@SVSI_FFT, x0, lb, ub, options, data, alpha, 100, 0);
end

model
estimated_parameters = x
SSE

end