clc
clear 
rng default
warning('off','all')
start_run_clock = clock;
start = tic;

%% -------------------第一階段:整理資料-----------------------------
D = xlsread('Data');   % 數值部分
data=sortrows(D,[1 6 5]);        % 排序目前資料, 交易日期最早的拉到最前面
clear D                          % 再依序用近月到遠越 & 價內到價平到價外排序

% 每一行資料給一個名稱 (只是介紹每一行是什麼資料)
% today = data(:,1);             % 交易日期
% strike = data(:,2);            % K, strike price
% option = data(:,3);            % 選擇權市場價格
% spot = data(:,4);              % S, spot price
% SoverK = data(:,5);            % S/K
% option_maturity = data(:,6);   % 選擇權距到期天數
% r = data(:,7);                 % r, interest rate

% 設定起始值
parameters = [];
% parameters.SV = [0.0008, 4, 0.4, 1.2, -0.2]'; 
                % mu,kappa,theta,sigma,rho 
% parameters.SI = [0.0008, 4, 0.4, 1.2]'; 
                % mu, kappa, theta, sigma, rho, lambda, mu_J, sigma_J
% parameters.SVSI = [0.0008, 4, 0.4, 1.2, 4, 0.4, 1.2, -0.2]';
%                %mu, kappa_v, theta_v, sigma_v, kappa_r. theta_r, sigma_r,
%                rho

% 篩選:刪除現貨(S)價格為 NaN 的資料
data(isnan(data(:,4)),:)=[];
% 篩選:刪除選擇權價格為 NaN 的資料
data(isnan(data(:,3)),:)=[];
% 篩選:刪除選擇權價格低於0.02的資料
data(data(:,3)<0.2,:) = [];
% 篩選:刪除選擇權價格大於500的資料
data(data(:,3)>500,:) = [];
% 篩選:刪除選擇權距到期日在365以上的資料
Max_days = 365;
data(data(:,6)>Max_days,:) = [];     
% 篩選:刪除選擇權距到期日低於3天的資料
data(data(:,6)<3,:) = [];    
% 篩選:刪除S/K > 2 and S/K < 0.5 的 option
data(data(:,5)<0.5,:) = [];    
data(data(:,5)>2,:) = [];    
% 篩選:刪除implied vol is nan 
BS_IV = blsimpv(data(:,4), data(:,2), data(:,7), data(:,6)/250, data(:,3));
data(isnan(BS_IV),:) = [];  


total_before=0;
j=0;                % 這裡的 j 是用來記錄有幾天有交易資料 
sum_h = 0;
unique_date = unique(data(:,1));
period_days = length(unique_date);
date_numbers_spot= zeros(period_days,4);  % 先設定一個大的紀錄區

% 假設2010/12/01開始的30天(period_days)天，其中只有22天,有交易，j最後就會是22
for m=1:period_days 
    % 依交易日期截取資料
    oneday_data = data( data(:,1)==unique_date(m), :);   % 7個欄位
    h = length(oneday_data(:,1));  % h 是用來記錄一天內有幾筆的資料
    if h < 15                      % 如果當天沒交易或筆數過少就跳過，找下一天
        continue
    else
        j = j+1;
        % 記錄前一天和當天有幾筆交易資料
        total_after = total_before+h;
        % 用來記錄哪一個交易日，有幾筆資料，累積幾筆資料，當日現貨價格  
        date_numbers_spot(j,:) = [oneday_data(1,1),h,total_after,oneday_data(1,4)];
        Date = datestr(oneday_data(1,1),26);
        fprintf('%s,\t\t%d,\t\t%d \n',Date,h,total_after);
        % 把篩過的每日資料，交易日排序(每天一組)存入period_data
        data((1+total_before):total_after,:) = oneday_data;
        total_before = total_after;
    end
end
date_numbers_spot(date_numbers_spot(:,1)==0,:) = []; % 刪掉多餘的紀錄區
available_data_days = j;  % 資料內有幾天是可以用來運算的

%% ---------------第二階段: 參數校估(calibration)，這一區塊最費時----------------------- 
% 資料整理好之後，不一定要全部的資料都跑，可以跑部分的日子
% 先預設calibration 2/17~2/22天資料，會再排除假日
Excution_start = datenum('2023/04/10');
Excution_end = datenum('2023/04/20');

start_index = find(date_numbers_spot(:,1)>=Excution_start, 1 );
end_index = find(date_numbers_spot(:,1)<=Excution_end,1,'last');

M = end_index - start_index + 1;
parameter_list_BS = zeros(M,3);   % 2+1個BS model參數
parameter_list_SV = zeros(M,7);   % 2+5個SV model參數
parameter_list_SI = zeros(M,6);   % 2+4個SI model參數
parameter_list_SVSI = zeros(M,10); % 2+8個SVSI model參數

for m = 1:M
    date1 = date_numbers_spot(start_index+m-1,1);
    datestr(date1, 26)    % 顯示樣本內的日期
    oneday_data = data( data(:,1)==date1, :); 
    
    % 0. Black-Scholes model (只能估implied vol. 了)
    [BS_implied_V, BS_para_se] = calibration(oneday_data, [], 0);

    % 1.SV model 
    % input的第二分量代表是哪個一個模型1:SV, 2:SVJ
    % output的第一分量：第一列放參數,第二列放估計標準誤.  第二分量：implied vol 的 SSE
    [SV_para, SV_impvol_SSE] = calibration(oneday_data, parameters, 1);

    % 2.SI model
    [SI_para, SI_impvol_SSE] = calibration(oneday_data, parameters, 2);

     % 2.SVSI model
    [SVSI_para, SVSI_impvol_SSE] = calibration(oneday_data, parameters, 3);

    % 排參數表. 先放日期.再放參數,...,最後放implied vol 的 SSE
    parameter_list_BS(m,:)   = [date1, reshape(BS_implied_V,1,[]), BS_para_se];
    parameter_list_SV(m,:)   = [date1, reshape(SV_para,1,[]), SV_impvol_SSE  ];
    parameter_list_SI(m,:)  = [date1, reshape(SI_para,1,[]), SI_impvol_SSE ];
    parameter_list_SVSI(m,:)  = [date1, reshape(SVSI_para,1,[]), SVSI_impvol_SSE ];
end

format short g
tEnd = toc(start);
disp('參數校估結束')
fprintf('Elapsed %d hours, %d minutes and %d seconds.\n',floor(tEnd/3600),floor(rem(tEnd/60,60)),round(rem(tEnd,60)));

%% -----------------------------第三階段:求每日的樣本"外"定價誤差------------------------------------
before = 0;
alpha = 0.5;
Excution_start = datenum('2023/04/10');
Excution_end = datenum('2023/04/20');

start_index = find(date_numbers_spot(:,1)>=Excution_start, 1 );
end_index = find(date_numbers_spot(:,1)<=Excution_end,1,'last');

M = end_index - start_index + 1;
for m = 2:M
    % 抽"昨天"的參數出來
    x_BS   = parameter_list_BS(m-1,  2);
    x_SV   = parameter_list_SV(m-1,  2:6);
    x_SI  = parameter_list_SI(m-1, 2:5);
    x_SVSI  = parameter_list_SVSI(m-1, 2:9);

    % 抽今天(1日)的資料出來
    date2 = date_numbers_spot(start_index+m-1,1);
    datestr(date2, 26)    % 顯示樣本外的日期
    oneday_data = data(data(:,1)==date2, :); 
    after = before + length(oneday_data(:,1));
    classify_info_out((1+before):after,1:4) = oneday_data(:,[1,5,6,3]);
    outsamp_BS((1+before):after,1)   = blsprice(oneday_data(:,4), oneday_data(:,2), oneday_data(:,7),...
                                                oneday_data(:,6)./250, x_BS);    
    outsamp_SV((1+before):after,1)   = SV_FFT(x_SV, oneday_data, alpha, 50, 1);
    outsamp_SI((1+before):after,1)  = SI_FFT(x_SI, oneday_data, alpha, 50, 1);
    outsamp_SVSI((1+before):after,1)  = SVSI_FFT(x_SVSI, oneday_data, alpha, 50, 1);
    before = after;
end

% The input data: 1.date, 2.SoverK, 3.Maturity, 4.market option price. & out-sample theoretical model price.
outsamp_BS_error   = error_classify([classify_info_out, outsamp_BS]);
outsamp_SV_error   = error_classify([classify_info_out, outsamp_SV]);
outsamp_SI_error  = error_classify([classify_info_out, outsamp_SI]);
outsamp_SVSI_error  = error_classify([classify_info_out, outsamp_SVSI]);

disp('定價誤差結束')
complete = clock
tEnd2 = toc(start);
fprintf('Elapsed time is %d hours, %d minutes and %d seconds.\n',...
    floor(tEnd2/3600),floor(rem(tEnd2/60,60)),round(rem(tEnd2,60)));
%% plot
A = categorical({'<30','30~60','60~90','90~120','120~180','>180'});
A = reordercats(A,{'<30','30~60','60~90','90~120','120~180','>180'});
figure(1) 
bar(A,[outsamp_BS_error(1,1:6);outsamp_SV_error(1,1:6);outsamp_SI_error(1,1:6);outsamp_SVSI_error(1,1:6)]) % 價內 relative error
legend('BS','SV','SI','SVSI','location','southeast')
figure(2)
bar(A,[outsamp_BS_error(3,1:6);outsamp_SV_error(3,1:6);outsamp_SI_error(3,1:6);outsamp_SVSI_error(3,1:6)]) % 價平 relative error
legend('BS','SV','SI','SVSI','location','southeast')
figure(3)
bar(A,[outsamp_BS_error(1,7:12);outsamp_SV_error(1,7:12);outsamp_SI_error(1,7:12);outsamp_SVSI_error(1,7:12)]) % 價內 absolute error
legend('BS','SV','SI','SVSI','location','northwest')
figure(4)
bar(A,[outsamp_BS_error(3,7:12);outsamp_SV_error(3,7:12);outsamp_SI_error(3,7:12);outsamp_SVSI_error(3,7:12)]) % 價平 absolute error
legend('BS','SV','SI','SVSI','location','northwest')
%% 讀取 Excel 檔案
filename = 'Data.xlsx'; % 請將 'your_filename.xlsx' 替換為實際的檔案名稱或路徑
data = xlsread(filename);

% 篩選大於1.05的資料
filteredData = data(data(:, 5) > 1.05, :);
%disp(filteredData)
count = sum(filteredData(:, 3) >500);
disp(count);
