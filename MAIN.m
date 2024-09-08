clc
clear 
rng default
warning('off','all')
start_run_clock = clock;
start = tic;

%% -------------------�Ĥ@���q:��z���-----------------------------
D = xlsread('Data');   % �ƭȳ���
data=sortrows(D,[1 6 5]);        % �Ƨǥثe���, �������̦����Ԩ�̫e��
clear D                          % �A�̧ǥΪ��컷�V & �������������~�Ƨ�

% �C�@���Ƶ��@�ӦW�� (�u�O���ШC�@��O������)
% today = data(:,1);             % ������
% strike = data(:,2);            % K, strike price
% option = data(:,3);            % ����v��������
% spot = data(:,4);              % S, spot price
% SoverK = data(:,5);            % S/K
% option_maturity = data(:,6);   % ����v�Z����Ѽ�
% r = data(:,7);                 % r, interest rate

% �]�w�_�l��
parameters = [];
% parameters.SV = [0.0008, 4, 0.4, 1.2, -0.2]'; 
                % mu,kappa,theta,sigma,rho 
% parameters.SI = [0.0008, 4, 0.4, 1.2]'; 
                % mu, kappa, theta, sigma, rho, lambda, mu_J, sigma_J
% parameters.SVSI = [0.0008, 4, 0.4, 1.2, 4, 0.4, 1.2, -0.2]';
%                %mu, kappa_v, theta_v, sigma_v, kappa_r. theta_r, sigma_r,
%                rho

% �z��:�R���{�f(S)���欰 NaN �����
data(isnan(data(:,4)),:)=[];
% �z��:�R������v���欰 NaN �����
data(isnan(data(:,3)),:)=[];
% �z��:�R������v����C��0.02�����
data(data(:,3)<0.2,:) = [];
% �z��:�R������v����j��500�����
data(data(:,3)>500,:) = [];
% �z��:�R������v�Z�����b365�H�W�����
Max_days = 365;
data(data(:,6)>Max_days,:) = [];     
% �z��:�R������v�Z�����C��3�Ѫ����
data(data(:,6)<3,:) = [];    
% �z��:�R��S/K > 2 and S/K < 0.5 �� option
data(data(:,5)<0.5,:) = [];    
data(data(:,5)>2,:) = [];    
% �z��:�R��implied vol is nan 
BS_IV = blsimpv(data(:,4), data(:,2), data(:,7), data(:,6)/250, data(:,3));
data(isnan(BS_IV),:) = [];  


total_before=0;
j=0;                % �o�̪� j �O�ΨӰO�����X�Ѧ������� 
sum_h = 0;
unique_date = unique(data(:,1));
period_days = length(unique_date);
date_numbers_spot= zeros(period_days,4);  % ���]�w�@�Ӥj��������

% ���]2010/12/01�}�l��30��(period_days)�ѡA�䤤�u��22��,������Aj�̫�N�|�O22
for m=1:period_days 
    % �̥������I�����
    oneday_data = data( data(:,1)==unique_date(m), :);   % 7�����
    h = length(oneday_data(:,1));  % h �O�ΨӰO���@�Ѥ����X�������
    if h < 15                      % �p�G��ѨS����ε��ƹL�ִN���L�A��U�@��
        continue
    else
        j = j+1;
        % �O���e�@�ѩM��Ѧ��X��������
        total_after = total_before+h;
        % �ΨӰO�����@�ӥ����A���X����ơA�ֿn�X����ơA���{�f����  
        date_numbers_spot(j,:) = [oneday_data(1,1),h,total_after,oneday_data(1,4)];
        Date = datestr(oneday_data(1,1),26);
        fprintf('%s,\t\t%d,\t\t%d \n',Date,h,total_after);
        % ��z�L���C���ơA�����Ƨ�(�C�Ѥ@��)�s�Jperiod_data
        data((1+total_before):total_after,:) = oneday_data;
        total_before = total_after;
    end
end
date_numbers_spot(date_numbers_spot(:,1)==0,:) = []; % �R���h�l��������
available_data_days = j;  % ��Ƥ����X�ѬO�i�H�ΨӹB�⪺

%% ---------------�ĤG���q: �ѼƮզ�(calibration)�A�o�@�϶��̶O��----------------------- 
% ��ƾ�z�n����A���@�w�n��������Ƴ��]�A�i�H�]��������l
% ���w�]calibration 2/17~2/22�Ѹ�ơA�|�A�ư�����
Excution_start = datenum('2023/04/10');
Excution_end = datenum('2023/04/20');

start_index = find(date_numbers_spot(:,1)>=Excution_start, 1 );
end_index = find(date_numbers_spot(:,1)<=Excution_end,1,'last');

M = end_index - start_index + 1;
parameter_list_BS = zeros(M,3);   % 2+1��BS model�Ѽ�
parameter_list_SV = zeros(M,7);   % 2+5��SV model�Ѽ�
parameter_list_SI = zeros(M,6);   % 2+4��SI model�Ѽ�
parameter_list_SVSI = zeros(M,10); % 2+8��SVSI model�Ѽ�

for m = 1:M
    date1 = date_numbers_spot(start_index+m-1,1);
    datestr(date1, 26)    % ��ܼ˥��������
    oneday_data = data( data(:,1)==date1, :); 
    
    % 0. Black-Scholes model (�u���implied vol. �F)
    [BS_implied_V, BS_para_se] = calibration(oneday_data, [], 0);

    % 1.SV model 
    % input���ĤG���q�N��O���Ӥ@�Ӽҫ�1:SV, 2:SVJ
    % output���Ĥ@���q�G�Ĥ@�C��Ѽ�,�ĤG�C����p�зǻ~.  �ĤG���q�Gimplied vol �� SSE
    [SV_para, SV_impvol_SSE] = calibration(oneday_data, parameters, 1);

    % 2.SI model
    [SI_para, SI_impvol_SSE] = calibration(oneday_data, parameters, 2);

     % 2.SVSI model
    [SVSI_para, SVSI_impvol_SSE] = calibration(oneday_data, parameters, 3);

    % �ưѼƪ�. ������.�A��Ѽ�,...,�̫��implied vol �� SSE
    parameter_list_BS(m,:)   = [date1, reshape(BS_implied_V,1,[]), BS_para_se];
    parameter_list_SV(m,:)   = [date1, reshape(SV_para,1,[]), SV_impvol_SSE  ];
    parameter_list_SI(m,:)  = [date1, reshape(SI_para,1,[]), SI_impvol_SSE ];
    parameter_list_SVSI(m,:)  = [date1, reshape(SVSI_para,1,[]), SVSI_impvol_SSE ];
end

format short g
tEnd = toc(start);
disp('�ѼƮզ�����')
fprintf('Elapsed %d hours, %d minutes and %d seconds.\n',floor(tEnd/3600),floor(rem(tEnd/60,60)),round(rem(tEnd,60)));

%% -----------------------------�ĤT���q:�D�C�骺�˥�"�~"�w���~�t------------------------------------
before = 0;
alpha = 0.5;
Excution_start = datenum('2023/04/10');
Excution_end = datenum('2023/04/20');

start_index = find(date_numbers_spot(:,1)>=Excution_start, 1 );
end_index = find(date_numbers_spot(:,1)<=Excution_end,1,'last');

M = end_index - start_index + 1;
for m = 2:M
    % ��"�Q��"���ѼƥX��
    x_BS   = parameter_list_BS(m-1,  2);
    x_SV   = parameter_list_SV(m-1,  2:6);
    x_SI  = parameter_list_SI(m-1, 2:5);
    x_SVSI  = parameter_list_SVSI(m-1, 2:9);

    % �⤵��(1��)����ƥX��
    date2 = date_numbers_spot(start_index+m-1,1);
    datestr(date2, 26)    % ��ܼ˥��~�����
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

disp('�w���~�t����')
complete = clock
tEnd2 = toc(start);
fprintf('Elapsed time is %d hours, %d minutes and %d seconds.\n',...
    floor(tEnd2/3600),floor(rem(tEnd2/60,60)),round(rem(tEnd2,60)));
%% plot
A = categorical({'<30','30~60','60~90','90~120','120~180','>180'});
A = reordercats(A,{'<30','30~60','60~90','90~120','120~180','>180'});
figure(1) 
bar(A,[outsamp_BS_error(1,1:6);outsamp_SV_error(1,1:6);outsamp_SI_error(1,1:6);outsamp_SVSI_error(1,1:6)]) % ���� relative error
legend('BS','SV','SI','SVSI','location','southeast')
figure(2)
bar(A,[outsamp_BS_error(3,1:6);outsamp_SV_error(3,1:6);outsamp_SI_error(3,1:6);outsamp_SVSI_error(3,1:6)]) % ���� relative error
legend('BS','SV','SI','SVSI','location','southeast')
figure(3)
bar(A,[outsamp_BS_error(1,7:12);outsamp_SV_error(1,7:12);outsamp_SI_error(1,7:12);outsamp_SVSI_error(1,7:12)]) % ���� absolute error
legend('BS','SV','SI','SVSI','location','northwest')
figure(4)
bar(A,[outsamp_BS_error(3,7:12);outsamp_SV_error(3,7:12);outsamp_SI_error(3,7:12);outsamp_SVSI_error(3,7:12)]) % ���� absolute error
legend('BS','SV','SI','SVSI','location','northwest')
%% Ū�� Excel �ɮ�
filename = 'Data.xlsx'; % �бN 'your_filename.xlsx' ��������ڪ��ɮצW�٩θ��|
data = xlsread(filename);

% �z��j��1.05�����
filteredData = data(data(:, 5) > 1.05, :);
%disp(filteredData)
count = sum(filteredData(:, 3) >500);
disp(count);
