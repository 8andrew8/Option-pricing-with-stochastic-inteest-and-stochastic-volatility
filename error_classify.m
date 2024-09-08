function table = error_classify(raw_data)
% raw_data 有5行
% 1.date, 2.SoverK, 3.Maturity, 4.market option price, 5.theoretical price.
data(:,1) = (raw_data(:,5)-raw_data(:,4))./raw_data(:,4);
data(:,2) = abs(raw_data(:,5)-raw_data(:,4));
data(:,3) = raw_data(:,3);
data(:,4) = raw_data(:,2);

% data有4行
% 1:relative error, 2:absolute error, 3:maturity days, 4:SoverK
%左邊六行是是相對差，右邊六行是絕對誤差
    
%----資料依價內(A)、價平(B)、價外(C)分類------------------------------------
% 先宣告A.B.C.矩陣，因為若無資料，下面程式會出現錯誤
    A=zeros(1,4);    B=zeros(1,4);    C=zeros(1,4);
    k=1;    %記錄每一組的個數
    kk=1;
    kkk=1;
    for i=1:length(data(:,4))
        if (0.01<data(i,4))&&(data(i,4)<0.95) %價內買權
            A(k,:)=data(i,:);  
            k=k+1;
        end
        if (0.95<=data(i,4))&&(data(i,4)<1.05)  %價平
            B(kk,:)=data(i,:);
            kk=kk+1;
        end          
        if (1.05<=data(i,4)) %價外
            C(kkk,:)=data(i,:);
            kkk=kkk+1;
        end
    end
    
%----將價內、價平、價外的資料再根據"距到期日"進行細分------------------------
    % 先宣告小矩陣，因為若該小類無資料，下面程式會出現錯誤
    A_1=zeros(1,4);    A_2=zeros(1,4);    A_3=zeros(1,4);   A_4=zeros(1,4);
    A_5=zeros(1,4);    A_6=zeros(1,4);    B_1=zeros(1,4);   B_2=zeros(1,4);
    B_3=zeros(1,4);    B_4=zeros(1,4);    B_5=zeros(1,4);   B_6=zeros(1,4);
    C_1=zeros(1,4);    C_2=zeros(1,4);    C_3=zeros(1,4);   C_4=zeros(1,4);
    C_5=zeros(1,4);    C_6=zeros(1,4);

    v=ones(6,1);                     %用途與k,kk,kkk一樣, 只是有6組
    for i=1:length(A(:,1))
        if A(i,3)<=30                %距到期日30天以內
            A_1(v(1),1:4)=A(i,1:4);
            v(1)=v(1)+1;
        end    
        if (A(i,3)<=60)&&(A(i,3)>30)
            A_2(v(2),1:4)=A(i,1:4);
            v(2)=v(2)+1;
        end
        if (A(i,3)<=90)&&(A(i,3)>60)
            A_3(v(3),1:4)=A(i,1:4);
            v(3)=v(3)+1;
        end
        if (A(i,3)<=120)&&(A(i,3)>90)
            A_4(v(4),1:4)=A(i,1:4);
            v(4)=v(4)+1;
        end
        if (A(i,3)<=180)&&(A(i,3)>120)
            A_5(v(5),1:4)=A(i,1:4);
            v(5)=v(5)+1;
        end
        if A(i,3)>180
            A_6(v(6),1:4)=A(i,1:4);
            v(6)=v(6)+1;
        end
    end
    % 將價平資料依到期日細分
    v=ones(6,1);                    %先重設index
    for i=1:length(B(:,1))
        if B(i,3)<=30                
            B_1(v(1),1:4)=B(i,1:4);
            v(1)=v(1)+1;
        end    
        if (B(i,3)<=60)&&(B(i,3)>30)
            B_2(v(2),1:4)=B(i,1:4);
            v(2)=v(2)+1;
        end
        if (B(i,3)<=90)&&(B(i,3)>60)
            B_3(v(3),1:4)=B(i,1:4);
            v(3)=v(3)+1;
        end
        if (B(i,3)<=120)&&(B(i,3)>90)
            B_4(v(4),1:4)=B(i,1:4);
            v(4)=v(4)+1;
        end
        if (B(i,3)<=180)&&(B(i,3)>120)
            B_5(v(5),1:4)=B(i,1:4);
            v(5)=v(5)+1;
        end
        if B(i,3)>180
            B_6(v(6),1:4)=B(i,1:4);
            v(6)=v(6)+1;
        end
    end
    % 將價外資料，依到期日細分
    v=ones(6,1); 
    for i=1:length(C(:,1))
        if C(i,3)<=30                
            C_1(v(1),1:4)=C(i,1:4);
            v(1)=v(1)+1;
        end    
        if (C(i,3)<=60)&&(C(i,3)>30)
            C_2(v(2),1:4)=C(i,1:4);
            v(2)=v(2)+1;
        end
        if (C(i,3)<=90)&&(C(i,3)>60)
            C_3(v(3),1:4)=C(i,1:4);
            v(3)=v(3)+1;
        end
        if (C(i,3)<=120)&&(C(i,3)>90)
            C_4(v(4),1:4)=C(i,1:4);
            v(4)=v(4)+1;
        end
        if (C(i,3)<=180)&&(C(i,3)>120)
            C_5(v(5),1:4)=C(i,1:4);
            v(5)=v(5)+1;
        end
        if C(i,3)>180
            C_6(v(6),1:4)=C(i,1:4);
            v(6)=v(6)+1;
        end
    end

%----各小組 pricing error 平均與標準差--------------------------------------

    % relative error-----------------------------
    % 價內 (A)
    table(1,1)=mean(A_1(:,1));      % 到期日小於30天
    table(2,1)=std(A_1(:,1))/sqrt(length(A_1(:,1)));
    table(1,2)=mean(A_2(:,1));      % 到期日30~60
    table(2,2)=std(A_2(:,1))/sqrt(length(A_2(:,1)));
    table(1,3)=mean(A_3(:,1));      % 到期日60~90
    table(2,3)=std(A_3(:,1))/sqrt(length(A_3(:,1)));
    table(1,4)=mean(A_4(:,1));      % 到期日90~120
    table(2,4)=std(A_4(:,1))/sqrt(length(A_4(:,1)));
    table(1,5)=mean(A_5(:,1));      % 到期日120~180
    table(2,5)=std(A_5(:,1))/sqrt(length(A_5(:,1)));       
    table(1,6)=mean(A_6(:,1));      % 到期日180以上
    table(2,6)=std(A_6(:,1))/sqrt(length(A_6(:,1)));    
    % 價平 (B)
    table(3,1)=mean(B_1(:,1));
    table(4,1)=std(B_1(:,1))/sqrt(length(B_1(:,1)));
    table(3,2)=mean(B_2(:,1));
    table(4,2)=std(B_2(:,1))/sqrt(length(B_2(:,1)));
    table(3,3)=mean(B_3(:,1));
    table(4,3)=std(B_3(:,1))/sqrt(length(B_3(:,1)));
    table(3,4)=mean(B_4(:,1));
    table(4,4)=std(B_4(:,1))/sqrt(length(B_4(:,1)));
    table(3,5)=mean(B_5(:,1));
    table(4,5)=std(B_5(:,1))/sqrt(length(B_5(:,1)));
    table(3,6)=mean(B_6(:,1));
    table(4,6)=std(B_6(:,1))/sqrt(length(B_6(:,1)));  
    % 價外 (C)
    table(5,1)=mean(C_1(:,1));
    table(6,1)=std(C_1(:,1))/sqrt(length(C_1(:,1)));
    table(5,2)=mean(C_2(:,1));
    table(6,2)=std(C_2(:,1))/sqrt(length(C_2(:,1)));
    table(5,3)=mean(C_3(:,1));
    table(6,3)=std(C_3(:,1))/sqrt(length(C_3(:,1)));
    table(5,4)=mean(C_4(:,1));
    table(6,4)=std(C_4(:,1))/sqrt(length(C_4(:,1)));
    table(5,5)=mean(C_5(:,1));
    table(6,5)=std(C_5(:,1))/sqrt(length(C_5(:,1)));
    table(5,6)=mean(C_6(:,1));
    table(6,6)=std(C_6(:,1))/sqrt(length(C_6(:,1)));  
    
    % absolute error--------------------------
    % 價內 (A)
    table(1,7)=mean(A_1(:,2));
    table(2,7)=std(A_1(:,2))/sqrt(length(A_1(:,2)));
    table(1,8)=mean(A_2(:,2));
    table(2,8)=std(A_2(:,2))/sqrt(length(A_2(:,2)));
    table(1,9)=mean(A_3(:,2));
    table(2,9)=std(A_3(:,2))/sqrt(length(A_3(:,2)));
    table(1,10)=mean(A_4(:,2));
    table(2,10)=std(A_4(:,2))/sqrt(length(A_4(:,2)));
    table(1,11)=mean(A_5(:,2));
    table(2,11)=std(A_5(:,2))/sqrt(length(A_5(:,2)));
    table(1,12)=mean(A_6(:,2));
    table(2,12)=std(A_6(:,2))/sqrt(length(A_6(:,2)));
    % 價平 (B)
    table(3,7)=mean(B_1(:,2));
    table(4,7)=std(B_1(:,2))/sqrt(length(B_1(:,2)));
    table(3,8)=mean(B_2(:,2));
    table(4,8)=std(B_2(:,2))/sqrt(length(B_2(:,2)));
    table(3,9)=mean(B_3(:,2));
    table(4,9)=std(B_3(:,2))/sqrt(length(B_3(:,2)));
    table(3,10)=mean(B_4(:,2));
    table(4,10)=std(B_4(:,2))/sqrt(length(B_4(:,2)));
    table(3,11)=mean(B_5(:,2));
    table(4,11)=std(B_5(:,2))/sqrt(length(B_5(:,2)));
    table(3,12)=mean(B_6(:,2));
    table(4,12)=std(B_6(:,2))/sqrt(length(B_6(:,2)));  
    % 價外 (C)
    table(5,7)=mean(C_1(:,2));
    table(6,7)=std(C_1(:,2))/sqrt(length(C_1(:,2)));
    table(5,8)=mean(C_2(:,2));
    table(6,8)=std(C_2(:,2))/sqrt(length(C_2(:,2)));
    table(5,9)=mean(C_3(:,2));
    table(6,9)=std(C_3(:,2))/sqrt(length(C_3(:,2)));
    table(5,10)=mean(C_4(:,2));
    table(6,10)=std(C_4(:,2))/sqrt(length(C_4(:,2)));
    table(5,11)=mean(C_5(:,2));
    table(6,11)=std(C_5(:,2))/sqrt(length(C_5(:,2)));
    table(5,12)=mean(C_6(:,2));
    table(6,12)=std(C_6(:,2))/sqrt(length(C_6(:,2)));     

%     %----report 無資料的組別-------------------------
%     if A_1(1,4)==0     
%         fprintf('A_1\n')    
%     end
%     if A_2(1,4)==0     
%         fprintf('A_2\n')    
%     end
%     if A_3(1,4)==0     
%         fprintf('A_3\n')    
%     end
%     if A_4(1,4)==0     
%         fprintf('A_4\n')    
%     end
%     if A_5(1,4)==0     
%         fprintf('A_5\n')    
%     end
%     if A_6(1,4)==0     
%         fprintf('A_6\n')    
%     end
%     if B_1(1,4)==0     
%         fprintf('B_1\n')    
%     end
%     if B_2(1,4)==0     
%         fprintf('B_2\n')    
%     end
%     if B_3(1,4)==0     
%         fprintf('B_3\n')    
%     end
%     if B_4(1,4)==0     
%         fprintf('B_4\n')    
%     end
%     if B_5(1,4)==0     
%         fprintf('B_5\n')    
%     end
%     if B_6(1,4)==0     
%         fprintf('B_6\n')    
%     end
%     if C_1(1,4)==0     
%         fprintf('C_1\n')    
%     end
%     if C_2(1,4)==0     
%         fprintf('C_2\n')    
%     end
%     if C_3(1,4)==0     
%         fprintf('C_3\n')    
%     end
%     if C_4(1,4)==0     
%         fprintf('C_4\n')    
%     end
%     if C_5(1,4)==0     
%         fprintf('C_5\n')    
%     end
%     if C_6(1,4)==0     
%         fprintf('C_6\n')    
%     end
    
end