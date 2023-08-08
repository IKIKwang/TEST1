function [Bestfitness, Leader_pos, Curve] = MRFO_CS(nPop, MaxIt, Rs, Re, L, data, PopPos, PopFit, Low, Up, Dim)
 
[bestfitness, bestindex] = max(PopFit);
Leader_pos = PopPos(bestindex, :);
BestX= PopPos(bestindex, :);
Leader_score = bestfitness; 
ub=Up;

%% 初始结果显示
x = Leader_pos(1:2:end);
y = Leader_pos(2:2:end);
disp('初始位置：' );
for i = 1:Dim/2
    disp([num2str(x(i)), '     ', num2str(y(i))]);
end
disp(['初始覆盖率：', num2str(Leader_score)]);
% 初始覆盖图
figure
for i = 1:Dim/2
    axis([0 Up 0 Up]);            % 限制坐标范围
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    fill(x(i)+Rs*cos(sita), y(i)+Rs*sin(sita), 'w');
    %plot(x(i)+Rs*cos(sita), y(i)+Rs*sin(sita), 'b');
end
plot(x, y, 'k+');
str_title= ['Initial deployment coverage rate：',num2str(Leader_score)];
title({str_title},'FontSize',12);

t = 0;       % 迭代次数计数器
%% 迭代寻优
for It=1:MaxIt  
    Coef=It/MaxIt; 
    omega = 3;
    a = 2*(1-Coef);
    P = randn*((sin(pi/2 * Coef))^omega + cos(pi/2 * Coef)-1);  % 动态扰动因子策略
    M = a*(2*rand-1) + P;     
    if rand<0.5
          r1=rand;   
          %r = 1-exp(-(4*(MaxIt-It)/It)^2);
          Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
          if  abs(M)<1    %局部搜索
              %newPopPos(1,:)=BestX+r*(BestX-PopPos(1,:))+Beta*(BestX-PopPos(1,:)); 
              newPopPos(1,:)=BestX+rand(1,Dim).*(BestX-PopPos(1,:))+Beta*(BestX-PopPos(1,:)); %Equation (4)
          else
              %IndivRand=rand(1,Dim).*(Up-Low)+Low; 
              %newPopPos(1,:)=get_cuckoos(PopPos(1,:),BestX,Low,Up);
              IndivRand=get_cuckoos(PopPos(1,:),BestX,Low,Up);
              newPopPos(1,:)=IndivRand+rand(1,Dim).*(IndivRand-PopPos(1,:))+Beta*(IndivRand-PopPos(1,:)); %Equation (7)         
          end              
     else 
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(1,:)=PopPos(1,:)+rand(1,Dim).*(BestX-PopPos(1,:))+Alpha.*(BestX-PopPos(1,:)); %Equation (1)
     end
     
    for i=2:nPop
        if rand<0.5
           r1=rand;       
           %r = 1-exp(-(4*(MaxIt-It)/It)^2);
           Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
             if  abs(M)<1  
                 %newPopPos(i,:)=BestX+r*(BestX-PopPos(i-1,:))+Beta*(BestX-PopPos(i,:));
                 newPopPos(i,:)=BestX+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(BestX-PopPos(i,:)); %Equation (4)
             else
                 %IndivRand=rand(1,Dim).*(Up-Low)+Low;  
                 %newPopPos(i,:)=get_cuckoos(PopPos(i-1,:),BestX,Low,Up);
                 IndivRand=get_cuckoos(PopPos(i-1,:),BestX,Low,Up);                
                 newPopPos(i,:)=IndivRand+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(IndivRand-PopPos(i,:));  %Equation (7)       
             end              
        else
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(i,:)=PopPos(i,:)+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Alpha.*(BestX-PopPos(i,:)); %Equation (1)
       end         
    end

    for i=1:nPop       
      newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
      newPopFit(i)=fun(newPopPos(i, :), L, Rs, Re, data);  
      if newPopFit(i)>PopFit(i)
         PopFit(i)=newPopFit(i);
         PopPos(i,:)=newPopPos(i,:);
      end
    end

       
  %---------------------------------        
     
    S=2;
    for i=1:nPop   
        newPopPos(i,:)=newPopPos(i,:)+S*(rand*BestX-rand*newPopPos(i,:)); %Equation (8)
    end

     
    for i=1:nPop    
        newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
        newPopFit(i)=fun(newPopPos(i, :), L, Rs, Re, data);   
        if newPopFit(i)>PopFit(i)
           PopFit(i)=newPopFit(i);
           PopPos(i,:)=newPopPos(i,:);
        end
    end
     

     
%       for i=1:nPop
%          if PopFit(i)<Leader_score
%             Leader_score=PopFit(i);
%             BestX=PopPos(i,:);            
%          end
%       end
%     Curve(It) = Leader_score;
  

    % 初始化领导者的位置和适应度值
    [bestfitness, bestindex] = max(PopFit);
    if bestfitness > Leader_score
        BestX = PopPos(bestindex, :);
        Leader_score = bestfitness;
    end
    Curve(It) = Leader_score;
    
    disp(['MRFO_CS: At iteration ', num2str(It), ' ,the best fitness is ', num2str(Leader_score)]);
end

Bestfitness = Curve(end);
%% 结果显示
x = BestX(1:2:end);
y = BestX(2:2:end);
disp('最优位置：' );
for i = 1:Dim/2
    disp([num2str(x(i)), '     ', num2str(y(i))]);
end
disp(['最优覆盖率：', num2str(Leader_score)]);
%% 绘图
% figure;
% plot(Curve, 'r', 'lineWidth', 2);          %  画出迭代图
% title('算法训练过程', 'fontsize', 12);
% xlabel('迭代次数', 'fontsize', 12);
% ylabel('粒子覆盖率', 'fontsize', 12);

figure
for i = 1:Dim/2
    axis([0 ub 0 ub]);            % 限制坐标范围
    sita = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    %fill(x(i)+Rs*cos(sita), y(i)+Rs*sin(sita), 'k');
    fill(x(i)+Rs*cos(sita), y(i)+Rs*sin(sita), 'w');
end
plot(x, y, 'r+');
str_title= ['AMRFOCS coverage rate：',num2str(Bestfitness)];
title({str_title},'FontSize',12);



