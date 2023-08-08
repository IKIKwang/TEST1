function [Bestfitness, Leader_pos, Curve] = MRFO(nPop, MaxIt, Rs, Re, L, data, PopPos, PopFit, Low, Up, Dim)

%% 
% 初始化节点位置
% for i = 1:N
%         fai = 2*pi*rand;
%     for j = 1:dim
%         sita = pi*rand;
%         Px(i, j) = cos(fai)*sin(sita);
%         Py(i, j) = sin(fai)*cos(sita);
%         Pz(i, j) = cos(sita);
%     end
%     Xx(i, :) = 0.5*(lb*(1+Px(i, :))+ub*(1-Px(i, :)));
%     Xy(i, :) = 0.5*(lb*(1+Py(i, :))+ub*(1-Py(i, :)));
%     Xz(i, :) = 0.5*(lb*(1+Pz(i, :))+ub*(1-Pz(i, :)));
% end
% X = [Xx; Xy; Xz];
% % 计算适应度值
% for i = 1:3*N
%     fitness(i) = fun(X(i, :), L, Rs, Re, data);
% end
% % 选择排名前N的种群
% [~, SortIndex] = sort(fitness, 'descend');
% fitness = fitness(SortIndex(1:N));
% X = X(SortIndex(1:N), :);
% % 初始化领导者的位置和适应度值
% [bestfitness, bestindex] = max(fitness);
% Leader_pos = X(bestindex, :);
% Leader_score = bestfitness; 
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
   % plot(x(i)+Rs*cos(sita), y(i)+Rs*sin(sita), 'b');
end
plot(x, y, 'k+');
str_title= ['Initial deployment coverage rate：',num2str(Leader_score)];
title({str_title},'FontSize',12);

t = 0;       % 迭代次数计数器
%% 迭代寻优
for It=1:MaxIt  
    Coef=It/MaxIt; 

       if rand<0.5
          r1=rand;                         
          Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
          if  Coef>rand                                                      
              newPopPos(1,:)=BestX+rand(1,Dim).*(BestX-PopPos(1,:))+Beta*(BestX-PopPos(1,:)); %Equation (4)
          else
              IndivRand=rand(1,Dim).*(Up-Low)+Low;                                
              newPopPos(1,:)=IndivRand+rand(1,Dim).*(IndivRand-PopPos(1,:))+Beta*(IndivRand-PopPos(1,:)); %Equation (7)         
          end              
       else 
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(1,:)=PopPos(1,:)+rand(1,Dim).*(BestX-PopPos(1,:))+Alpha.*(BestX-PopPos(1,:)); %Equation (1)
       end
     
    for i=2:nPop
        if rand<0.5
           r1=rand;                         
           Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
             if  Coef>rand                                                      
                 newPopPos(i,:)=BestX+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(BestX-PopPos(i,:)); %Equation (4)
             else
                 IndivRand=rand(1,Dim).*(Up-Low)+Low;                                
                 newPopPos(i,:)=IndivRand+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(IndivRand-PopPos(i,:));  %Equation (7)       
             end              
        else
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(i,:)=PopPos(i,:)+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Alpha.*(BestX-PopPos(i,:)); %Equation (1)
       end         
    end

%     for i=1:nPop       
%         newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
%         newPopFit(i)=fun(newPopPos(i, :), L, Rs, Re, data);  
%         if newPopFit(i)>PopFit(i)
%            PopFit(i)=newPopFit(i);
%            PopPos(i,:)=newPopPos(i,:);
%         end
%     end        
    
  

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
%          if PopFit(i)<BestF
%             Leader_score=PopFit(i);
%             BestX=PopPos(i,:);            
%          end
%       end
% 
%       Curve(It) = Leader_score;
  


    % 初始化领导者的位置和适应度值
     [bestfitness, bestindex] = max(PopFit);
     if bestfitness > Leader_score
         BestX = PopPos(bestindex, :);
         Leader_score = bestfitness;
     end
     Curve(It) = Leader_score;
     
    disp(['MRFO: At iteration ', num2str(It), ' ,the best fitness is ', num2str(Leader_score)]);
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
    fill(x(i)+Rs*cos(sita), y(i)+Rs*sin(sita), 'w');
end
plot(x, y, 'r+');
str_title= ['MRFO coverage rate：',num2str(Bestfitness)];
title({str_title},'FontSize',12);



