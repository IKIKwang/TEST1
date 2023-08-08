%%主程序
clc;
clear ;
format long;
close all;




global N;
global M;
global L;
global W;
global Grid_cen_x;
global Grid_cen_y;
global Grid_cen_x_and_y;
global sur_area_grid;
global ger;
global r;
global R;
L = 10;% 曲面的长
W = 10;%曲面的宽
%假设1平方米一个网格
M = (L/0.5)*(W/0.5);%网格总数

r = 1; %先一种感知半径
R = 2*r;%通信半径
%假设大、中为5，剩下为小
N = 100;%30个传感器节点
sizepop = 50;%种群规模
dimension = 3;% 空间维数x、y，z
ger = 300;% 最大迭代次数
Up=1;
Low=-1;
limit = [0, 10];            % 设置位置参数限制
%个数限制

load fix_serson_public.mat;
fix_serson = fix_serson_public;
fix_serson_index = (N-6+1:1:N);%共6个

%实验次数
Exper_num = 1;
Exper_best_fitness_ger = zeros(Exper_num,ger);%保存20次实验中每代最好的值

for Exper_index = 1:Exper_num

    load wolf_pop_public.mat;%加载该种群
    load first_init_public.mat;%加载第一个个初始方案
    wolf_pop = wolf_pop_public;


    wolf_pop_temp = zeros(dimension,N,sizepop);%临时种群





    %画个曲面先
    figure(1);
    x_range_point = (0:0.2:10);%配合画曲面
    Dreama_an_grid(x_range_point);
    title('马鞍形曲面山坡');
    hold on;


    %画下初始化部署的三维图  可查看二维
    figure(2);
    x_range_point = (0:0.2:10);%配合画曲面
    Dreama_an_grid(x_range_point);
    hold on;
    %画第一种部署方案
    for i=1:N
        draw_ball(r,wolf_pop(1,i,1),wolf_pop(2,i,1),wolf_pop(3,i,1));
        hold on;
    end


    %查看三维图的
    figure(3);
    x_range_point = (0:0.2:10);%配合画曲面
    Dreama_an_grid(x_range_point);
    hold on;
    %画第一种部署方案
    for i=1:N
        plot3(wolf_pop(1,i,1),wolf_pop(2,i,1),wolf_pop(3,i,1),'.','MarkerSize',20,'color','r');%画球心 
        hold on;
    end


    %接下来二维面划分2500个网格
    %求网格中心坐标
    X_mat = (0:0.5:10);%x矩阵
    Y_mat = (0:0.5:10);%y矩阵


    %计算每个网格的表面积
    sur_area_grid = zeros(L/0.5,W/0.5);%用于存每个小曲面的表面积 最底下的为第一行
    for i=1:L/0.5
        for j=1:W/0.5
            sur_area_grid(i,j) = get_sur_area(X_mat(1,j),X_mat(1,j+1),Y_mat(1,i),Y_mat(1,i+1));
        end
    end

    all_sur_area = sum(sum(sur_area_grid));


    Grid_cen_x = zeros(1,L/0.5);%网格中心点x坐标
    Grid_cen_y = zeros(1,W/0.5);%网格中心点y坐标
    %前后两者相加之和除以2
    for i=1:L/0.5
        Grid_cen_x(1,i) = (X_mat(i)+X_mat(i+1))/2;
        Grid_cen_y(1,i) = (Y_mat(i)+Y_mat(i+1))/2;
    end

    %监测点
    %最下面那一行为第一行
    Grid_cen_x_and_y = zeros(L,W,3);%共L*W个网格中心，但是每个网格中心有x,y,z
    for i=1:L/0.5
        for j=1:W/0.5
            Grid_cen_x_and_y(i,j,1) = Grid_cen_x(1,j);%1代表x
            Grid_cen_x_and_y(i,j,2) = Grid_cen_y(1,i);%把y坐标放到第二行
            Grid_cen_x_and_y(i,j,3) = ma_an([Grid_cen_x_and_y(i,j,1),Grid_cen_x_and_y(i,j,2)]);%把z坐标放到第三行
        end
    end

    %计算种群的覆盖率
    wolf_fitness = zeros(1,sizepop);

    best_fitness = -inf;%求最大的
    record_ger = zeros(1,ger);%每次迭代中最好的fitness
    record_pop_ave = zeros(1,ger);
    wolf_leader_pops = zeros(dimension,N,sizepop);%整合成三维

    %灰狼算法更新
    iter = 1;
    while iter <= ger

        %计算该种群中每个个体的适应值

        %计算种群覆盖率
        for i=1:sizepop
            wolf_fitness(1,i) = get_Grid_cover_unit_and_rate(wolf_pop(:,:,i),all_sur_area);
        end

        %保存初始化对应的覆盖率
        if iter == 1
            init_fitness = wolf_fitness(1,1);
        end
        
        [fitness_max,max_index] = max(wolf_fitness);%对适应值进行排序
        wolf_leader = wolf_pop(:,:,max_index);%领头狼

        for i=1:sizepop
            wolf_leader_pops(:,:,i) = wolf_leader;%化成三维
        end


        %更新最好的的适应值
        if fitness_max > best_fitness
            best_fitness = fitness_max;
            best_indivi = wolf_leader;
        end
        record_ger(1,iter) = best_fitness;%保存该代中最优的适应度值

        disp(['实验次数：',num2str(Exper_index)]);
        disp(['迭代次数：',num2str(iter)]);
        disp(['最佳覆盖率：',num2str(best_fitness)]);

        sum_fitness = sum(wolf_fitness);
        record_pop_ave(1,iter) = sum_fitness/sizepop;%算出该种群的平均适应度

        
        Coef=iter/ger;
        omega = 3;
        a = 2*(1-Coef);
        P = randn*((sin(pi/2 * Coef))^omega + cos(pi/2 * Coef)-1);  % 动态扰动因子策略
        M = a*(2*rand-1) + P;     
        if rand<0.5
          r1=rand;                         
          Beta=2*exp(r1*((ger-iter+1)/ger))*(sin(2*pi*r1));    
          if  abs(M)<1   
              pop_new = wolf_leader_pops+rand(dimension,N,sizepop).*(wolf_leader_pops-wolf_pop(:,1))+Beta*(wolf_leader_pops-wolf_pop(:,1));
%               pop_new(j,k,i) = wolf_leader_pops(j,k,i);
%               newPopPos(1,:)=BestX+rand(1,Dim).*(BestX-PopPos(1,:))+Beta*(BestX-PopPos(1,:)); %Equation (4)
          else              
             % IndivRand=rand(1,Dim).*(Up-Low)+Low;   
              IndivRand=get_cuckoos(wolf_pop(:,1),wolf_leader_pops,Low,Up);
              pop_new = IndivRand+rand(dimension,N,sizepop).*(IndivRand-wolf_pop(:,1))+Beta*(IndivRand-wolf_pop(:,1));
%              newPopPos(1,:)=IndivRand+rand(1,Dim).*(IndivRand-PopPos(1,:))+Beta*(IndivRand-PopPos(1,:)); %Equation (7)         
          end              
        else 
            Alpha=2*rand(dimension,N,sizepop).*(-log(rand(dimension,N,sizepop))).^0.5;    
            pop_new = wolf_pop(:,1)+rand(dimension,N,sizepop).*(wolf_leader_pops-wolf_pop(:,1))+Alpha.*(wolf_leader_pops-wolf_pop(:,1));
  %          newPopPos(1,:)=PopPos(1,:)+rand(1,Dim).*(BestX-PopPos(1,:))+Alpha.*(BestX-PopPos(1,:)); %Equation (1)
        end

     for i=2:sizepop
        if abs(M)<1
             r1=rand;                         
             Beta=2*exp(r1*((ger-iter+1)/ger))*(sin(2*pi*r1));    
             if  Coef>rand  
                 pop_new = wolf_leader_pops+rand(dimension,N,sizepop).*(wolf_pop(:,i-1)-wolf_pop(:,i))+Beta*(wolf_leader_pops-wolf_pop(:,i));
                % newPopPos(i,:)=BestX+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(BestX-PopPos(i,:)); %Equation (4)
             else
                 %IndivRand=rand(1,Dim).*(Up-Low)+Low;   
                 IndivRand=get_cuckoos(wolf_pop(:,i-1),wolf_leader_pops,Low,Up);
                 pop_new = IndivRand+rand(dimension,N,sizepop).*(wolf_pop(:,i-1)-wolf_pop(:,i))+Beta*(IndivRand-wolf_pop(:,i));
               %  newPopPos(i,:)=IndivRand+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(IndivRand-PopPos(i,:));  %Equation (7)       
             end              
        else
            Alpha=2*rand(dimension,N,sizepop).*(-log(rand(dimension,N,sizepop))).^0.5;     
            pop_new = wolf_pop(:,i)+rand(dimension,N,sizepop).*(wolf_pop(:,i-1)-wolf_pop(:,i))+Alpha.*(wolf_leader_pops-wolf_pop(:,i));
           % newPopPos(i,:)=PopPos(i,:)+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Alpha.*(BestX-PopPos(i,:)); %Equation (1)
       end         
     end
    
%      S=2;
%     for i=1:sizepop    
%         pop_new = wolf_pop(:,i)+S*(rand*(wolf_leader_pops-wolf_pop(:,i)));
%         %newPopPos(i,:)=newPopPos(i,:)+S*(rand*BestX-rand*newPopPos(i,:)); %Equation (8)
%     end

    
    wolf_pop = wolf_leader_pops;
    
         %让两个节点变化位置
        swap_node_num = 1;
        for i=1:sizepop
            node_index = randperm(N,swap_node_num);
            for k=1:swap_node_num
                wolf_pop(:,node_index(1,k),i) = pop_new(:,node_index(1,k),i);
                wolf_pop(3,node_index(1,k),i) = ma_an([wolf_pop(1,node_index(1,k),i),wolf_pop(2,node_index(1,k),i)]);%得更新第三维
            end
            
            %把固定维度加上
%             for z = 1:6
%                 wolf_pop(:,fix_serson_index(1,z),i) = fix_serson(:,z);
%             end
        end


       


        iter = iter+1;
    end
    clc;
    %保存数据
    Exper_best_fitness_ger(Exper_index,:) = record_ger;
end


%最后求平均值 然后保存
HGP_best_fitness_ger = sum(Exper_best_fitness_ger,1)/Exper_num;
save HGP_best_fitness_ger.mat HGP_best_fitness_ger;


figure(6);
plot(record_ger);
title('每代中最优的适应值')

figure(7);
plot(record_pop_ave);
title('每个种群的平均适应值')


%最后的部署图
%画下初始化部署的三维图  可查看二维
figure(8);
x_range_point = (0:0.2:10);%配合画曲面
Dreama_an_grid(x_range_point);
hold on;
%画第一种部署方案
for i=1:N
    draw_ball(r,best_indivi(1,i),best_indivi(2,i),best_indivi(3,i));
    hold on;
end


%查看三维图的
figure(9);
x_range_point = (0:0.2:10);%配合画曲面
Dreama_an_grid(x_range_point);
hold on;
%画第一种部署方案
for i=1:N
    plot3(best_indivi(1,i),best_indivi(2,i),best_indivi(3,i),'.','MarkerSize',20,'color','r');%画球心 
    hold on;
end



disp('进入连通性的判断...');
pause(1);
%连通性的判断
is_connec = get_3D_connection(best_indivi);
if is_connec == 1
    disp('整个网络是连通的');
else
    disp('整个网络不是连通的');
end





%最小生出树的产生
figure(10);
x_range_point = (0:0.2:10);%配合画曲面
Dreama_an_grid(x_range_point);
hold on;
axis([0,10,0,10,-1,1]);%横纵坐标最小和最大值 
hold on;
for i=1:N
    plot3(best_indivi(1,i),best_indivi(2,i),best_indivi(3,i),'.','MarkerSize',20,'color','r');%画球心 
    hold on;
end
grid on
draw_MST(best_indivi(1,:),best_indivi(2,:),best_indivi(3,:));
hold on;



figure(11);
axis([0,10,0,10,-1,1]);%横纵坐标最小和最大值 
hold on;
for i=1:N
    plot3(best_indivi(1,i),best_indivi(2,i),best_indivi(3,i),'.','MarkerSize',20,'color','r');%画球心 
    hold on;
end
grid on
draw_MST(best_indivi(1,:),best_indivi(2,:),best_indivi(3,:));
hold on;



%绘初始化的覆盖图
figure(12);
Grid_cover_unit = draw_init_end_cover(first_init_public);
x_range_point = (0:0.2:10);%配合画曲面
Dreama_an_grid(x_range_point);
hold on;
axis([0,10,0,10,-1,1]);%横纵坐标最小和最大值 
hold on;
for i=1:N
    plot3(first_init_public(1,i),first_init_public(2,i),first_init_public(3,i),'.','MarkerSize',20,'color','r');%画球心 
    hold on;
end
grid on
%画覆盖的区域
for i=1:L/0.5
    for j=1:W/0.5
        if Grid_cover_unit(i,j) == 1 %该监测点被覆盖
            plot3(Grid_cen_x_and_y(i,j,1),Grid_cen_x_and_y(i,j,2),Grid_cen_x_and_y(i,j,3),'.','MarkerSize',20,'color','k');%画球心
            hold on;
        end
    end
end
hold on;



%绘最终的覆盖图
figure(13);
Grid_cover_unit = draw_init_end_cover(best_indivi);
x_range_point = (0:0.2:10);%配合画曲面
Dreama_an_grid(x_range_point);
hold on;
axis([0,10,0,10,-1,1]);%横纵坐标最小和最大值 
hold on;
for i=1:N
    plot3(best_indivi(1,i),best_indivi(2,i),best_indivi(3,i),'.','MarkerSize',20,'color','r');%画球心 
    hold on;
end
grid on
%画覆盖的区域
for i=1:L/0.5
    for j=1:W/0.5
        if Grid_cover_unit(i,j) == 1 %该监测点被覆盖
            plot3(Grid_cen_x_and_y(i,j,1),Grid_cen_x_and_y(i,j,2),Grid_cen_x_and_y(i,j,3),'.','MarkerSize',20,'color','k');%画球心
            hold on;
        end
    end
end
hold on;

save best_indivi.mat best_indivi;

%调用最后一个种群连通的判断函数
last_pop_connec(wolf_pop,sizepop,all_sur_area);
