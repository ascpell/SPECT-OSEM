%{
计算系统矩阵采用的是siddon方法
程序大致结构为：
1、将图像视为网格，首先得到网格每条线的x坐标与y坐标
2、遍历每个角度每条射线，根据射线的直线方程得到与网格y方向直线的交点 x_cross
3、将x_cross中超出图像范围的值滤除，再滤除与x相等的点
4、将x与x_cross组合成intersect_x起来再进行排序
5、计算射线上intersect_x对应的intersect_y的值
6、遍历以上的点，计算每两点之间的距离及穿过的像素编号，赋值给system matrix

%}

clear all;
N = 128;    %size of img (unit pixel)
size = 520;  %size of img (unit mm)

theta = 0:6:354;  %degree of each projection
ds = 4.0625;  %interval of detector cell(unit pixel)

ray_num = 128;    %cell number of detector/ray

[sys_m] = proj(ray_num,theta,ds,N);
save('sys_m.mat','sys_m');

function [sys_m] = proj(ray_num,theta,ds,N)

sys_m = zeros(N*N,ray_num*size(theta,2));%system matrix

for i = theta %for every degree
    beta = (i)*pi/180;%convert degree to rad 
    index = floor(i/6)+1;%proj index of projection in i
    
    %sin(i)=0 while i = 180 and i = 0
    if(i == 180)
        for j = 1:ray_num
            n = (1:N)+(ray_num-j)*N;%pixel index the j-th ray cross
            sys_m(n,(index-1)*(ray_num)+j) = ds*ones(N,1);
        end
    elseif(i == 0)
        for j = 1:ray_num
            n = (1:N)+(ray_num-j)*N;%pixel index the j-th ray cross
            sys_m(n,(index-1)*(ray_num)+j) = ds*ones(N,1);
            projection(index,j) = sum(ds*img(:,ray_num-j+1));
        end
    else%for other degree
        for j = 1:ray_num
            rad = beta;
            t = (j-(ray_num+1)/2)*ds;%index of each detector cell while the center of detector is 0
            
            %the grid's x and y 
            x = (-N/2:N/2)*ds; 
            y = (-N/2:N/2)*ds;
    
            x_cross = y*tan(rad) + t / cos(rad);%the intersection's x of j-th ray to grid's y
            y_cross = x/tan(rad) - t / sin(rad);%the intersection's y of j-th ray to grid's x
    
            x_cross(abs(x_cross) > N/2*ds) = N*100; % note the point exceed img as 100*N
            y_cross(abs(y_cross) > N/2*ds) = N*100; 
    
            x(y_cross == N*100) = N*100;% in fact we only need x_cross and x
            
            for k = 1:size(x,2) %remove same point
                for h = 1:size(x_cross,2)
                    if(x(k) == x_cross(h))
                        x_cross(h) = N*100;% note the same point as 100*N
                    end
                end
            end
            
            p = zeros(1,2*(N+1));%unite x and x_cross into one vector
            p(1,1:N+1) = x;
            p(1,N+2:2*(N+1)) = x_cross;
    
            intersect_x = p(p < N*100);%filter noted point
            intersect_x = sort(intersect_x,2);%sort the point
           
            intersect_y = intersect_x/tan(rad) - t / sin(rad);%calculate y through x value

            for m = 1:size(intersect_x,2)-1
                current_x = intersect_x(m);
                current_y = intersect_y(m);
    
                next_x = intersect_x(m+1);
                next_y = intersect_y(m+1);
    
                posy = ceil(next_y/ds)+N/2;%the index of the pixel between two point
                posx = ceil(next_x/ds)+N/2;
                
                
                if(posy > N)
                    posy = N;
                end
                if posy <= 0
                    posy = 1;
                end

                if(posx > N)
                    posx = N;
                end
                if posx <= 0
                    posx = 1;
                end
                
                %calculate system matrix
                sys_m((posx-1)*N+posy,(index-1)*ray_num + j) = sys_m((posx-1)*N+posy,(index-1)*ray_num + j)+...
                double(sqrt((next_x - current_x)^2 + (next_y-current_y)^2));    
                
            end
        end
    end
end

end