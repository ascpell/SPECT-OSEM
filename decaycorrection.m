%{
程序由系统矩阵计算衰减矫正后的系统矩阵
1、首先计算背景衰减系数分布，再遍历每个像素
2、再遍历穿过这个像素的射线
3、计算射线穿过此像素后到探测器距离内的衰减系数积分
4、最后赋值给衰减矫正矩阵
%}
%读取系统矩阵
sys_m = load('sys_m.mat');
sys_m = sys_m.sys_m;

N=128;%图像大小
ds = 4.0625;%像素尺寸
ray_num = 128;%射线数目

decayCor_mat = zeros(size(sys_m,1),size(sys_m,2));%衰减矫正矩阵

mu_water = double(0.01538);%unit /mm
mu = zeros(N,N);
for i = 1:N
    for j = 1:N
        x = i-(N+1)/2;
        y = j-(N+1)/2;
        d = sqrt((x*ds)^2+(y*ds)^2);
        if(d<200)
            mu(i,j) = mu_water;%对于半径200mm内的像素赋值为mu_water
        end
    end
end
theta = 0:6:354;%每个投影对应角度

for m=1:N
    for n = 1:N
        x = N-m+1;%系统矩阵的重建结果与实际图像水平方向需要翻转
        y = n;

        index = (y-1)*N+x;%图像中坐标(x,y)在系统矩阵中的编号

        temp = find(sys_m(index,:) > 0);
        ray_crossxy = sys_m(:,temp);%取穿过(x,y)这点的射线，即相交长度不为0的射线

        for ray = temp%遍历这些射线
            deg = floor(ray/ray_num)*6;%射线对应角度
            ray_idx = ray-(deg/6)*ray_num;%在对应角度下射线编号

            dis = 0;
            %不同角度到探测器的方向不一样
            if(deg >= 0 && deg<90)
                d =sys_m(:,ray);
                idx=find(d>0);%只取该射线穿过的像素

                for id = idx'
                    ix = mod(id,N);
                    ix(ix == 0) = N;
                    iy = ceil((id)/N);
                    
                    %0-90°下，射线在(x,y)右上方的像素到达探测器
                    if((iy== y) && (ix >= x))
                        lm = sys_m(id,ray);
                        dis= dis+sys_m(id,ray)*mu(iy,ix);
                    elseif(iy<y)
                        dis = dis + sys_m(id,ray)*mu(iy,ix);
                    end
                end
            elseif(deg>90 && deg<=180)
                d =sys_m(:,ray);
                idx=find(d>0);%只取该射线穿过的像素
                for id = idx'
                    ix = mod(id,N);
                    ix(ix == 0) = N;
                    iy = ceil((id)/N);
                    
                    %90-180°下，射线在(x,y)左上方的像素到达探测器
                    if(iy== y && ix <= x)
                        dis= dis+sys_m(id,ray)*mu(iy,ix);
                    elseif(iy<y)
                        dis = dis + sys_m(id,ray)*mu(iy,ix);
                    end
                end
            elseif(deg>180&&deg<270)
                d =sys_m(:,ray);
                idx=find(d>0);%只取该射线穿过的像素
                for id = idx'
                    ix = mod(id,N);
                    ix(ix == 0) = N;
                    iy = ceil((id)/N);
                    
                    %180-270°下，射线在(x,y)左下方的像素到达探测器
                    if(iy== y && ix <= x)
                        dis= dis+sys_m(id,ray)*mu(iy,ix);
                    elseif(iy>y)
                        dis = dis + sys_m(id,ray)*mu(iy,ix);
                    end
                end
            else
                d =sys_m(:,ray);
                idx=find(d>0);%只取该射线穿过的像素
                for id = idx'
                    ix = mod(id,N);
                    ix(ix == 0) = N;
                    iy = ceil((id)/N);
                    
                    %180-270°下，射线在(x,y)右下方的像素到达探测器
                    if(iy== y && ix >= x)
                        dis= dis+sys_m(id,ray)*mu(iy,ix);
                    elseif(iy>y)
                        dis = dis + sys_m(id,ray)*mu(iy,ix);
                    end
                end
            end

            decayCor_mat(index,ray) = sys_m(index,ray)*exp(-dis);
        end
        
    end
end
save('decayCor.mat','decayCor_mat');