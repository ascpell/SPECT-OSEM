clear all;
%读取数据
sys_m = load('sys_m.mat');
sys_m = sys_m.sys_m;

a=fopen('Proj_1e5Counts');
p=fread(a,'float');
proj=reshape(p,[128,60]);
fclose(a);

N=128;%图像大小
ds = 4.0625;%像素宽度
ray_num = 128;%射线数
theta = 0:6:354;%投影对应的角度
img = ones(N*N,1);

iter =0;
os_num = 6;%有序子集数
iternum = 10;%完成60次MLEM迭代需要迭代次数，方便比较不同有序子集数目下图像的变化
count = 1;%有序子集编号，从count开始，每隔iternum个取一个投影作为有序子集

background_mean= zeros(4*iternum+1,1);
std = zeros(4*iternum+1,1);%背景区域标准差
%背景区域CRC及CNR
CNR1 = zeros(4*iternum+1,1);
CRC1 = zeros(4*iternum+1,1);

CNR2 = zeros(4*iternum+1,1);
CRC2 = zeros(4*iternum+1,1);

CNR3 = zeros(4*iternum+1,1);
CRC3 = zeros(4*iternum+1,1);

while(iter < 4*iternum+1)

    if(count>os_num)
        count = 1;%每次迭代都换一组有序子集，当子集取完从头开始取
    end
    %选取有序子集，每隔iternum*个角度取一个以使得有序子集分布均匀
    proj_os = proj(:,count:os_num:os_num*(iternum-1)+count);
    proj_os = reshape(proj_os,ray_num*iternum,1);
    
    %选取有序子集对应的系统矩阵
    index = [];
    for os_idx = count:os_num:count+(iternum-1)*os_num
        index = [index,(os_idx-1)*ray_num+1:os_idx*ray_num];
    end
    sys_os = sys_m(:,index);
    
    %OSEM迭代
    i = 1:N*N;
    temp = img.*sys_os; 
    s = sum(temp)';
    ratio = proj_os./sum(temp)';
    k = sum(sys_os(i,:).*(ratio'),2);

    m = sum(sys_os(i,:),2); 
    m(m==0)=4.0625*iternum/2;%对于某个像素i，有可能所有射线都与这个像素无相交，为使迭代继续设定一个平均值
    img(i) = img(i)./m .* k;%更新图像
    
    %所选取的三个热区的圆心及半径，以像素为单位
    roi1 = [N/2-50/ds,N/2+86.6/ds];
    R1 = 48/ds/2;
    roi2 = [N/2+50/ds,N/2+86.6/ds];
    R2 = 40/ds/2;
    roi3 = [N/2+100/ds,N/2];
    R3 = 32/ds/2;
    
    %取背景区域记为img_bc
    img1 = reshape(img,N,N);
    img_bc = img1(54:74,54:74);
    
    %计数背景区域均值及方差
    background_mean(iter+1) = sum(sum(img_bc))/40/15;
    std(iter+1) = std2(img_bc);
    
    %计数CRC及CNR
    sum1 = 0;
    c = 0;
    for y = round(N/2-50/ds-R1+1):round(N/2-50/ds+R1)
        for x = round(N/2+86.6/ds-R1+1):round(N/2+86.6/ds+R1)
            if((y-roi1(1))^2+(x-roi1(2))^2 <= R1^2)
                sum1 = sum1 + img1(y,x);
                c = c+1;
            end
        end
    end
    mean_roi1 = sum1/c;
    CRC1(iter+1) = abs(mean_roi1-background_mean(iter+1))/background_mean(iter+1)/4;
    CNR1(iter+1) = abs(mean_roi1-background_mean(iter+1))/std(iter+1);

    sum1 = 0;
    c = 0;
    for y = round(roi2(1)-R2+1):round(roi2(1)+R2)
        for x = round(roi2(2)-R2+1):round(roi2(2)+R2)
            if((y-roi2(1))^2+(x-roi2(2))^2 <= R2^2)
                sum1 = sum1 + img1(y,x);
                c = c+1;
            end
        end
    end
    mean_roi2 = sum1/c;
    CRC2(iter+1) = abs(mean_roi2-background_mean(iter+1))/background_mean(iter+1)/4;
    CNR2(iter+1) = abs(mean_roi2-background_mean(iter+1))/std(iter+1);

    sum1 = 0;
     c = 0;
    for y = round(roi3(1)-R3+1):round(roi3(1)+R3)
        for x = round(roi3(2)-R3+1):round(roi3(2)+R3)
            if((y-roi3(1))^2+(x-roi3(2))^2 <= R3^2)
                sum1 = sum1 + img1(y,x);
                 c = c+1;
            end
        end
    end
    mean_roi3 = sum1/c;
    CRC3(iter+1) = abs(mean_roi3-background_mean(iter+1))/background_mean(iter+1)/4;
    CNR3(iter+1) = abs(mean_roi3-background_mean(iter+1))/std(iter+1);
    
    count = count+1;
    iter = iter + 1;
end

img1 = reshape(img,N,N);
img1 = rot90(img1,2);
imshow(img1/max(max(img1)));


figure;
plot(CRC1);
hold on
plot(CRC2);
hold on
plot(CRC3);
legend('CRC1','CRC2','CRC3');
title('1e6CRC');

figure;
plot(CNR1);
hold on
plot(CNR2);
hold on
plot(CNR3);
legend('CNR2','CNR2','CNR3');
title('1e6CNR');
%保存数据
a = fopen('1e5NoCor','a');
fwrite(a,2*img1);
fclose(a);