%% 视频预处理 
hsrc = vision.VideoFileReader('MVI_8229.avi', ...
                                  'VideoOutputDataType', 'uint8');
    hfg = vision.ForegroundDetector(...
        'NumTrainingFrames', 150, ... % 5 because of short video
        'InitialVariance', (30)^2,...
        'NumGaussians',4,...
        'MinimumBackgroundRatio',0.7); % initial standard deviation of 30
    hfg1 = vision.ForegroundDetector(...
        'NumTrainingFrames', 150, ... % 5 because of short video
        'InitialVariance', (30)^2,...
        'NumGaussians',4,...
        'MinimumBackgroundRatio',0.7); % initial standard deviation of 30
    hblob = vision.BlobAnalysis(...
        'CentroidOutputPort', false, 'AreaOutputPort', false, ...%不输出质心坐标，不输出面积
        'BoundingBoxOutputPort', true, 'MinimumBlobArea', 500);
    hsi = vision.ShapeInserter('BorderColor','White');
  ti=20;
  n=1;
%     hsnk = vision.VideoPlayer();
    while ti~=1%~isDone(hsrc)
      step(hsrc);%获得视频中的下一帧
       step(hsrc);%获得视频中的下一帧
      frame  = step(hsrc);%获得视频中的下一帧
     %%
    ti=ti-1;
    if(ti==0)%短期模型       
        reset(hfg1);
        ti=20;
    end

     %% 
      fgMask = step(hfg, frame); %获得前景的mask，格式是布尔量
      fgMask1 = step(hfg1, frame);
     % bbox= step(hblob, fgMask);
      mask=fgMask.*ones(480,640)*255;
      mask1=fgMask1.*ones(480,640)*255;
      maskr=frame;
      maskr(:,:,1)=mask;
      maskr(:,:,2)=mask1;
      maskr(:,:,3)=zeros(480,640);
       if(ti==1)%
        PR(n,:,:)=mask;
        PG(n,:,:)=mask1;
%         n=n+1;
        break
%      out   = frame+maskr; % draw
%       step(hsnk, out); % view results in the video player
    end
    
    end
%     release(hsnk);
    release(hsrc);
%% 视频处理结束 
PG2=zeros(480,640);
PG2(:,:)=PG(1,:,:).*PR(n,:,:);%结合短期学习与长期学习获得的前景
x=-90:640;
%% 录入图像并显示
f=PG2;
figure(1);
subplot(2,2,1);imshow(frame+maskr);title('原图');

%% 提取图像边缘
[m,n]=size(f);%得到图像矩阵行数m，列数n
for i=3:m-2
    for j=3:n-2%处理领域较大，所以从图像（3,3）开始，在（m-2,n-2）结束 
        l(i,j)=-f(i-2,j)-f(i-1,j-1)-2*f(i-1,j)-f(i-1,j+1)-f(i,j-2)-2*f(i,j-1)+16*f(i,j)-2*f(i,j+1)-f(i,j+2)-f(i+1,j-1)-2*f(i+1,j)-f(i+1,j+1)-f(i+2,j);%LoG算子
    end
end
subplot(2,2,2);imshow(l);title('LoG算子提取图像边缘');

%%  滤波
[m,n]=size(l);
for i=2:m-1    
 for j=2:n-1   
    y(i,j)=l(i-1,j-1)+l(i-1,j)+l(i-1,j+1)+l(i,j-1)+l(i,j)+l(i,j+1)+l(i+1,j-1)+l(i+1,j)+l(i+1,j+1);
    y(i,j)=y(i,j)/9;  %LoG算子提取边缘后，对结果进行均值滤波以去除噪声，为下一步hough变换提取直线作准备
 end
end
subplot(2,2,3);imshow(y);title('均值滤波器处理后')

%% 二值化
q=im2uint8(y);
[m,n]=size(q);
for i=1:m    
 for j=1:n  
    if q(i,j)>80;   %设置二值化的阈值为80
        q(i,j)=255; %对图像进行二值化处理，使图像边缘更加突出清晰
    else
        q(i,j)=0; 
    end
 end
end
subplot(2,2,4);imshow(q);title('二值化处理后');

%% 检测直线
 
%Hough变换检测直线，使用（a，p）参数空间，a∈[0,180],p∈[0,2d]
a=180; %角度的值为0到180度
d=round(sqrt(m^2+n^2)); %图像对角线长度为p的最大值
s=zeros(a,2*d); %存储每个(a,p)个数
[aa,bb]=find(q==255);
z=zeros(1,2);   %bu用元胞存储每个被检测的点的坐标
zz=zeros(1,2);
zzz=zeros(1,2);%用于最后的简化
ll=length(aa);%ll个
for loop=1:ll
        i=aa(loop);
     j=bb(loop);
            for k=1:a
                p = round(i*cos(pi*k/180)+j*sin(pi*k/180));%对每个点从1到180度遍历一遍，取得经过该点的所有直线的p值(距离)（取整）
                if(p > 0)%若p大于0，则将点存储在（d，2d）空间
                    s(k,d+p)=s(k,d+p)+1;%（a，p）相应的累加器单元加一
                    if(s(k,d+p) == 125)
                       z=[z;k,p];
                    end
                else
                    ap=abs(p)+1;%若p小于0，则将点存储在（0，d）空间
                    s(k,ap)=s(k,ap)+1;
                   if(s(k,ap) ==125)
                    z=[z;k-180,ap]; 
                    end 
                end
            end
       
end


%% 显示效果
figure(22);                                                        
plot(0,0);
axis([0 640 0 480]);
title('直线检测');
hold;

% xx=z(:,2)./sin((pi.*z(:,1)./180));
yy=z(:,2)./sin(pi.*z(:,1)./180);
k=-cot(pi.*z(:,1)./180);

a=find(z(:,1)~=0);
for i=2:(1+length(a))    
 y=k(i).*x+yy(i);
 plot(y,480-x,'-r','LineWidth',1); %%%%%%%%%%%%%%%55
  D(i)=y(1);
end

figure(201)
for i=-800:800
    DD(i+801)=sum(D(D>((i-1)*10+1)&D<i*10));
end

bar(-800:800,DD)
title('按视平线上的截距统计');
xlabel('视平线截距')
ylabel('直线数目统计')

%%
figure(202)
plot(z(:,1),z(:,2),'r+')
hold on
i=find(DD==max(DD(801:801+70)))-801;%
zz=[z((D>((i-1)*10+1)&D<i*10),1),z((D>((i-1)*10+1)&D<i*10),2)];
plot(zz(:,1),zz(:,2),'go')%圈出消失点在前方的复合条件的点
axis([-30 180 0 480]);
title('符合条件直线簇特征值');
xlabel('角度')
xlabel('距离')
%% 识别结果
figure(203)
a=find(zz(:,1)~=0);
xx=zz(:,2)./sin((pi.*zz(:,1)./180));
yy=zz(:,2)./sin(pi.*zz(:,1)./180);
k=-cot(pi.*zz(:,1)./180);
title('筛选后的直线');
imshow(frame);hold on;
for i=1:(length(a))    
 y=k(i).*x+yy(i);
 plot(y,x,'g')
end
axis([0 640 0 480]);
%%