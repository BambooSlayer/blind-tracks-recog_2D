function varargout = mydemoGUI(varargin)
% MYDEMOGUI MATLAB code for mydemoGUI.fig
%      MYDEMOGUI, by itself, creates a new MYDEMOGUI or raises the existing
%      singleton*.
%
%      H = MYDEMOGUI returns the handle to a new MYDEMOGUI or the handle to
%      the existing singleton*.
%
%      MYDEMOGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYDEMOGUI.M with the given input arguments.
%
%      MYDEMOGUI('Property','Value',...) creates a new MYDEMOGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mydemoGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mydemoGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mydemoGUI

% Last Modified by GUIDE v2.5 12-Jan-2016 18:43:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mydemoGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mydemoGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mydemoGUI is made visible.
function mydemoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mydemoGUI (see VARARGIN)

% Choose default command line output for mydemoGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mydemoGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mydemoGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in READ.
function READ_Callback(hObject, eventdata, handles)
% hObject    handle to READ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% I = imread('文件名.格式后缀');
global f f1 l y q z zz
set(handles.text11,'string','');
%选择图片路径
[filename,pathname,filterindex]=uigetfile({'*.jpg';'*.bmp';'*.gif'},'选择图片');
%合成路径+文件名
str=[pathname filename]; 
 %读取图片
f=imread(str); 
%使用第一个axes 
x=-90:640;
%% 显示图片
axes(handles.axes1)   
imshow(f);%原图
%% 灰度
f1=rgb2gray(f);
f1=histeq(f1);                  %直方图均衡化
axes(handles.axes2);imshow(f1);
f1=im2double(f1);
%% 提取图像边缘
[m,n]=size(f1);%得到图像矩阵行数m，列数n
for i=3:m-2
    for j=3:n-2%处理领域较大，所以从图像（3,3）开始，在（m-2,n-2）结束 
        l(i,j)=-f1(i-2,j)-f1(i-1,j-1)-2*f1(i-1,j)-f1(i-1,j+1)-f1(i,j-2)-2*f1(i,j-1)+16*f1(i,j)-2*f1(i,j+1)-f1(i,j+2)-f1(i+1,j-1)-2*f1(i+1,j)-f1(i+1,j+1)-f1(i+2,j);%LoG算子
    end
end
axes(handles.axes3);imshow(l);%title('LoG算子提取图像边缘');
%%  滤波
[m,n]=size(l);
for i=2:m-1    
 for j=2:n-1   
    y(i,j)=l(i-1,j-1)+l(i-1,j)+l(i-1,j+1)+l(i,j-1)+l(i,j)+l(i,j+1)+l(i+1,j-1)+l(i+1,j)+l(i+1,j+1);
    y(i,j)=y(i,j)/9;  %LoG算子提取边缘后，对结果进行均值滤波以去除噪声，为下一步hough变换提取直线作准备
 end
end
axes(handles.axes4);imshow(y);%title('均值滤波器处理后')
%% 二值化
q=im2uint8(y);
[m,n]=size(q);
for i=1:m    
 for j=1:n  
    if q(i,j)>100;   %设置二值化的阈值为80
        q(i,j)=255; %对图像进行二值化处理，使图像边缘更加突出清晰
    else
        q(i,j)=0; 
    end
 end
end
axes(handles.axes5);imshow(q);%title('二值化处理后');
%% 检测直线
%Hough变换检测直线，使用（a，p）参数空间，a∈[0,180],p∈[0,2d]
a=180; %角度的值为0到180度
d=round(sqrt(m^2+n^2)); %图像对角线长度为p的最大值
s=zeros(a,2*d); %存储每个(a,p)个数
[aa,bb]=find(q==255);
z=zeros(1,2);   %bu用元胞存储每个被检测的点的坐标

ll=length(aa);%ll个
for loop=1:ll
        i=aa(loop);
     j=bb(loop);
            for k=1:a
                p = round(i*cos(pi*k/180)+j*sin(pi*k/180));%对每个点从1到180度遍历一遍，取得经过该点的所有直线的p值(距离)（取整）
                if(p > 0)%若p大于0，则将点存储在（d，2d）空间
                    s(k,d+p)=s(k,d+p)+1;%（a，p）相应的累加器单元加一
                    if(s(k,d+p) == 100)
                       z=[z;k,p];
                    end
                else
                    ap=abs(p)+1;%若p小于0，则将点存储在（0，d）空间
                    s(k,ap)=s(k,ap)+1;
                   if(s(k,ap) ==100)
                    z=[z;k-180,ap]; 
                    end 
                end
            end
       
end
%% 显示效果
axes(handles.axes6);%figure(22);    
hold off;
plot(0,0);
axis([0 640 0 480]);
% title('直线检测');
hold on;

% xx=z(:,2)./sin((pi.*z(:,1)./180));
yy=z(:,2)./sin(pi.*z(:,1)./180);
k=-cot(pi.*z(:,1)./180);

a=find(z(:,1)~=0);
% D=[];
for i=2:(1+length(a))    
 yn=k(i).*x+yy(i);
 plot(yn,480-x,'-r','LineWidth',1); %%%%%%%%%%%%%%%55  z=-90
  D(i)=yn(1);
end
%%
axes(handles.axes7);% figure(201)
hold off;
for i=-800:800
    DD(i+801)=abs(sum(D(D>((i-1)*10+1)&D<i*10)));
end

bar(-800:800,DD)
axis([-200 200 0 max(DD)+10]);
% title('按视平线上的截距统计');
% xlabel('视平线截距')
% ylabel('直线数目统计')
%%
axes(handles.axes8);%figure(202)
hold off
plot(z(:,1),z(:,2),'r+')
hold on
i=find(DD==max(DD(801:801+70)))-801;%
zz=[z((D>((i-1)*10+1)&D<i*10),1),z((D>((i-1)*10+1)&D<i*10),2)];
plot(zz(:,1),zz(:,2),'go')%圈出消失点在前方的复合条件的点
axis([-30 180 0 480]);
% title('符合条件直线簇特征值');
% xlabel('角度')
% ylabel('距离')
%% 识别结果
axes(handles.axes9);%figure(203)
a=find(zz(:,1)~=0);
xx=zz(:,2)./sin((pi.*zz(:,1)./180));
yy=zz(:,2)./sin(pi.*zz(:,1)./180);
k=-cot(pi.*zz(:,1)./180);
% title('筛选后的直线');
imshow(f);hold on;
for i=1:(length(a))    
 yn=k(i).*x+yy(i);
 plot(yn,x,'g')
end
 ii=find((abs(yy-mean(yy))-min(abs(yy-mean(yy))))<15);
for i=1:length(ii)%%max(1,ii-ceil(0.25*length(k))):min(length(k),ii+ceil(0.25*length(k)))  
 yn=k(ii(i)).*x+yy(ii(i));
 plot(yn,x,'r--','linewidth',5)
end

axis([0 640 0 480]);




% --- Executes on button press in yt.
function yt_Callback(hObject, eventdata, handles)
% hObject    handle to yt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f
figure(1);imshow(f);%原图

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f1 
figure(1);imshow(f1);
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global l
figure(1);imshow(l);
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global y 
figure(1);imshow(y);
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global q 
figure(1);imshow(q);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global z zz f
figure(1);
x=-130:640;
%% 显示效果

hold off;


% title('直线检测');

imshow(f);hold on;
axis([-150 640+100 -150 480+100]);
% xx=z(:,2)./sin((pi.*z(:,1)./180));
yy=z(:,2)./sin(pi.*z(:,1)./180);
k=-cot(pi.*z(:,1)./180);

a=find(z(:,1)~=0);
% D=[];
for i=2:(1+length(a))    
 yn=k(i).*x+yy(i);
 plot(yn,x,'g','LineWidth',1); %%%%%%%%%%%%%%%55
  D(i)=yn(1);
end
% plot([0,0,640,640],[0,480,480,0],'b--','LineWidth',3)
% plot([0,640],[0,0],'b--','LineWidth',3)
plot([0,640],[-90,-90],'r--','LineWidth',5)

%%
a=find(zz(:,1)~=0);
yy=zz(:,2)./sin(pi.*zz(:,1)./180);
k=-cot(pi.*zz(:,1)./180);
% title('筛选后的直线');

% for i=1:(length(a))    
%  yn=k(i).*x+yy(i);
%  plot(yn,x,'g')
% end
 ii=find((abs(yy-mean(yy))-min(abs(yy-mean(yy))))<15);
for i=1:length(ii)%%max(1,ii-ceil(0.25*length(k))):min(length(k),ii+ceil(0.25*length(k)))  
 yn=k(ii(i)).*x+yy(ii(i));
 plot(yn,x,'r--','linewidth',5)
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('result.fig');  


% --- Executes on button press in BF.
function BF_Callback(hObject, eventdata, handles)
% hObject    handle to BF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f f1 l y q z zz
x=-90:640;
for prpr=[1 20 40 88 97 103 121]
set(handles.text11,'string',['正在处理第',num2str(prpr),'帧']);
 %读取图片
f=imread([num2str(prpr),'.jpg']); 
%使用第一个axes 

%% 显示图片
axes(handles.axes1)   
imshow(f);%原图
%% 灰度
f1=rgb2gray(f);
f1=histeq(f1);                  %直方图均衡化
axes(handles.axes2);imshow(f1);
f1=im2double(f1);
%% 提取图像边缘
[m,n]=size(f1);%得到图像矩阵行数m，列数n
for i=3:m-2
    for j=3:n-2%处理领域较大，所以从图像（3,3）开始，在（m-2,n-2）结束 
        l(i,j)=-f1(i-2,j)-f1(i-1,j-1)-2*f1(i-1,j)-f1(i-1,j+1)-f1(i,j-2)-2*f1(i,j-1)+16*f1(i,j)-2*f1(i,j+1)-f1(i,j+2)-f1(i+1,j-1)-2*f1(i+1,j)-f1(i+1,j+1)-f1(i+2,j);%LoG算子
    end
end
axes(handles.axes3);imshow(l);%title('LoG算子提取图像边缘');
%%  滤波
[m,n]=size(l);
for i=2:m-1    
 for j=2:n-1   
    y(i,j)=l(i-1,j-1)+l(i-1,j)+l(i-1,j+1)+l(i,j-1)+l(i,j)+l(i,j+1)+l(i+1,j-1)+l(i+1,j)+l(i+1,j+1);
    y(i,j)=y(i,j)/9;  %LoG算子提取边缘后，对结果进行均值滤波以去除噪声，为下一步hough变换提取直线作准备
 end
end
axes(handles.axes4);imshow(y);%title('均值滤波器处理后')
%% 二值化
q=im2uint8(y);
q=(q>100)*255;
% [m,n]=size(q);
% for i=1:m    
%  for j=1:n  
%     if q(i,j)>100;   %设置二值化的阈值为80
%         q(i,j)=255; %对图像进行二值化处理，使图像边缘更加突出清晰
%     else
%         q(i,j)=0; 
%     end
%  end
% end
axes(handles.axes5);imshow(q);%title('二值化处理后');
%% 检测直线
%Hough变换检测直线，使用（a，p）参数空间，a∈[0,180],p∈[0,2d]
a=180; %角度的值为0到180度
d=round(sqrt(m^2+n^2)); %图像对角线长度为p的最大值
s=zeros(a,2*d); %存储每个(a,p)个数
[aa,bb]=find(q==255);
z=zeros(1,2);   %bu用元胞存储每个被检测的点的坐标

ll=length(aa);%ll个
for loop=1:ll
        i=aa(loop);
     j=bb(loop);
            for k=1:a
                p = round(i*cos(pi*k/180)+j*sin(pi*k/180));%对每个点从1到180度遍历一遍，取得经过该点的所有直线的p值(距离)（取整）
                if(p > 0)%若p大于0，则将点存储在（d，2d）空间
                    s(k,d+p)=s(k,d+p)+1;%（a，p）相应的累加器单元加一
                    if(s(k,d+p) == 100)
                       z=[z;k,p];
                    end
                else
                    ap=abs(p)+1;%若p小于0，则将点存储在（0，d）空间
                    s(k,ap)=s(k,ap)+1;
                   if(s(k,ap) ==100)
                    z=[z;k-180,ap]; 
                    end 
                end
            end
       
end
%% 显示效果
axes(handles.axes6);%figure(22);    
hold off;
plot(0,0);
axis([0 640 0 480]);
% title('直线检测');
hold on;

% xx=z(:,2)./sin((pi.*z(:,1)./180));
yy=z(:,2)./sin(pi.*z(:,1)./180);
k=-cot(pi.*z(:,1)./180);

a=find(z(:,1)~=0);
% D=[];
for i=2:(1+length(a))    
 yn=k(i).*x+yy(i);
 plot(yn,480-x,'-r','LineWidth',1); %%%%%%%%%%%%%%%55  z=-90
  D(i)=yn(1);
end
%%
axes(handles.axes7);% figure(201)
hold off;
for i=-800:800
    DD(i+801)=abs(sum(D(D>((i-1)*10+1)&D<i*10)));
end

bar(-800:800,DD)
axis([-200 200 0 max(DD)+10]);
% title('按视平线上的截距统计');
% xlabel('视平线截距')
% ylabel('直线数目统计')
%%
axes(handles.axes8);%figure(202)
hold off;
plot(z(:,1),z(:,2),'r+')
hold on
i=find(DD==max(DD(801:801+70)))-801;%
zz=[z((D>((i-1)*10+1)&D<i*10),1),z((D>((i-1)*10+1)&D<i*10),2)];
plot(zz(:,1),zz(:,2),'go')%圈出消失点在前方的复合条件的点
axis([-30 180 0 480]);
% title('符合条件直线簇特征值');
% xlabel('角度')
% ylabel('距离')
%% 识别结果
axes(handles.axes9);%figure(203)
a=find(zz(:,1)~=0);
xx=zz(:,2)./sin((pi.*zz(:,1)./180));
yy=zz(:,2)./sin(pi.*zz(:,1)./180);
k=-cot(pi.*zz(:,1)./180);
% title('筛选后的直线');
imshow(f);hold on;
for i=1:(length(a))    
 yn=k(i).*x+yy(i);
 plot(yn,x,'g')
end
 ii=find((abs(yy-mean(yy))-min(abs(yy-mean(yy))))<15);
for i=1:length(ii)%%max(1,ii-ceil(0.25*length(k))):min(length(k),ii+ceil(0.25*length(k)))  
 yn=k(ii(i)).*x+yy(ii(i));
 plot(yn,x,'r--','linewidth',5)
end

axis([0 640 0 480]);
clear D
end


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
hold off
plot(0,0)
axes(handles.axes2)
hold off
plot(0,0)
axes(handles.axes3)
hold off
plot(0,0)
axes(handles.axes4)
hold off
plot(0,0)
axes(handles.axes5)
hold off
plot(0,0)
axes(handles.axes6)
hold off
plot(0,0)
axes(handles.axes7)
hold off
plot(0,0)
axes(handles.axes8)
hold off
plot(0,0)
axes(handles.axes9)
hold off
plot(0,0)