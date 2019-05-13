function varargout = main_all(varargin)
% MAIN_ALL MATLAB code for main_all.fig
%      MAIN_ALL, by itself, creates a new MAIN_ALL or raises the existing
%      singleton*.
%
%      H = MAIN_ALL returns the handle to a new MAIN_ALL or the handle to
%      the existing singleton*.
%
%      MAIN_ALL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_ALL.M with the given input arguments.
%
%      MAIN_ALL('Property','Value',...) creates a new MAIN_ALL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_all_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_all_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_all

% Last Modified by GUIDE v2.5 19-Feb-2018 23:38:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_all_OpeningFcn, ...
                   'gui_OutputFcn',  @main_all_OutputFcn, ...
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


% --- Executes just before main_all is made visible.
function main_all_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_all (see VARARGIN)

% Choose default command line output for main_all
handles.output = hObject;



tic
[filename pathname] = uigetfile({'*.jpg';'*.bmp'},'File Selector');
image = strcat(pathname, filename);
i=imread(image);
osize=dir('test2.jpg');
osize=osize.bytes;
axes(handles.axes1)
imshow(i);
i=imresize(i,[396 692]);
[m,n]=size(i);

j=zeros(round(m/3),round(n/3));
x=1;
y=1;

for ax=2:3:m
    for ax1=2:3:n
        j(x,y)=i(ax,ax1);
        y=y+1;
    end
    x=x+1;
    y=1;
end

j=uint8(j);
axes(handles.axes2)
imshow(j);
j=imresize (j,[132,231]);

imwrite(j, 'compressed.jpg');

csize=dir('compressed.jpg');
csize=csize.bytes;




e=imread('compressed.jpg');
[rw,cl]=size(e)
e=double(e);


rw1=1;
while rw1<rw 

    o=2;
    clm=1;
    cl1=1;

    while cl1<cl
        a=e(rw1,cl1); %Read first pixel from compressed image
        cl1=cl1+1;
        b=e(rw1,cl1); % Read second pixel from compressed image
        cl1=cl1+1;
        c=e(rw1,cl1); % Read second pixel from compressed image
        cl1=cl1+1;

        % now calculate ip1,ip2,ip3,ip4,ip5,ip6 using previopus values

        ip1=round((((a-(2*b)+c)*(o*o))+(((16*b)-(11*a)-(5*c))*o)+((28*a)-(14*b)+(4*c)))/18);
        o=o+1;
        ip2=round((((a-(2*b)+c)*(o*o))+(((16*b)-(11*a)-(5*c))*o)+((28*a)-(14*b)+(4*c)))/18);
        o=o+2;
        ip3=round((((a-(2*b)+c)*(o*o))+(((16*b)-(11*a)-(5*c))*o)+((28*a)-(14*b)+(4*c)))/18);
        o=o+1;
        ip4=round((((a-(2*b)+c)*(o*o))+(((16*b)-(11*a)-(5*c))*o)+((28*a)-(14*b)+(4*c)))/18);
        o=o+2;
        ip5=round((((a-(2*b)+c)*(o*o))+(((16*b)-(11*a)-(5*c))*o)+((28*a)-(14*b)+(4*c)))/18);
        o=o+1;
        ip6=round((((a-(2*b)+c)*(o*o))+(((16*b)-(11*a)-(5*c))*o)+((28*a)-(14*b)+(4*c)))/18);
        o=o+2;
        %now we have all nine values to be saved in decompressed image
        de(rw1,clm)=a;
        clm=clm+1;
        de(rw1,clm)=ip1;
        clm=clm+1;
        de(rw1,clm)=ip2;
        clm=clm+1;
        de(rw1,clm)=b;
        clm=clm+1;
        de(rw1,clm)=ip3;
        clm=clm+1;
        de(rw1,clm)=ip4;
        clm=clm+1;
        de(rw1,clm)=c;
        clm=clm+1;
        de(rw1,clm)=ip5;
        clm=clm+1;
        de(rw1,clm)=ip6;
        clm=clm+1;
        o=2;

    end
    rw1=rw1+1;
end



%de=uint8(de);

de=resizem(de,[132 693]);

[rw22,cl22]=size(de);
cl2=1;
while cl2<cl22
    oo=2;
    rw2=1;
    row=1;
    while rw2<rw22
        oo=2;
        a2=de(rw2,cl2); %Read first pixel from compressed image
        rw2=rw2+1;
        b2=de(rw2,cl2); % Read second pixel from compressed image
        rw2=rw2+1;
        c2=de(rw2,cl2); % Read second pixel from compressed image
        rw2=rw2+1;

        % now calculate ip1,ip2,ip3,ip4,ip5,ip6 using previopus values

        ip12=round((((a2-(2*b2)+c2)*(oo*oo))+(((16*b2)-(11*a2)-(5*c2))*oo)+((28*a2)-(14*b2)+(4*c2)))/18);
        oo=oo+1;
        ip22=round((((a2-(2*b2)+c2)*(oo*oo))+(((16*b2)-(11*a2)-(5*c2))*oo)+((28*a2)-(14*b2)+(4*c2)))/18);
        oo=oo+2;
        ip32=round((((a2-(2*b2)+c2)*(oo*oo))+(((16*b2)-(11*a2)-(5*c2))*oo)+((28*a2)-(14*b2)+(4*c2)))/18);
        oo=oo+1;
        ip42=round((((a2-(2*b2)+c2)*(oo*oo))+(((16*b2)-(11*a2)-(5*c2))*oo)+((28*a2)-(14*b2)+(4*c2)))/18);
        oo=oo+2;
        ip52=round((((a2-(2*b2)+c2)*(oo*oo))+(((16*b2)-(11*a2)-(5*c2))*oo)+((28*a2)-(14*b2)+(4*c2)))/18);
        oo=oo+1;
        ip62=round((((a2-(2*b2)+c2)*(oo*oo))+(((16*b2)-(11*a2)-(5*c2))*oo)+((28*a2)-(14*b2)+(4*c2)))/18);
        oo=oo+2;
        %now we have all nine values to be saved in decompressed image
        de2(row,cl2)=a2;
        row=row+1;
        de2(row,cl2)=ip12;
        row=row+1;
        de2(row,cl2)=ip22;
        row=row+1;
        de2(row,cl2)=b2;
        row=row+1;
        de2(row,cl2)=ip32;
        row=row+1;
        de2(row,cl2)=ip42;
        row=row+1;
        de2(row,cl2)=c2;
        row=row+1;
        de2(row,cl2)=ip52;
        row=row+1;
        de2(row,cl2)=ip62;
        row=row+1;
        
    end
    cl2=cl2+1;
end

de2=uint8(de2);
imwrite(de2,'decompressed.jpg');


%de2=imread('decompressed.jpg');
%axes(hand

figure()
imshow(de2);

de2=double(de2(:));
de2a=max(de2(:));
de2i=min(de2(:));
mse=std(de2(:));
snr=20*log10((de2a-de2i)./mse)

toc

compressedratio=osize/csize;
compressedratio=round(compressedratio);
h=msgbox(sprintf('compression ratio is %d',compressedratio));

i=double(i);
dd=imread('decompressed.jpg');
dd=imresize(dd,[396 692]);
dd=double(dd);
SE=(i-dd).^2;
mse=sum(sum(SE))/(m*n);
PSNR=10*log10(255^2/mse);
%h=msgbox(sprintf('PSNR is %d',PSNR));
sprintf('PSNR is %d',PSNR)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_all wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_all_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
