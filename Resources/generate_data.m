% Generate field map

load('imagePhantom.mat')
maskobj = imagePhantom>0;
A = 20; % FM value in Hertz 
B = 3; % FM value in Hertz

x0=0;
y0=-22;
radius=8;
[n m]=size(imagePhantom);
x = (0:n-1)-n/2+1/2;
y = (0:m-1)-m/2+1/2;
circ = zeros(m,n);
[X Y] = meshgrid(x,y);

circ = (sqrt((X-x0).^2+(Y-y0).^2)<radius);

FM = A.*circ.*maskobj+maskobj.*(2*pi)*B;

% Generate k-space data
load kx_sp
load ky_sp

tsamp = 5e-6;  %
tt=col([0: tsamp : tsamp*(size(kx,1)-1)]);

A = fast_mr_v2(col(kx), col(ky), 0,  24, 64,64,1, 2*64,2*64,1, 5, tt, FM, 0, 5, 1,0,0,logical(ones(64,64)));

load sen 
As = sense(A,reshape(sen,64*64*1,4));

data = A* col(imagePhantom);
datas = As* col(imagePhantom);




