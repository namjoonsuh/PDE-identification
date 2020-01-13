x = linspace(-pi,pi,100);
fineX = linspace(-pi,pi,1000);
dx = diff(x);
dx = dx(1);
f = @(x) sin(11*(x-1)+0.8*pi);
sigma = 0.5;
rng(5)
u = f(x)+ normrnd(0,sigma,size(x));
w = 0.1;
K     = @(x) subplus(3/4*(1-(x/w).^2));
halfWindowNum = ceil(w/dx);
gridList = 1:halfWindowNum;
gridList = [-fliplr(gridList),0,gridList]*dx;
weight = K(gridList);
X = ones(size(gridList,2),1);

Yextend  = padarray(u,[0,halfWindowNum],'symmetric');
Y = circshift(flipud(gallery('circul',Yextend)),1);
Y = Y(1:2*halfWindowNum+1,1:length(u));
beta = X'*(weight'.*Y);
smoothedData = (X'*(weight'.*X))\beta;

figure
fY = abs(fft(f(fineX),251));
Pyy = fY.*conj(fY)/251;
F = 1000/251*(0:127);
plot(F,Pyy(1:128))
hold on 
fY = angle(fft(smoothedData,251));
Pyy = fY.*conj(fY)/251;
F = 1000/251*(0:127);
plot(F,Pyy(1:128))


