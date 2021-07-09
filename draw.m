t = cputime;
I = imread('Doraemon.jfif');
% I = imread('fourier.jpeg');
% I = imgaussfilt(I,2);

% imshow(I);
%imshow(rgb2gray(I));
% imshow(edge(rgb2gray(I),'canny'));
% imshow(edge(rgb2gray(I),'sobel'));

canny_img = edge(rgb2gray(I),'canny');
sobel_img = edge(rgb2gray(I),'sobel');
img = sobel_img;
global count;
global pixel_distance;
global record;
global drop;
global interpolationM;
count = 1;
pixel_distance = 1;
record = zeros(size(img));

drop = 10;
interpolationM = 1;
figure;
hold on;
camroll(-90);
bfs(img);
s = sprintf('drop:%d\ninterpolationM:%d\nCPUtime:%.2fsecs', drop, interpolationM, cputime - t);
text(0, 3.5, s);
s = sprintf('./pic/Doraemon_%d_%d.png', drop, interpolationM);
% s = sprintf('./pic/fourier_%d_%d.png', drop, interpolationM);
saveas(gcf,s);

function plot_fft_point(x,y,m)
    N = length(x);
    z = x+1i*y;
    c = fft(z)/N;
%     disp(N);
    c(N/2+1) = c(N/2+1)/2;
    M = length(z)*m;
    n = 0:M;
    t = 2*pi*n/M;
    for n=1:length(t)
        x0 = real(c(1));
        y0 = imag(c(1));
        for k=1:N/2
            x0 = x0 + real(c(k+1)*exp(1i*k*t(n)))+ real(c(N+1-k)*exp(-1i*k*t(n)));
            y0 = y0 + imag(c(k+1)*exp(1i*k*t(n)))+imag(c(N+1-k)*exp(-1i*k*t(n)));
        end
        plot(x0,y0,'*')
        axis equal
%         axis([-5,5,-5,5])
        pause(0.01)
    end
end
function bfs(img)
    global count;
    global record;
    img_size= size(img);
    for i = 1:img_size(1)
        for j = 1:img_size(2)
          if img(i,j) == 1 && record(i,j) == 0
              mark(img, i, j);
              count = count + 1;
          end
        end
    end
end
function mark(img, x, y)
    global count;
    global record;
    global pixel_distance;
    global drop;
    global interpolationM;
    sz = size(img);
    x_stack = [x];
    y_stack = [y];
    xn = [x];
    yn = [y];
    while true
        if isempty(x_stack)
            break
        end
        x = x_stack(end);
        y = y_stack(end);
        x_stack(end) = [];
        y_stack(end) = [];

        for i = x - pixel_distance: x + pixel_distance
            for j = y - pixel_distance: y +pixel_distance
               if i >= 1 && j >= 1 && i <= sz(1) && j <= sz(2) && record(i, j) == 0 && img(i, j) == 1
                   record(i, j) = count;
                   x_stack(end + 1) = i;
                   y_stack(end + 1) = j;
                   xn(end+1) = i;
                   yn(end+1) = j;
               end
            end            
        end
    end
    
    k = 1: drop:length(xn);
        xn = xn(k)./sz(1) * 5;
        yn = yn(k)./sz(2) * 5;
        if mod(length(xn), 2) == 1
            xn(end) = [];
            yn(end) = [];
        end
        if length(xn) > 1
            plot_fft_point(xn, yn, interpolationM);
        end
end