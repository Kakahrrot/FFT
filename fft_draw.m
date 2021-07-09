N = 16;
n = 0:N-1;  
t = 2*pi*n/N;
xn = 3*cos(t) + sin(8*t);
yn = 3*sin(t) + cos(5*t);
plot_fft_with_circle_and_arrow(xn,yn,4);
function plot_fft_with_circle_and_arrow(x,y,m)
    N = length(x);
    z = x+1i*y;
    c = fft(z)/N;
    c(N/2+1) = c(N/2+1)/2;
    zi = interp_fft(z,m);
    xi = real(zi);
    yi = imag(zi);
    M = length(z)*m;
    n = 0:M;
    t = 2*pi*n/M;
    plot(x,y,'o',xi,yi,'-');
    for n=1:length(t)
        x0 = real(c(1));
        y0 = imag(c(1));
        plot(x,y,'o',xi,yi,'-');
        hold on;
        for k=1:N/2
            plot(x0+abs(c(k+1))*cos(t),y0+abs(c(k+1))*sin(t),'--');
            hold on;
            xt = x0; yt = y0;
            x0 = x0 + real(c(k+1)*exp(1i*k*t(n)));
            y0 = y0 + imag(c(k+1)*exp(1i*k*t(n)));
            draw_arrow([xt, yt], [x0, y0]);
            plot(x0+abs(c(N+1-k))*cos(t),y0+abs(c(N+1-k))*sin(t),'--');
            hold on;
            xt = x0; yt = y0;
            x0 = x0 + real(c(N+1-k)*exp(-1i*k*t(n)));
            y0 = y0 + imag(c(N+1-k)*exp(-1i*k*t(n))); 
            draw_arrow([xt, yt], [x0, y0]);
        end
        plot(x0,y0,'*')
        axis equal
        axis([-5,5,-5,5])
        hold off;
        pause(0.001)
    end
end
function draw_arrow(p1, p2)
    dp = p2-p1;
    quiver(p1(1),p1(2),dp(1),dp(2),0)
end
function cm = expand_fft_m(c,m)
    N = length(c);
    cm = [m*c(1:N/2) m*c(N/2+1)/2 zeros(1,(m-1)*N-1) m*c(N/2+1)/2 m*c(N/2+2:end)];
end
function yi = interp_fft(y,m)
    yi = ifft(expand_fft_m(fft(y),m));
end