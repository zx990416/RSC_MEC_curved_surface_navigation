function v = speed2D_cy(x,y,t)
v = zeros(length(t),1);
Lx = (max(x)-min(x))/2;
for i = 2:length(t)-1
    deltax = x(i+1)-x(i-1);
    if abs(deltax) > Lx
        deltax = 2*Lx - Lx;
    end
    v(i) = sqrt((deltax)^2+(y(i+1)-y(i-1))^2)/(t(i+1)-t(i-1));
end
v(1) = v(2);
v(end) = v(end-1);
end
