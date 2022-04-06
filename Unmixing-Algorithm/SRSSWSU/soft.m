function y = soft(x,T)

T = T + eps;%指的是1和离他最近的浮点数之间的距离。
y = max(abs(x) - T, 0);
y = y./(y+T) .* x;

