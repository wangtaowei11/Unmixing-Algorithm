function y = soft(x,T)

T = T + eps;%ָ����1����������ĸ�����֮��ľ��롣
y = max(abs(x) - T, 0);
y = y./(y+T) .* x;

