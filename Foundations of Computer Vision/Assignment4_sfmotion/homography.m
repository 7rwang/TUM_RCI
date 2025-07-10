function [R,T,H] = homography(u, v)

n=size(u,2)/2;
%construct H

h1=[ zeros(3,n);...
     -u(1:2:end); -v(1:2:end); -ones(1,n);...
     u(1:2:end).*v(2:2:end); v(1:2:end).*v(2:2:end); v(2:2:end);];
h2=[ u(1:2:end); v(1:2:end); ones(1,n);...
      zeros(3,n);...
     -u(1:2:end).*u(2:2:end); -v(1:2:end).*u(2:2:end); -u(2:2:end);];
h3=[ -u(1:2:end).*v(2:2:end); -v(1:2:end).*v(2:2:end); -v(2:2:end);...
      u(1:2:end).*u(2:2:end); v(1:2:end).*u(2:2:end); u(2:2:end);...
      zeros(3,n);];

H=[h1 h2 h3]';

[U S V] =svd(H);
H=reshape(V(:,9),3,3)';

% calculate solutions
[U S V] =svd(H'*H);
H=H/sqrt(S(2,2));
[U S V] =svd(H'*H);
u1=(sqrt(1-S(3,3))*V(:,1)+sqrt(S(1,1)-1)*V(:,3))/sqrt(S(1,1)-S(3,3));
u2=(sqrt(1-S(3,3))*V(:,1)-sqrt(S(1,1)-1)*V(:,3))/sqrt(S(1,1)-S(3,3));
U1=[V(:,2) u1 cross(V(:,2),u1)];
U2=[V(:,2) u2 cross(V(:,2),u2)];
W1=[H*V(:,2) H*u1 cross(H*V(:,2),H*u1)];
W2=[H*V(:,2) H*u2 cross(H*V(:,2),H*u2)];

R=W1*U1';
T=(H-R)*cross(u1,V(:,2));
T=T/norm(T);
%R=W2*U2'
