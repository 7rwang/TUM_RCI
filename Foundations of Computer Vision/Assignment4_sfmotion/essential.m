% [T,E]=essential(u,v)
%
% T - transformation matrix, E - essential matrix
% u,v - coordinate pairs (see below)
%
% calculates the essential matrix and
% the Transformation matrix
% input uni-focal point coordinates u11,u21,u12,u22,u13,u23,etc

function [T,E]=essential(u,v)

p=zeros(length(u)/2,9);
for i=1:length(u)/2
  p(i,1)=u(2*i-1)*u(2*i);
  p(i,2)=u(2*i-1)*v(2*i);
  p(i,3)=u(2*i-1);
  p(i,4)=v(2*i-1)*u(2*i);
  p(i,5)=v(2*i-1)*v(2*i);
  p(i,6)=v(2*i-1);
  p(i,7)=u(2*i);
  p(i,8)=v(2*i);
  p(i,9)=1;
end


[u s g]=svd(p);
E=reshape(g(:,9),3,3);
[u s g]=svd(E);
abs(1-s(1,1)/s(2,2))
s1=zeros(3,3);
s1(1,1)=(s(1,1)+s(2,2))/2;
s1(2,2)=s1(1,1);
%E=u*s1*g';
R=u*[0 -1 0;1 0 0; 0 0 1]*g';
if(det(R)<0)
  g(:,3)=-g(:,3); 
  R=u*[0 -1 0;1 0 0; 0 0 1]*g';
end
if(R(1,1)<0 | R(2,2)<0 | R(3,3)<0)
  R=u*[0 1 0;-1 0 0; 0 0 1]*g';
  if(det(R)<0)
    g(:,3)=-g(:,3); 
    R=u*[0 1 0;-1 0 0; 0 0 1]*g';
  end
end

T1=u*[0 -1 0;1 0 0; 0 0 1]*s1*u';
T2=u*[0 1 0;-1 0 0; 0 0 1]*s1*u';
if(T1(1,2)<0)
 T=[-T1(2,3),T1(1,3),-T1(1,2)]';
else
 T=[-T2(2,3),T2(1,3),-T2(1,2)]';
end

T=[R,T/norm(T)];

