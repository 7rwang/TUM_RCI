% create a transformation matrix Tr with 
% the 3 angles of rotation as parameter and the translation vector 
Tr=set_rot(20,30,50,[2,3,4]);

% generate random 3D points in plane in 3m distance to camera
P=[randn(1,30)*2;randn(1,30)*2;ones(1,30)*3;ones(1,30)];

% calculate the position of the points after applying transformation R
P1=Tr*P;

% calculate image coordinates for a unifocal camera f/px=1 of the 
% generated points
u(1:2:60)=P(1,:)./P(3,:);
u(2:2:60)=P1(1,:)./P1(3,:);
v(1:2:60)=P(2,:)./P(3,:);
v(2:2:60)=P1(2,:)./P1(3,:);

%add noise to image detection.
%u=u+(randn(1,60)*0.002);
%v=v+(randn(1,60)*0.002);

[R,T,H]=homography(u,v);
R
Tr(1:3,1:3)

[ T Tr(1:3,4)/norm(Tr(1:3,4))]


%calculate angular error in rotation around the 3 axes
%rodrigues(Tr(1:3,1:3)*T(1:3,1:3)')*180/pi
%calculate angular errror of the estimated translation vector T
%acos(T(1:3,4)'*Tr(1:3,4)/(norm(T(1:3,4))*norm(Tr(1:3,4))))*180/pi
