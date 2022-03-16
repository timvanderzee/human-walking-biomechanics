function yp=grad5(y,dx)
% function yp=grad5(y,dx)
% 031091 KvS
% purpose: calculation of 5 point derivative
% inputs : y : vector containing signal as a function of x
%          dx: change in x between steps
% output : yp: derivative of y with respect to x
 
if nargin~=2
  disp('OOPS: grad5 expects 2 input arguments')
  disp('      hit any key to continue');pause
  return
end

%y=y(:);
[nrow,ncol]=size(y);


yp(1,1:ncol)=zeros(1,ncol);
yp(2,1:ncol )=zeros(1,ncol);
yp(nrow-1,1:ncol)=zeros(1,ncol);
yp(nrow,1:ncol)=zeros(1,ncol);

yp(1,1:ncol)=(y(2,:)-y(1,:))/dx;
yp(2,1:ncol)=(y(3,:)-y(1,:))/(2*dx);
yp(nrow-1,1:ncol)=(y(nrow,:)-y(nrow-2,:))/(2*dx);
yp(nrow,1:ncol)=(y(nrow,:)-y(nrow-1,:))/dx;

coef=[1 -8 0 8 -1];
for i=3:nrow-2;
  yp(i,:)=(coef*y(i-2:i+2,:))/(12*dx);
end;

%if flip; y=y';yp=yp';end;
