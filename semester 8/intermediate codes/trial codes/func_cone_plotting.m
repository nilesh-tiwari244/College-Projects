function []= func_cone_plotting(Aov,S1x,S2x,S3x,colour)
rad=tan((Aov/2));
p1=[0 0 0];
p2=1.5*[S1x S2x S3x];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'k')
hold on;
axis_length=1.5;
p2=[0 0 axis_length];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'b')
hold on;
p2=[0 axis_length 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'g')
hold on;
p2=[axis_length 0 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'r')
hold on;


p2=[0 0 -1*axis_length];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'b')
hold on;
p2=[0 -1*axis_length 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'g')
hold on;
p2=[-1*axis_length 0 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'r')
hold on;
pre=0.01;
for xx=-2:pre:2
    for yy=-2:pre:2
        for zz=-2:pre:2
            pp=S1x*(xx-S1x)+S2x*(yy-S2x)+S3x*(zz-S3x);
            if abs(pp)<0.01
                dis=sqrt((xx-S1x)^2+(yy-S2x)^2+(zz-S3x)^2);
                if abs(dis-rad)<0.01
                    p2=[xx yy zz];
                    v=[p2;p1];
                    plot3(v(:,1),v(:,2),v(:,3),colour)
                    hold on;
                end
            end
        end
    end
end
legend('Direction vector','Z axis','Y axis','X axis');
xlabel('X coordinate');
ylabel('Y coordinate');
zlabel('Z coordinate');
grid on;
hold on;
