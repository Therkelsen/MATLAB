clear;

robot = rigidBodyTree('DataFormat', 'column');

body = rigidBody('body1');
jnt = rigidBodyJoint('q1', 'revolute');
jnt.JointAxis = [0 0 1];
jnt.HomePosition = 0;
tform = trvec2tform([0, 0, 1.320]);
setFixedTransform(jnt, tform);
body.Joint = jnt;
addBody(robot, body, 'base');

body = rigidBody('body2');
jnt = rigidBodyJoint('q2', 'revolute');
jnt.JointAxis = [0 1 0];
jnt.HomePosition = 0;
tform = trvec2tform([0.16, 0, 0]);
setFixedTransform(jnt, tform);
body.Joint = jnt;
addBody(robot, body, 'body1');

body = rigidBody('body3');
jnt = rigidBodyJoint('q3', 'revolute');
jnt.JointAxis = [0 1 0];
jnt.HomePosition = 0;
tform = trvec2tform([0, 0, 1.120]);
setFixedTransform(jnt, tform);
body.Joint = jnt;
addBody(robot, body, 'body2');

body = rigidBody('body4');
jnt = rigidBodyJoint('q4', 'revolute');
jnt.JointAxis = [1 0 0];
jnt.HomePosition = 0;
tform = trvec2tform([1.6492, 0, 0.40]);
setFixedTransform(jnt, tform);
body.Joint = jnt;
addBody(robot, body, 'body3');

body = rigidBody('body5');
jnt = rigidBodyJoint('q5', 'revolute');
jnt.JointAxis = [0 1 0];
jnt.HomePosition = 0;
tform = trvec2tform([0, 0, 0]);
setFixedTransform(jnt, tform);
body.Joint = jnt;
addBody(robot, body, 'body4');

body = rigidBody('body6');
jnt = rigidBodyJoint('q6', 'revolute');
jnt.JointAxis = [1 0 0];
jnt.HomePosition = 0;
tform = trvec2tform([0, 0, 0]);
setFixedTransform(jnt, tform);
body.Joint = jnt;
addBody(robot, body, 'body5');

body = rigidBody('endeffector');
tform = trvec2tform([0.44, 0, 0]);
setFixedTransform(body.Joint, tform);
addBody(robot, body, 'body6');
  
config = homeConfiguration(robot);
tform = getTransform(robot, config, 'endeffector', 'base')
 
%showdetails(robot);
%show(robot);
%vizztree = interactiveRigidBodyTree(robot);
 
t = (0:.5:20)';
count = length(t);
center = [0.5 0.5 0];
radius = 1.75;
th = t*(6*pi/t(end));
points = center + radius*[cos(th) sin(th) ones(size(th))];
for n = 1:count
    points(n,3) = 0.5+n*0.075;
end

q0 = homeConfiguration(robot);
ndof = length(q0);
qs = zeros(count, ndof);
  
ik = inverseKinematics('RigidBodyTree', robot);
weights = [0, 0, 0, 1, 1, 1];
endEffector = 'endeffector';

qInit = q0;
for i = 1:count
    % Solve for each circle point
    point = points(i,:);
    qSol = ik(endEffector, trvec2tform(point), weights, qInit);
    % Store the configuration
    qs(i,:) = qSol;
    % Start from prior solution
    qInit = qSol;
end

% figure
% % show(robot,qs(1,:)');
% % for n = 1:count
% %     show(robot,qs(n,:)');
% %     hold on
% % end
% view(2)
% hold on
% plot3(points(:,1),points(:,2),points(:,3))
% axis([-2.5 2.5 -2.5 2.5 0 4])
 
framesPerSecond = 5;
r = rateControl(framesPerSecond);
for n = 1:10
   for i = 1:count
       show(robot,qs(i,:)','PreservePlot',false);
       drawnow
       waitfor(r);
   end
end

