#include <iostream>
#include <math.h>
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/PositionSensor.hpp>
#include <webots/Supervisor.hpp>
#include <eigen/Eigen/Dense>

using namespace webots;
using namespace Eigen;

const double pi = 3.14159;

//Creates rotational matrix e^(w_hat*theta)
Matrix3d makeRotMat(Vector3d w, double theta)
{
  //Skew matrix
  Matrix3d w_hat;
  w_hat<<0, -w(2), w(1),
          w(2), 0, -w(0),
          -w(1), w(0), 0;
  
  Matrix3d I;
  I = Matrix3d::Identity();
  MatrixXd result= I + w_hat*sin(theta)+w_hat*w_hat*(1-cos(theta));
  return result;
}

//Creates translation vector t
Vector3d makeTransVec(Matrix3d e,Vector3d w, Vector3d v,double theta)
{
  Matrix3d I = Matrix3d::Identity();
  MatrixXd t = (I-e)*(w.cross(v))+w*w.transpose()*v*theta;
  return t;
}

//Creates 4x4 matrix
MatrixXd matExp(Matrix3d e, Vector3d t)
{
  Matrix4d res = Matrix4d::Zero();
  res<<e;
  res(0,3) = t(0);
  res(1,3) = t(1);
  res(2,3) = t(2);
  res(3,3) = 1;
  return res;
}

Matrix4d ForwardKinematics(Vector3d theta, Matrix4d M)
{

  //Origin point for joint 1
  Vector3d q1;
  q1<<0, 0.67183, 0;
  //Rotation along axis
  Vector3d w1;
  w1<<0,1,0;
  
  //Linear motion along axis
  VectorXd v1= -w1.cross(q1);
  
  //Rotational matrix (per Roudriges' rotational formula)
  Matrix3d e1 = makeRotMat(w1,theta(0));
  //Translation vector
  Vector3d t1 = makeTransVec(e1,w1,v1, theta(0));
  //Compose exponential matrix
  Matrix4d exp1 = matExp(e1,t1);
  
  
  //Same goes for joints 2 and 3
  Vector3d q2;
  q2<<0, 0.67183, -0.23368;
  Vector3d w2;
  w2<<0,0,-1;
  
  VectorXd v2 = -w2.cross(q2);
  Matrix3d e2 = makeRotMat(w2,theta(1));
  Vector3d t2 = makeTransVec(e2,w2,v2,theta(1));
  Matrix4d exp2 = matExp(e2,t2);
  
  Vector3d q3;
  q3<<0.4, 0.67183, -0.148;
  Vector3d w3;
  w3<<0,0,-1;
  
  VectorXd v3 = -w3.cross(q3);
  Matrix3d e3 = makeRotMat(w3,theta(2));
  Vector3d t3 = makeTransVec(e3,w3,v3,theta(2));
  Matrix4d exp3 = matExp(e3,t3);

  //Exponential formula
  Matrix4d final = exp1*exp2*exp3*M;
  return final;
}

Vector3d InverseKinematics(double px, double py, double pz)
{
  py-=0.67183;
  //Shoulder joint offset
  const double d1 = 0.152;
  const double r  = sqrt(px*px+pz*pz);
  double fi = atan2(px,pz) - pi/2;
  double alpha = atan2(d1, sqrt(r*r-d1*d1));
  double theta1 = fi - alpha;
  
  //Lengths of links
  const double a2 = 0.4318;
  const double a3 = 0.48895;
  double D = -(r*r-d1*d1+py*py-a2*a2-a3*a3)/(2*a2*a3);
  
  double theta3 = atan2(D,-sqrt(1-D*D)) - pi;
  
  const double s3 = sin(theta3);
  const double c3 = cos(theta3);
  double theta2 = atan2(a3*c3,a2+a3*s3)-atan2(py,sqrt(r*r-d1*d1));
  
  Vector3d result;
  result<<theta1, theta2, theta3;
  return result;
}

int main(int argc, char **argv) {
  Robot *robot = new Robot();
  
  //Desired configuration
  const Vector3d theta(0.6,-0.7,2);
  
  //Zero configuration matrix
  Matrix4d M;
  M<< 1, 0, 0, 0.4318,
      0, -1, 0, 1.16078,
      0, 0, -1, -0.152,
      0, 0, 0, 1;
  
  std::cout<<"M\n"<<M<<std::endl<<std::endl;
  
  Matrix4d T = ForwardKinematics(theta, M);
  std::cout<<"T\n"<<T<<std::endl;
  std::cout<<"\nINV\n"<<InverseKinematics(0.562537, 1.02967, -0.56902)<<std::endl;
  std::cout<<std::endl
  <<"Fwd(Inv)\n"<<ForwardKinematics(InverseKinematics(0.562537,1.02967,-0.56902), M)
  <<std::endl;
  //Initialise motors
  Motor* motor1 = robot->getMotor("joint1");
  Motor* motor2 = robot->getMotor("joint2"); 
  Motor *motor3 = robot->getMotor("joint3");
  
  //Set to desired configuration
  motor1->setPosition(theta(0));
  motor2->setPosition(theta(1));
  motor3->setPosition(theta(2));
  
  delete robot;
  return 0;
}
