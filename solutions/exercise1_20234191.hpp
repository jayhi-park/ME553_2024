//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <string.h>
#include <iostream>

using namespace std;

class Joint {
  private:
    string name_;
    Joint* parent_;

    // position
    Eigen::Vector3d r_12_; // position vector between parent joint(1) and current joint(2)

    // orientation
    Eigen::Matrix3d R_12_;   // Rotation matrix between parent joint(1) and current joint(2). 
    Eigen::Matrix3d R_22p_ = Eigen::Matrix3d::Identity(); // Rotation matrix by revolution joint. If it is fixed joint, it should be identity matrix.

    Eigen::Matrix3d CalRotationMatrix(double r, double p, double y){
      Eigen::Matrix3d R, Rx, Ry, Rz;

      Rx << 1.0d,   0.0d,    0.0d,
            0.0d, cos(r), -sin(r),
            0.0d, sin(r),  cos(r);

      Ry <<  cos(p), 0.0d, sin(p),
               0.0d, 1.0d,   0.0d,
            -sin(p), 0.0d, cos(p);

      Rz << cos(y), -sin(y), 0.0d,
            sin(y),  cos(y), 0.0d,
              0.0d,    0.0d, 1.0d;  

      R = Rx * Ry * Rz;
      return R;
    }

  public:
    Joint(string name, Joint* parent, const double rpy[], const double xyz[], double theta) : name_(name), parent_(parent){
      // make positon vector
      r_12_ << xyz[0], xyz[1], xyz[2];

      // cal rotation matrix
      R_22p_ = CalRotationMatrix(theta, 0, 0);
      R_12_ = CalRotationMatrix(rpy[0], rpy[1], rpy[2]);
    }

    string GetName(){
      return name_;
    }

    Joint* GetParent() {
      return parent_;
    }

    void ChangeFramefrom2toB(Eigen::Vector3d& position_vector_2) {
      // the frame of position_vector_2 and the 2 frame of joint_12 should be same.
      position_vector_2 = r_12_ + R_12_ * R_22p_ * position_vector_2;

      if (parent_ == nullptr){ // current joint is root joint.
        return ;
      }
      parent_->ChangeFramefrom2toB(position_vector_2);
    }
};


/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  // define ANYmal joints
  vector<Joint*> joints;

  joints.push_back(new Joint("base_LH_HAA", nullptr, (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){-0.2999, 0.104, 0.0}, 0));
  joints.push_back(new Joint("LH_HAA", joints.back(), (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, -gc[6+7])); // axis = [-1 0 0]
  joints.push_back(new Joint("LH_HIP_LH_hip_fixed", joints.back(), (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){0, 0, 0}, 0));
  joints.push_back(new Joint("LH_hip_fixed_LH_HFE", joints.back(), (const double[]){0, 0, 1.57079632679}, (const double[]){-0.0599, 0.08381, 0.0}, 0));
  joints.push_back(new Joint("LH_HFE", joints.back(), (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[6+8])); // axis = [1 0 0]
  joints.push_back(new Joint("LH_THIGH_LH_thigh_fixed", joints.back(), (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0));
  joints.push_back(new Joint("LH_thigh_fixed_LH_KFE", joints.back(), (const double[]){0, 0, 1.57079632679}, (const double[]){-0.0, 0.1003, -0.285}, 0));
  joints.push_back(new Joint("LH_KFE", joints.back(), (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[6+9])); // axis = [1 0 0]
  joints.push_back(new Joint("LH_shank_LH_shank_fixed", joints.back(), (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0));
  joints.push_back(new Joint("LH_Shank_fixed_LH_FOOT", joints.back(), (const double[]){0, 0, 0}, (const double[]){-0.08795, 0.01305, -0.33797}, 0));

  Eigen::Vector3d position_vector_be = Eigen::Vector3d::Zero();

  joints.back()->ChangeFramefrom2toB(position_vector_be);

  // r_wb
  Eigen::Vector3d position_vector_wb = gc.segment(0, 2);

  // R_wb
  Eigen::Quaternion<double> q(gc[3], gc[4], gc[5], gc[6]);
  Eigen::Matrix3d R_wb = q.normalized().toRotationMatrix();

  for (auto ptr : joints){
    delete ptr;
  }

  return position_vector_wb + R_wb * position_vector_be; // position_vector_be; /// replace this
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
