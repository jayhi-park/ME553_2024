#pragma once

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
    Eigen::Vector3d r_12_1_; // position vector between parent joint(1) and current joint(2) in parent frame(1)
    Eigen::Vector3d r_12_b_; // position vector between parent joint(1) and current joint(2) in base frame
    Eigen::Vector3d r_b2_b_; // position vector between root joint(0) and current joint(2) in base frame

    // orientation
    Eigen::Matrix3d R_12_;   // Rotation matrix between parent joint(1) and current joint(2). 
    Eigen::Matrix3d R_22p_ = Eigen::Matrix3d::Identity(); // Rotation matrix by revolution joint. If it is fixed joint, it should be identity matrix.
    Eigen::Matrix3d R_b1_ = Eigen::Matrix3d::Zero();

    // rotation axis
    Eigen::Vector3d P_2_1_; // Rotation axis in parent frame(1)
    Eigen::Vector3d P_2_2_; // Rotation axis in current frame(2)

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

      R = Rz * Ry * Rx;
      return R;
    }

  public:
    Joint(string name, Joint* parent, const double rpy[], const double xyz[], double theta, const double p[]) : name_(name), parent_(parent){
      // make positon vector
      r_12_1_ << xyz[0], xyz[1], xyz[2];

      // cal rotation matrix
      R_22p_ = CalRotationMatrix(theta, 0, 0);
      R_12_ = CalRotationMatrix(rpy[0], rpy[1], rpy[2]);

      // set rotation axis
      P_2_2_ << p[0], p[1], p[2];
    }

    string GetName() {
      return name_;
    }

    Joint* GetParent() {
      return parent_;
    }

    Eigen::Matrix3d GetR_12() {
      return R_12_;
    }

    Eigen::Matrix3d GetR_22p() {
      return R_22p_;
    }

    Eigen::Matrix3d GetR_B1() {
      if (parent_ == nullptr){ // current joint is root joint.
        return Eigen::Matrix3d::Identity();
      } 

      if (R_b1_ == Eigen::Matrix3d::Zero()) {
        R_b1_ = parent_->GetR_B1() * parent_->GetR_12() * parent_->GetR_22p();
      }

      return R_b1_;
    }

    void Calr_12_B() {
      r_12_b_ = GetR_B1() * r_12_1_;
    }

    Eigen::Vector3d Getr_12_B() {
      return r_12_b_;
    }

    Eigen::Vector3d GetP_2_B() {
      // return GetR_B1() * P_2_1_;
      return GetR_B1() * R_12_ * R_22p_ * P_2_2_;
    }
};

class ArticulatedSystem {
  private:
    vector<Joint*> joints;

    // root velocity
    Eigen::Vector3d root_v_;
    
    // root orientation
    Eigen::Vector3d root_w_;

    // position vectors in root(base) frame
    vector<Eigen::Vector3d> r_b_vec_; // r_b2_b

    // rotation axis in root(base) frame
    vector<Eigen::Vector3d> P_b_vec_; // P_2_b

    // world frame
    Eigen::Vector3d r_wb_;
    Eigen::Matrix3d R_wb_;
    vector<Eigen::Vector3d> r_w_vec_; // r_w2_w
    vector<Eigen::Vector3d> P_w_vec_; // P_2_w

    
  public:
    ArticulatedSystem(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
      joints.push_back(new Joint("base_LH_HAA", nullptr, (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){-0.2999, 0.104, 0.0}, 0, (const double[]){0, 0, 0}));
      joints.push_back(new Joint("LH_HAA", joints.back(), (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, -gc[6+7], (const double[]){1, 0, 0})); // axis = [-1 0 0]
      joints.push_back(new Joint("LH_HIP_LH_hip_fixed", joints.back(), (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
      joints.push_back(new Joint("LH_hip_fixed_LH_HFE", joints.back(), (const double[]){0, 0, 1.57079632679}, (const double[]){-0.0599, 0.08381, 0.0}, 0, (const double[]){0, 0, 0}));
      joints.push_back(new Joint("LH_HFE", joints.back(), (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[6+8], (const double[]){1, 0, 0})); // axis = [1 0 0]
      joints.push_back(new Joint("LH_THIGH_LH_thigh_fixed", joints.back(), (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
      joints.push_back(new Joint("LH_thigh_fixed_LH_KFE", joints.back(), (const double[]){0, 0, 1.57079632679}, (const double[]){-0.0, 0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
      joints.push_back(new Joint("LH_KFE", joints.back(), (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[6+9], (const double[]){1, 0, 0})); // axis = [1 0 0]
      joints.push_back(new Joint("LH_shank_LH_shank_fixed", joints.back(), (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
      joints.push_back(new Joint("LH_Shank_fixed_LH_FOOT", joints.back(), (const double[]){0, 0, 0}, (const double[]){-0.08795, 0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));

      // calculate position vectors in root(base) frame
      for (auto& joint : joints) {
        joint->Calr_12_B();
        r_b_vec_.push_back(joint->Getr_12_B()); // r_12_b
      }

      for (int i = 1; i < r_b_vec_.size(); i++) {
        r_b_vec_[i] += r_b_vec_[i - 1]; // r_b2_b
      }  

      // get rotation axis vectors in root(base) frame
      for (auto& joint : joints) {
        P_b_vec_.push_back(joint->GetP_2_B()); // P_2_b
      }
    }

    vector<Eigen::Vector3d> Getr_b_vec() {
      return r_b_vec_;
    }  

    vector<Eigen::Vector3d> GetP_b_vec() {
      return P_b_vec_;
    }  

    void SetWorldFrame(Eigen::Vector3d r_wb, Eigen::Matrix3d R_wb) {
      r_wb_ = r_wb;
      R_wb_ = R_wb;

      CalPoseFramefromBtoW();
    }

    void CalPoseFramefromBtoW() {
      r_w_vec_.clear();
      P_w_vec_.clear();

      // cal world frame case
      for (int i = 0; i < r_b_vec_.size(); i++) {
        r_w_vec_.push_back(r_wb_ + R_wb_ * r_b_vec_[i]);
        P_w_vec_.push_back(R_wb_ * P_b_vec_[i]);
      }
    }

    Eigen::MatrixXd GetPositionalJ() {
      Eigen::MatrixXd J = Eigen::Matrix<double, 3, 9>::Zero();

      // I
      J.block<3,3>(0, 0) = Eigen::Matrix3d::Identity();

      // -r_we_w
      J(0, 4) =  (r_w_vec_.back() - r_wb_)(2);
      J(0, 5) = -(r_w_vec_.back() - r_wb_)(1);
      J(1, 3) = -(r_w_vec_.back() - r_wb_)(2);
      J(1, 5) =  (r_w_vec_.back() - r_wb_)(0);
      J(2, 3) =  (r_w_vec_.back() - r_wb_)(1);
      J(2, 4) = -(r_w_vec_.back() - r_wb_)(0);

      // revolute joint is 1, 4, 7
      J.block<3, 1>(0, 6) = P_w_vec_[1].cross(r_w_vec_.back() - r_w_vec_[1]);
      J.block<3, 1>(0, 7) = P_w_vec_[4].cross(r_w_vec_.back() - r_w_vec_[4]);
      J.block<3, 1>(0, 8) = P_w_vec_[7].cross(r_w_vec_.back() - r_w_vec_[7]); 

      return J;
    } 

    Eigen::MatrixXd GetAngularJ() {
      Eigen::MatrixXd J = Eigen::Matrix<double, 3, 9>::Zero();

      // I
      J.block<3,3>(0, 3) = Eigen::Matrix3d::Identity();

      // revolute joint is 1, 4, 7
      J.block<3, 1>(0, 6) = P_w_vec_[1];
      J.block<3, 1>(0, 7) = P_w_vec_[4];
      J.block<3, 1>(0, 8) = P_w_vec_[7];

      return J;
    } 
};

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  // r_wb
  Eigen::Vector3d r_wb = gc.segment(0, 3);

  // R_wb
  Eigen::Quaternion<double> q(gc[3], gc[4], gc[5], gc[6]);
  Eigen::Matrix3d R_wb = q.normalized().toRotationMatrix();

  // create system and set world frame
  ArticulatedSystem system(gc, gv);
  system.SetWorldFrame(r_wb, R_wb);

  // get generalized coor
  Eigen::VectorXd generalized_vel = Eigen::Vector<double, 9>::Zero();
  generalized_vel.block<3,1>(0, 0) = gv.segment(0, 3); // v
  generalized_vel.block<3,1>(3, 0) = gv.segment(3, 3); // w
  generalized_vel(6) = -gv[6+6];
  generalized_vel(7) = gv[6+7];
  generalized_vel(8) = gv[6+8];

  return system.GetPositionalJ() * generalized_vel;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  // r_wb
  Eigen::Vector3d r_wb = gc.segment(0, 3);

  // R_wb
  Eigen::Quaternion<double> q(gc[3], gc[4], gc[5], gc[6]);
  Eigen::Matrix3d R_wb = q.normalized().toRotationMatrix();

  // create system and set world frame
  ArticulatedSystem system(gc, gv);
  system.SetWorldFrame(r_wb, R_wb);

  // get generalized coor
  Eigen::VectorXd generalized_vel = Eigen::Vector<double, 9>::Zero();
  generalized_vel.block<3,1>(0, 0) = gv.segment(0, 3); // v
  generalized_vel.block<3,1>(3, 0) = gv.segment(3, 3); // w
  generalized_vel(6) = -gv[6+6];
  generalized_vel(7) = gv[6+7];
  generalized_vel(8) = gv[6+8];

  return system.GetAngularJ() * generalized_vel;
}