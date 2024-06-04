#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <set>
#include <vector>
#include <string.h>
#include <iostream>

using namespace std;

class Link
{
  public:
    int idx_;
    string name_;
    int parent_joint_idx_ = -1;

    // local
    Eigen::Vector3d r_po_p_; // position vector between parent joint(p) and origin(o) in parent frame(p)
    double mass_;
    Eigen::Matrix3d inertia_o_; // in origin frame(o)

    // world
    Eigen::Vector3d r_po_w_; // in world frame(w)
    Eigen::Vector3d r_wo_w_; // in world frame(w)
    Eigen::Matrix3d inertia_w_; // in world frame(w)

    // rotation matrix
    Eigen::Matrix3d R_po_; // rotation matrix from origin frame(o) to parent frame(p)
    Eigen::Matrix3d R_wp_; // rotation matrix from parent frame(p) to world frame(w)

    // Jacobian
    Eigen::MatrixXd Jacobian_p_; // positional Jacobian of origin
    Eigen::MatrixXd Jacobian_a_; // angular Jacobian of origin

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

    Link(int idx, string name, int parent_joint_idx, const double rpy[], const double xyz[], double mass, const double inertia[]) : 
      idx_(idx), name_(name), parent_joint_idx_(parent_joint_idx), mass_(mass){
        
      // origin in parent frame
      r_po_p_ << xyz[0], xyz[1], xyz[2];

      // inertia in origin frame
      inertia_o_ << inertia[0], inertia[1], inertia[2],
                    inertia[1], inertia[3], inertia[4], 
                    inertia[2], inertia[4], inertia[5];

      // cal rotation matrix
      R_po_ = CalRotationMatrix(rpy[0], rpy[1], rpy[2]);
    }

    Link(int idx, string name, int parent_joint_idx, Eigen::Matrix3d R_po, Eigen::Vector3d r_po_p, double mass, Eigen::Matrix3d inertia_o) : 
      idx_(idx), name_(name), parent_joint_idx_(parent_joint_idx), R_po_(R_po), r_po_p_(r_po_p), mass_(mass), inertia_o_(inertia_o) {}
};

class Joint
{
  public:
    int idx_;
    int moving_idx_ = -1;
    string name_;
    int parent_link_idx_ = -1;
    int child_link_idx_ = -1;

    // local
    Eigen::Vector3d r_pc_p_; // position vector between parent joint(p) and current joint(c) in parent frame(p)

    // world
    Eigen::Vector3d r_pc_w_; // position vector between parent joint(p) and current joint(c) in world frame(w)
    Eigen::Vector3d r_wc_w_; // position vector between world(w) and current joint(c) in world frame(w)

    Eigen::Vector3d v_c_w_; // velocity of current frame in world frame(w)
    Eigen::Vector3d v_cp_w_; // velocity of current prime frame in world frame(w)
    Eigen::Vector3d omega_c_w_; // angular velocity of current frame in world frame(w)
    Eigen::Vector3d omega_cp_w_; // angular velocity of current prime frame in world frame(w)

    Eigen::Vector3d a_c_w_; // acceleration of current frame in world frame(w)
    Eigen::Vector3d a_cp_w_; // acceleration of current prime frame in world frame(w)

    Eigen::Vector3d alpha_c_w_; // angular acceleration of current frame in world frame(w)
    Eigen::Vector3d alpha_cp_w_; // acceleration of current prime frame in world frame(w)

    // rotation matrix
    Eigen::Matrix3d R_pc_;                                // rotation matrix between parent joint(p) and current joint(c). 
    Eigen::Matrix3d R_ccp_ = Eigen::Matrix3d::Identity(); // rotation matrix by revolution joint. If it is fixed joint, it should be identity matrix.
    Eigen::Matrix3d R_wp_;                                // rotation matrix from parent frame(p) to world frame(w)

    // rotation axis
    Eigen::Vector3d P_c_; // rotation axis of joint in current frame(c)
    Eigen::Vector3d P_w_; // rotation axis of joint in world frame(w)

    // motion subspace matrix
    Eigen::VectorXd S_w_ = Eigen::VectorXd::Zero(6); // motion subspace matrix of joint in world frame(w)

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
    Joint(int idx, int moving_idx, string name, int parent_link_idx, int child_link_idx, const double rpy[], const double xyz[], double theta, const double p[]) : 
      idx_(idx), moving_idx_(moving_idx), name_(name), parent_link_idx_(parent_link_idx), child_link_idx_(child_link_idx){

      // current pos in parent frame
      r_pc_p_ << xyz[0], xyz[1], xyz[2];

      // set rotation axis
      P_c_ << p[0], p[1], p[2];
      Eigen::Vector3d rotation_vector = P_c_ * theta;

      // cal rotation matrix
      R_pc_ = CalRotationMatrix(rpy[0], rpy[1], rpy[2]);
      R_ccp_ = CalRotationMatrix(rotation_vector[0], rotation_vector[1], rotation_vector[2]);
    }
};

class ArticulatedSystem
{
private:
  // root
  Link* root; // body

  // world
  Eigen::Vector3d r_wb_w_; // position vector between world(w) and body(b) in world frame(w)
  Eigen::Vector3d v_b_w_; // velocity of body frame in world frame(w)
  Eigen::Vector3d omega_b_w_; // angular velocity of body frame in world frame(w)
  Eigen::Vector3d a_b_w_; // acceleration of body frame in world frame(w)
  Eigen::Vector3d alpha_b_w_; // angular acceleration of body frame in world frame(w)

  // rotation matrix
  Eigen::Matrix3d R_wb_; // rotation matrix from body frame(b) to world frame(w)

  // joint and link
  vector<Joint*> joints_;
  vector<Link*> links_;

  // generalized coordinate and velocity 
  Eigen::VectorXd gc_;
  Eigen::VectorXd gv_;

  // moving joints list
  vector<Joint*> moving_joints_;

  // num_of_moving_joints
  int num_of_moving_joints_;

  // mass matrix
  Eigen::MatrixXd M_;

  // non-linearities
  Eigen::VectorXd non_linearities_;
  Eigen::Vector3d gravity_ = Eigen::Vector3d(0, 0, -9.81);

public:
  ArticulatedSystem(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv)
  {
    // set root
    root = new Link(0, "base", -1, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.0, (const double[]){0, 0, 0, 0, 0, 0});
    links_.push_back(root);
    v_b_w_ = gv.segment(0, 3);
    omega_b_w_ = gv.segment(3, 3);
    a_b_w_ = Eigen::Vector3d::Zero();
    alpha_b_w_ = Eigen::Vector3d::Zero();

    gc_ = gc.segment(6, gc.size() - 6);
    gv_ = gv.segment(6, gv.size() - 6);

    // r_wb
    r_wb_w_ = gc.segment(0, 3);

    // R_wb
    Eigen::Quaternion<double> q(gc[3], gc[4], gc[5], gc[6]);
    R_wb_ = q.normalized().toRotationMatrix();

    // set links
    links_.push_back(new Link(links_.size(), "base_inertia", 0, (const double[]){0, 0, 0}, (const double[]){-0.018, -0.002, 0.024}, 6.222, (const double[]){0.017938806, 0.00387963, 0.001500772, 0.370887745, 6.8963e-05, 0.372497653}));
    links_.push_back(new Link(links_.size(), "face_front", 1, (const double[]){0, 0, 0}, (const double[]){0.042, -0.001, 0.004}, 0.73, (const double[]){0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938}));
    links_.push_back(new Link(links_.size(), "face_rear", 2, (const double[]){0, 0, 0}, (const double[]){0.042, -0.001, 0.004}, 0.73, (const double[]){0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938}));
    links_.push_back(new Link(links_.size(), "battery", 3, (const double[]){0, 0, 0}, (const double[]){-0.00067, -0.00023, -0.03362}, 5.53425, (const double[]){0.00749474794, 0.00016686282, 7.82763e-05, 0.0722338913, 1.42902e-06, 0.07482717535}));
    links_.push_back(new Link(links_.size(), "docking_hatch_cover", 4, (const double[]){0, 0, 0}, (const double[]){-0.003, 0.0, 0.005}, 0.065, (const double[]){0.00063283, 0.0, 3.45e-07, 0.00110971, 0.0, 0.00171883}));
    links_.push_back(new Link(links_.size(), "lidar_cage", 5, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0, 0, 0, 0}));
    links_.push_back(new Link(links_.size(), "lidar", 6, (const double[]){0, 0, 0}, (const double[]){-0.012, 0.001, -0.008}, 0.695, (const double[]){0.000846765, 6.9565e-05, 0.00027111, 0.001367583, 5.8984e-05, 0.001363673}));
    links_.push_back(new Link(links_.size(), "hatch", 7, (const double[]){0, 0, 0}, (const double[]){0.116, 0.0, 0.0758}, 0.142, (const double[]){0.001, 0.001, 0.001, 0.001, 0.001, 0.001}));
    
    links_.push_back(new Link(links_.size(), "LF_HAA", 8, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LF_HIP", 9, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LF_hip_fixed", 10, (const double[]){0, 0, 0}, (const double[]){0.048, 0.008, -0.003}, 0.74, (const double[]){0.001393106, 8.4012e-05, 2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "LF_HFE", 11, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LF_THIGH", 12, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LF_thigh_fixed", 13, (const double[]){0, 0, 0}, (const double[]){0.0, 0.018, -0.169}, 1.03, (const double[]){0.018644469, 5.2e-08, 1.0157e-05, 0.019312599, 0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "LF_KFE", 14, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LF_SHANK", 15, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LF_shank_fixed", 16, (const double[]){0, 0, 0}, (const double[]){0.03463, 0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, 2.142561e-05, 1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "LF_FOOT", 17, (const double[]){0, 0, 0}, (const double[]){0.00948, -0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, 2.63048e-06, 6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05}));

    links_.push_back(new Link(links_.size(), "RF_HAA", 18, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RF_HIP", 19, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RF_hip_fixed", 20, (const double[]){0, 0, 0}, (const double[]){0.048, -0.008, -0.0030}, 0.74, (const double[]){0.001393106, -8.4012e-05, 2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "RF_HFE", 21, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RF_THIGH", 22, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RF_thigh_fixed", 23, (const double[]){0, 0, 0}, (const double[]){0.0, -0.018, -0.169}, 1.03, (const double[]){0.018644469, -5.2e-08, 1.0157e-05, 0.019312599, -0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "RF_KFE", 24, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RF_SHANK", 25, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RF_shank_fixed", 26, (const double[]){0, 0, 0}, (const double[]){0.03463, -0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, -2.142561e-05, 1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "RF_FOOT", 27, (const double[]){0, 0, 0}, (const double[]){0.00948, 0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, -2.63048e-06, 6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05}));

    links_.push_back(new Link(links_.size(), "LH_HAA", 28, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LH_HIP", 29, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LH_hip_fixed", 30, (const double[]){0, 0, 0}, (const double[]){-0.048, 0.008, -0.003}, 0.74, (const double[]){0.001393106, -8.4012e-05, -2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "LH_HFE", 31, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LH_THIGH", 32, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LH_thigh_fixed", 33, (const double[]){0, 0, 0}, (const double[]){0.0, 0.018, -0.169}, 1.03, (const double[]){0.018644469, -5.2e-08, -1.0157e-05, 0.019312599, 0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "LH_KFE", 34, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LH_SHANK", 35, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LH_shank_fixed", 36, (const double[]){0, 0, 0}, (const double[]){-0.03463, 0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, -2.142561e-05, -1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "LH_FOOT", 37, (const double[]){0, 0, 0}, (const double[]){-0.00948, -0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, -2.63048e-06, -6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05}));

    links_.push_back(new Link(links_.size(), "RH_HAA", 38, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RH_HIP", 39, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RH_hip_fixed", 40, (const double[]){0, 0, 0}, (const double[]){-0.048, -0.008, -0.003}, 0.74, (const double[]){0.001393106, 8.4012e-05, -2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "RH_HFE", 41, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RH_THIGH", 42, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RH_thigh_fixed", 43, (const double[]){0, 0, 0}, (const double[]){0.0, -0.018, -0.169}, 1.03, (const double[]){0.018644469, 5.2e-08, -1.0157e-05, 0.019312599, -0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "RH_KFE", 44, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RH_SHANK", 45, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RH_shank_fixed", 46, (const double[]){0, 0, 0}, (const double[]){-0.03463, -0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, 2.142561e-05, -1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "RH_FOOT", 47, (const double[]){0, 0, 0}, (const double[]){-0.00948, 0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, 2.63048e-06, -6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05}));

    // set joints
    joints_.push_back(new Joint(joints_.size(), -1, "base_to_base_inertia", 0, 1, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_face_front", 0, 2, (const double[]){0, 0, 0}, (const double[]){0.4145, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_face_rear", 0, 3, (const double[]){0, 0, 3.14159265359}, (const double[]){-0.4145, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_battery", 0, 4, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_to_docking_hatch_cover", 0, 5, (const double[]){0, 0, 0}, (const double[]){0.343, 0.0, -0.07}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_to_lidar_cage", 0, 6, (const double[]){0, 0, 0}, (const double[]){-0.364, 0.0, 0.0735}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "lidar_cage_to_lidar", 6, 7, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0.0687}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_hatch", 0, 8, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));

    joints_.push_back(new Joint(joints_.size(), -1, "base_LF_HAA", 0, 9, (const double[]){2.61799387799, 0, 0}, (const double[]){0.2999, 0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LF_HAA", 9, 10, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[7], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LF_HIP_LF_hip_fixed", 10, 11, (const double[]){-2.61799387799, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LF_hip_fixed_LF_HFE", 11, 12, (const double[]){0, 0, 1.57079632679}, (const double[]){0.0599, 0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LF_HFE", 12, 13, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[8], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LF_THIGH_LF_thigh_fixed", 13, 14, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LF_thigh_fixed_LF_KFE", 14, 15, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LF_KFE", 15, 16, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[9], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LF_shank_LF_shank_fixed", 16, 17, (const double[]){0, 0, -1.5707963267}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LF_shank_fixed_LF_FOOT", 17, 18, (const double[]){0, 0, 0}, (const double[]){0.08795, 0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));
    
    joints_.push_back(new Joint(joints_.size(), -1, "base_RF_HAA", 0, 19, (const double[]){-2.61799387799, 0, 0}, (const double[]){0.2999, -0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RF_HAA", 19, 20, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[10], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RF_HIP_RF_hip_fixed", 20, 21, (const double[]){2.61799387799, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RF_hip_fixed_RF_HFE", 21, 22, (const double[]){0, 0, -1.57079632679}, (const double[]){0.0599, -0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RF_HFE", 22, 23, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[11], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RF_THIGH_RF_thigh_fixed", 23, 24, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RF_thigh_fixed_RF_KFE", 24, 25, (const double[]){0, 0, -1.57079632679}, (const double[]){0, -0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RF_KFE", 25, 26, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[12], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RF_shank_RF_shank_fixed", 26, 27, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RF_shank_fixed_RF_FOOT", 27, 28, (const double[]){0, 0, 0}, (const double[]){0.08795, -0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));
    
    joints_.push_back(new Joint(joints_.size(), -1, "base_LH_HAA", 0, 29, (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){-0.2999, 0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LH_HAA", 29, 30, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[13], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LH_HIP_LH_hip_fixed", 30, 31, (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LH_hip_fixed_LH_HFE", 31, 32, (const double[]){0, 0, 1.57079632679}, (const double[]){-0.0599, 0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LH_HFE", 32, 33, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[14], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LH_THIGH_LH_thigh_fixed", 33, 34, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LH_thigh_fixed_LH_KFE", 34, 35, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LH_KFE", 35, 36, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[15], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LH_shank_LH_shank_fixed", 36, 37, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LH_shank_fixed_LH_FOOT", 37, 38, (const double[]){0, 0, 0}, (const double[]){-0.08795, 0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));

    joints_.push_back(new Joint(joints_.size(), -1, "base_RH_HAA", 0, 39, (const double[]){2.61799387799, 0, -3.14159265359}, (const double[]){-0.2999, -0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RH_HAA", 39, 40, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[16], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RH_HIP_RH_hip_fixed", 40, 41, (const double[]){2.61799387799, 0, -3.14159265359}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RH_hip_fixed_RH_HFE", 41, 42, (const double[]){0, 0, -1.57079632679}, (const double[]){-0.0599, -0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RH_HFE", 42, 43, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[17], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RH_THIGH_RH_thigh_fixed", 43, 44, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RH_thigh_fixed_RH_KFE", 44, 45, (const double[]){0, 0, -1.57079632679}, (const double[]){0, -0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RH_KFE", 45, 46, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[18], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RH_shank_RH_shank_fixed", 46, 47, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RH_shank_fixed_RH_FOOT", 47, 48, (const double[]){0, 0, 0}, (const double[]){-0.08795, -0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));

    num_of_moving_joints_ = moving_joints_.size();
  }

  Eigen::Matrix3d VecToSkew(const Eigen::Vector3d& v) {
      Eigen::Matrix3d skew;
      skew <<     0,    -v.z(),     v.y(),
              v.z(),         0,    -v.x(),
             -v.y(),     v.x(),         0;
      return skew;
  }

  void MakeCompositeBody() {
    for (int i = 0; i < joints_.size(); i++){
      Joint* joint = joints_[i];

      if (joint->P_c_ != Eigen::Vector3d::Zero()) continue;

      // consider only fixed joint
      Link* parent_link = links_[joint->parent_link_idx_];
      Link* child_link = links_[joint->child_link_idx_];

      Eigen::Vector3d r_p1o1_p1_ = parent_link->r_po_p_;
      Eigen::Vector3d r_p2o2_p2_ = child_link->r_po_p_;
      Eigen::Vector3d r_p1o2_p1_ = joint->r_pc_p_ + joint->R_pc_ * joint->R_ccp_ * r_p2o2_p2_;
      Eigen::Vector3d r_p1o3_p1_ = (parent_link->mass_ * r_p1o1_p1_ + child_link->mass_ * r_p1o2_p1_) / (parent_link->mass_ + child_link->mass_);

      Eigen::Vector3d r_o3o1_p1_ = r_p1o1_p1_ - r_p1o3_p1_;
      Eigen::Vector3d r_o3o2_p1_ = r_p1o2_p1_ - r_p1o3_p1_;

      Eigen::Matrix3d inertia_o1_o1_ = parent_link->inertia_o_;
      Eigen::Matrix3d inertia_o2_o2_ = child_link->inertia_o_;
      Eigen::Matrix3d inertia_o1_p1_ = parent_link->R_po_* inertia_o1_o1_ * parent_link->R_po_.transpose();
      Eigen::Matrix3d inertia_o2_p1_ = (joint->R_pc_ * joint->R_ccp_ * child_link->R_po_) * inertia_o2_o2_ * (joint->R_pc_ * joint->R_ccp_ * child_link->R_po_).transpose();
      Eigen::Matrix3d inertia_o3_p1_ = inertia_o1_p1_ + inertia_o2_p1_ - parent_link->mass_ * VecToSkew(r_o3o1_p1_) * VecToSkew(r_o3o1_p1_) - child_link->mass_ * VecToSkew(r_o3o2_p1_) * VecToSkew(r_o3o2_p1_);
      // inertia_o3_p1_ is ame to inertia_o3_o3_

      Eigen::Matrix3d R_po = Eigen::Matrix3d::Identity();

      int new_link_idx = links_.size();

      Link* new_link = new Link(new_link_idx, parent_link->name_ + "_" + child_link->name_, parent_link->parent_joint_idx_, R_po, r_p1o3_p1_, parent_link->mass_ + child_link->mass_, inertia_o3_p1_);
      links_.push_back(new_link);

      // update other joints
      for (int j = 0; j < joints_.size(); j++){
        if (joints_[j]->parent_link_idx_ == child_link->idx_)
        {
          joints_[j]->parent_link_idx_ = new_link_idx;
          joints_[j]->R_pc_ = joint->R_pc_ * joints_[j]->R_pc_;
          joints_[j]->r_pc_p_ = joint->r_pc_p_ + joint->R_pc_ * joints_[j]->r_pc_p_;
        }  
        else if (joints_[j]->parent_link_idx_ == parent_link->idx_)
        {
          joints_[j]->parent_link_idx_ = new_link_idx;
        }  
        else if (joints_[j]->child_link_idx_ == parent_link->idx_)
        {
          joints_[j]->child_link_idx_ = new_link_idx;

        }        
      }

      // update root
      if (new_link->name_.substr(0, 5) == "base_") root = new_link;
    }

    // for (int i = 0; i < moving_joints_.size(); i++){
    //   Joint* joint = moving_joints_[i];
    //   Link* parent_link = links_[joint->parent_link_idx_];
    //   Link* child_link = links_[joint->child_link_idx_];

    //   cout << "joint " << joint->name_ << endl;
    //   cout << "joint R_pc " << joint->R_pc_ << endl;
    //   cout << "joint r_pc_p " << joint->r_pc_p_ << endl;
    //   cout << "parent " << parent_link->name_ << endl;
    //   cout << "parent parent joint idx " << parent_link->parent_joint_idx_ << endl;
    //   cout << "parent r_po_p " << parent_link->r_po_p_ << endl;
    //   cout << "parent R_po " << parent_link->R_po_ << endl;
    //   cout << "parent mass " << parent_link->mass_ << endl;
    //   cout << "child " << child_link->name_ << endl;
    //   cout << "child parent joint idx " << child_link->parent_joint_idx_ << endl;
    //   cout << "child r_po_p " << child_link->r_po_p_ << endl;
    //   cout << "child R_po " << child_link->R_po_ << endl;
    //   cout << "child mass " << child_link->mass_ << endl << endl;
    // }
    // exit(1);
  }

  void ForwardKinematics() {
    // set root
    root->R_wp_ = R_wb_;
    root->r_po_w_ = root->R_wp_ * root->r_po_p_;
    root->r_wo_w_ = r_wb_w_ + root->r_po_w_ ;
    root->inertia_w_ = (root->R_wp_ * root->R_po_) * root->inertia_o_ * (root->R_wp_ * root->R_po_).transpose();

    // set joints and links
    for (int i = 0; i < moving_joints_.size(); i++){
      Joint* joint = moving_joints_[i];
      Link* parent_link = links_[joint->parent_link_idx_];
      Link* child_link = links_[joint->child_link_idx_];

      // joint
      joint->R_wp_ = parent_link->R_wp_;
      joint->P_w_ = joint->R_wp_ * joint->R_pc_ * joint->P_c_;
      joint->S_w_.segment(3, 3) = joint->P_w_;
      joint->r_pc_w_ = joint->R_wp_ * joint->r_pc_p_;

      joint->r_wc_w_ = parent_link->r_wo_w_ - parent_link->r_po_w_ + joint->r_pc_w_;

      if (parent_link->parent_joint_idx_ == -1){ // it is connected to root
        // velocity, omega
        joint->v_c_w_ = v_b_w_ + omega_b_w_.cross(joint->r_pc_w_);
        joint->omega_c_w_ = omega_b_w_;
        
        joint->v_cp_w_ = joint->v_c_w_ + (joint->S_w_ * gv_[i]).segment(0, 3);
        joint->omega_cp_w_ = joint->omega_c_w_ + (joint->S_w_ * gv_[i]).segment(3, 3);

        // acceleration, alpha
        joint->a_c_w_ = a_b_w_ + alpha_b_w_.cross(joint->r_pc_w_) + VecToSkew(omega_b_w_) * VecToSkew(omega_b_w_) * joint->r_pc_w_;
        joint->alpha_c_w_ = alpha_b_w_;
        joint->a_cp_w_ = joint->a_c_w_; // TODO: assuem only revolute joint
        joint->alpha_cp_w_ = joint->alpha_c_w_ + joint->omega_c_w_.cross(joint->P_w_) * gv_[i]; // TODO: assuem only revolute joint
      }
      else {
        Joint* parent_joint = joints_[parent_link->parent_joint_idx_];

        // velocity, omega
        joint->v_c_w_ = parent_joint->v_cp_w_ + parent_joint->omega_cp_w_.cross(joint->r_pc_w_);
        joint->omega_c_w_ = parent_joint->omega_cp_w_;
        joint->v_cp_w_ = joint->v_c_w_ + (joint->S_w_ * gv_[i]).segment(0, 3);
        joint->omega_cp_w_ = joint->omega_c_w_ + (joint->S_w_ * gv_[i]).segment(3, 3);

        // acceleration, alpha
        joint->a_c_w_ = parent_joint->a_cp_w_ + parent_joint->alpha_cp_w_.cross(joint->r_pc_w_) + VecToSkew(parent_joint->omega_cp_w_) * VecToSkew(parent_joint->omega_cp_w_) * joint->r_pc_w_;
        joint->alpha_c_w_ = parent_joint->alpha_cp_w_;
        joint->a_cp_w_ = joint->a_c_w_; // TODO: assuem only revolute joint
        joint->alpha_cp_w_ = joint->alpha_c_w_ + joint->omega_c_w_.cross(joint->P_w_) * gv_[i]; // TODO: assuem only revolute joint
      }

      // link
      child_link->R_wp_ = joint->R_wp_ * joint->R_pc_ * joint->R_ccp_;
      child_link->r_po_w_ = child_link->R_wp_ * child_link->r_po_p_;
      child_link->r_wo_w_ = joint->r_wc_w_ + child_link->r_po_w_ ;
      child_link->inertia_w_ = (child_link->R_wp_ * child_link->R_po_) * child_link->inertia_o_ * (child_link->R_wp_ * child_link->R_po_).transpose();
    }

    // for (int i = 0; i < moving_joints_.size(); i++){
    //   Joint* joint = moving_joints_[i];
    //   cout << "joint " << joint->name_ << endl;
    //   cout << "velocity " << joint->v_cp_w_ << endl;
    //   cout << "omega " << joint->omega_cp_w_ << endl;
    //   cout << "acceleration " << joint->a_cp_w_ << endl;
    //   cout << "alpha " << joint->alpha_cp_w_ << endl << endl;
    // }
  }

  void CalMassMatrix() {
    M_ = Eigen::MatrixXd::Zero(6 + num_of_moving_joints_, 6 + num_of_moving_joints_);
    double m_comp[num_of_moving_joints_][num_of_moving_joints_] = {0,};
    Eigen::Vector3d r_po_w_comp[num_of_moving_joints_][num_of_moving_joints_] = {Eigen::Vector3d::Zero(),};
    Eigen::Matrix3d I_o_w_comp[num_of_moving_joints_][num_of_moving_joints_] = {Eigen::Matrix3d::Zero(),};

    int leaf = moving_joints_.size() - 1;

    for (int j = moving_joints_.size() - 1; j >= 0; j--){
      Joint* joint_j = moving_joints_[j];
      Link* parent_link = links_[joint_j->parent_link_idx_];
      Link* child_link = links_[joint_j->child_link_idx_];

      // cout << "j " << j << endl;
      // cout << "parent_link " << parent_link->name_ << endl;
      // cout << "parent_link parent_joint_idx " << parent_link->parent_joint_idx_ << endl;
      // cout << "child_link " << child_link->name_ << endl;

      // calculate composite body properties
      if (j != moving_joints_.size() - 1){
        int k = j + 1;
        Joint* joint_k = moving_joints_[k];
        Link* parent_link_k = links_[joint_k->parent_link_idx_];
        Link* child_link_k = links_[joint_k->child_link_idx_];

        if (parent_link_k->idx_ != child_link->idx_){ // j is leaf joint
          leaf = j;
        }
        else { // j is not leaf joint
          m_comp[j][leaf] = child_link->mass_ + m_comp[k][leaf];
          r_po_w_comp[j][leaf] = (child_link->mass_ * (child_link->r_po_w_) + m_comp[k][leaf] * (joint_k->r_pc_w_ + r_po_w_comp[k][leaf]))/ m_comp[j][leaf];

          Eigen::Vector3d r_1 = child_link->r_po_w_ - r_po_w_comp[j][leaf];
          Eigen::Vector3d r_2 = (joint_k->r_pc_w_ + r_po_w_comp[k][leaf]) - r_po_w_comp[j][leaf];

          I_o_w_comp[j][leaf] = child_link->inertia_w_ + I_o_w_comp[k][leaf] - child_link->mass_ * VecToSkew(r_1) * VecToSkew(r_1) - m_comp[k][leaf] * VecToSkew(r_2) * VecToSkew(r_2);
        }
      }

      if (j == leaf){ // j is leaf joint
        m_comp[j][j] = child_link->mass_;
        r_po_w_comp[j][j] = child_link->r_po_w_;
        I_o_w_comp[j][j] = child_link->inertia_w_;
      }

      // calculate M_jj
      Eigen::Vector3d J_p = joint_j->P_w_.cross(r_po_w_comp[j][leaf]);
      Eigen::Vector3d J_a = joint_j->P_w_;

      Eigen::VectorXd S_j_ = Eigen::VectorXd::Zero(6);
      S_j_.block<3,1>(3, 0) = joint_j->P_w_;

      M_.block<1,1>(6 + j, 6 + j) = (J_p.transpose() * m_comp[j][leaf] * J_p + J_a.transpose() * I_o_w_comp[j][leaf] * J_a);

      Eigen::MatrixXd r = Eigen::MatrixXd::Identity(6, 6);
      Eigen::MatrixXd I = Eigen::MatrixXd::Identity(6, 6);
      while (parent_link->parent_joint_idx_ != -1){
        Joint* joint_i = joints_[parent_link->parent_joint_idx_];
        parent_link = links_[joint_i->parent_link_idx_];
        child_link = links_[joint_i->child_link_idx_];

        r.block<3,3>(0, 3) = -VecToSkew(joint_j->r_wc_w_ - joint_i->r_wc_w_ );
        
        I.block<3,3>(0, 0) = m_comp[j][leaf] * I.block<3,3>(0, 0);
        I.block<3,3>(0, 3) = -m_comp[j][leaf] * VecToSkew(r_po_w_comp[j][leaf]);
        I.block<3,3>(3, 0) = m_comp[j][leaf] * VecToSkew(r_po_w_comp[j][leaf]);
        I.block<3,3>(3, 3) = I_o_w_comp[j][leaf] - m_comp[j][leaf] * VecToSkew(r_po_w_comp[j][leaf]) * VecToSkew(r_po_w_comp[j][leaf]);

        M_.block<1,1>(6 + j, 6 + joint_i->moving_idx_) = S_j_.transpose() * I * r * joint_i->S_w_;
        M_.block<1,1>(6 + joint_i->moving_idx_, 6 + j) = M_.block<1,1>(6 + j, 6 + joint_i->moving_idx_);
      }    

      // calculate M_base_leg coupling      
      I.block<3,3>(0, 0) = m_comp[j][leaf] * I.block<3,3>(0, 0);
      I.block<3,3>(0, 3) = -m_comp[j][leaf] * VecToSkew(r_po_w_comp[j][leaf]);
      I.block<3,3>(3, 0) = m_comp[j][leaf] * VecToSkew(r_po_w_comp[j][leaf]);
      I.block<3,3>(3, 3) = I_o_w_comp[j][leaf] - m_comp[j][leaf] * VecToSkew(r_po_w_comp[j][leaf]) * VecToSkew(r_po_w_comp[j][leaf]);

      r.block<3,3>(0, 3) = -VecToSkew(joint_j->r_wc_w_ - r_wb_w_);
      M_.block<1,6>(6 + j, 0) = S_j_.transpose() * I * r;
      M_.block<6,1>(0, 6 + j) = M_.block<1,6>(6 + j, 0);
    }

    // calculate base mass matrix
    double m_base = root->mass_;
    Eigen::Vector3d r_po_w_base = root->r_po_w_; 
    Eigen::Matrix3d I_o_w_base = root->inertia_w_;

    vector<int> leg_moving_joint_idx[2] = {{0, 3, 6, 9}, {2, 5, 8, 11}};

    for (int i = 0; i < leg_moving_joint_idx[0].size(); i++){
      int leg_start_idx = leg_moving_joint_idx[0][i];
      int leg_end_idx = leg_moving_joint_idx[1][i];

      Eigen::Vector3d new_r_po_w_base = (m_base * r_po_w_base + m_comp[leg_start_idx][leg_end_idx] * (moving_joints_[leg_start_idx]->r_pc_w_ + r_po_w_comp[leg_start_idx][leg_end_idx])) / (m_base + m_comp[leg_start_idx][leg_end_idx]);
      double new_m_base = m_base + m_comp[leg_start_idx][leg_end_idx];

      Eigen::Vector3d r_1 = r_po_w_base - new_r_po_w_base;
      Eigen::Vector3d r_2 = (moving_joints_[leg_start_idx]->r_pc_w_ + r_po_w_comp[leg_start_idx][leg_end_idx]) - new_r_po_w_base;

      I_o_w_base = I_o_w_base + I_o_w_comp[leg_start_idx][leg_end_idx] - m_base * VecToSkew(r_1) * VecToSkew(r_1) - m_comp[leg_start_idx][leg_end_idx] * VecToSkew(r_2) * VecToSkew(r_2);
      m_base = new_m_base;
      r_po_w_base = new_r_po_w_base;
    }

    Eigen::MatrixXd I_base_ = Eigen::MatrixXd::Identity(6, 6);
    I_base_.block<3,3>(0, 0) = m_base * I_base_.block<3,3>(0, 0);
    I_base_.block<3,3>(0, 3) = -m_base * VecToSkew(r_po_w_base);
    I_base_.block<3,3>(3, 0) = m_base * VecToSkew(r_po_w_base);
    I_base_.block<3,3>(3, 3) = I_o_w_base - m_base * VecToSkew(r_po_w_base) * VecToSkew(r_po_w_base);

    M_.block<6,6>(0, 0) = I_base_;
  }

  void CalNonlinearities() {
    non_linearities_ = Eigen::VectorXd::Zero(6 + num_of_moving_joints_);
    vector<Eigen::VectorXd> FT(num_of_moving_joints_, Eigen::VectorXd::Zero(6));
    Eigen::VectorXd b_j = Eigen::VectorXd::Zero(6);

    int leaf = moving_joints_.size() - 1;

    for (int j = moving_joints_.size() - 1; j >= 0; j--){
      Joint* joint_j = moving_joints_[j];
      Link* parent_link = links_[joint_j->parent_link_idx_];
      Link* child_link = links_[joint_j->child_link_idx_];

      // check it is leaf or not
      if (j != moving_joints_.size() - 1){
        int k = j + 1;
        Joint* joint_k = moving_joints_[k];
        Link* parent_link_k = links_[joint_k->parent_link_idx_];
        Link* child_link_k = links_[joint_k->child_link_idx_];

        if (parent_link_k->idx_ != child_link->idx_){ // j is leaf joint
          leaf = j;
        }
      }

      // calculate b
      b_j.segment(0, 3) = child_link->mass_ * VecToSkew(joint_j->omega_cp_w_) * VecToSkew(joint_j->omega_cp_w_) * child_link->r_po_w_;
      b_j.segment(3, 3) = VecToSkew(joint_j->omega_cp_w_) * (child_link->inertia_w_ - child_link->mass_ * VecToSkew(child_link->r_po_w_) * VecToSkew(child_link->r_po_w_)) * joint_j->omega_cp_w_;

      Eigen::Matrix<double, 6, 1> FT_by_gravity;
      FT_by_gravity << child_link->mass_ * gravity_, VecToSkew(child_link->r_po_w_) * (child_link->mass_ * gravity_);

      b_j = b_j - FT_by_gravity;

      // calculate M_j
      Eigen::MatrixXd M_j = Eigen::MatrixXd::Zero(6, 6);
      M_j.block<3,3>(0, 0) = child_link->mass_ * Eigen::MatrixXd::Identity(3, 3);
      M_j.block<3,3>(0, 3) = -child_link->mass_ * VecToSkew(child_link->r_po_w_);
      M_j.block<3,3>(3, 0) = child_link->mass_ * VecToSkew(child_link->r_po_w_);
      M_j.block<3,3>(3, 3) = child_link->inertia_w_ - child_link->mass_ * VecToSkew(child_link->r_po_w_) * VecToSkew(child_link->r_po_w_);

      // make twist_dot_j
      Eigen::Matrix<double, 6, 1> twist_dot_j;
      twist_dot_j << joint_j->a_cp_w_ , joint_j->alpha_cp_w_;

      if (j == leaf){
        FT[j] = M_j * twist_dot_j + b_j;
        non_linearities_.segment(6 + j, 1) = joint_j->S_w_.transpose() * FT[j];
      }
      else {
        Eigen::Matrix<double, 6, 1> FT_by_child;
        FT_by_child << FT[j + 1].segment(0, 3), FT[j + 1].segment(3, 3) + VecToSkew(moving_joints_[j + 1]->r_pc_w_) * (FT[j + 1].segment(0, 3));
        twist_dot_j << joint_j->a_cp_w_ , joint_j->alpha_cp_w_;
        FT[j] = M_j * twist_dot_j + b_j + FT_by_child;
        non_linearities_.segment(6 + j, 1) = joint_j->S_w_.transpose() * FT[j];
      }
    }

    // calculate root
    // calculate b_j
    b_j.segment(0, 3) = root->mass_ * VecToSkew(omega_b_w_) * VecToSkew(omega_b_w_) * root->r_po_w_;
    b_j.segment(3, 3) = VecToSkew(omega_b_w_) * (root->inertia_w_ - root->mass_ * VecToSkew(root->r_po_w_) * VecToSkew(root->r_po_w_)) * omega_b_w_;
    Eigen::Matrix<double, 6, 1> FT_by_gravity;
    FT_by_gravity << root->mass_ * gravity_, VecToSkew(root->r_po_w_) * (root->mass_ * gravity_);

    b_j = b_j - FT_by_gravity;

    // calculate M_j
    Eigen::MatrixXd M_j = Eigen::MatrixXd::Zero(6, 6);
    M_j.block<3,3>(0, 0) = root->mass_ * Eigen::MatrixXd::Identity(3, 3);
    M_j.block<3,3>(0, 3) = -root->mass_ * VecToSkew(root->r_po_w_);
    M_j.block<3,3>(3, 0) = root->mass_ * VecToSkew(root->r_po_w_);
    M_j.block<3,3>(3, 3) = root->inertia_w_ - root->mass_ * VecToSkew(root->r_po_w_) * VecToSkew(root->r_po_w_);

    Eigen::Matrix<double, 6, 1> FT_root = b_j;

    Eigen::Matrix<double, 6, 1> FT_by_child;
    vector<int> leg_moving_joint_idx[2] = {{0, 3, 6, 9}, {2, 5, 8, 11}};
    for (int i = 0; i < leg_moving_joint_idx[0].size(); i++){
      int leg_start_idx = leg_moving_joint_idx[0][i];
      FT_by_child << FT[leg_start_idx].segment(0, 3), FT[leg_start_idx].segment(3, 3) + VecToSkew(moving_joints_[leg_start_idx]->r_pc_w_) * (FT[leg_start_idx].segment(0, 3));
      FT_root += FT_by_child;
    }
    non_linearities_.segment(0, 6) = FT_root;
  }
  
  Eigen::MatrixXd GetMassMatrix(){
    CalMassMatrix();
    return M_;
  }

  Eigen::MatrixXd GetNonlinearities(){
    CalNonlinearities();
    return non_linearities_;
  }
};
/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  ArticulatedSystem anymal(gc, gv);
  anymal.MakeCompositeBody();
  anymal.ForwardKinematics();

  return anymal.GetNonlinearities();
}

