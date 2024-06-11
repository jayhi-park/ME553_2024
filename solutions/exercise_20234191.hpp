#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <set>
#include <vector>
#include <string.h>
#include <iostream>

using namespace std;

enum JointType
{
  REVOLUTE,
  PRISMATIC,
  FIXED
};

class Link
{
  public:
    int idx_;
    string name_;
    int parent_joint_idx_ = -1;
    vector<int> children_joint_idx_;

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

    // mass matrix & non-linearities for single rigid body
    Eigen::MatrixXd M_single_ = Eigen::MatrixXd::Zero(6, 6);
    Eigen::VectorXd b_single_ = Eigen::VectorXd::Zero(6);

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
    JointType joint_type_;

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
    Eigen::VectorXd S_dot_w_ = Eigen::VectorXd::Zero(6);    

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
    Joint(int idx, int moving_idx, string name, int parent_link_idx, int child_link_idx, JointType joint_type, const double rpy[], const double xyz[], double theta, const double p[]) : 
      idx_(idx), moving_idx_(moving_idx), name_(name), parent_link_idx_(parent_link_idx), child_link_idx_(child_link_idx), joint_type_(joint_type){

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
  Joint* root_joint; // base
  Link* root; // body

  // world
  Eigen::VectorXd gf_root_w_;

  // joint and link
  vector<Joint*> joints_;
  vector<Link*> links_;

  // generalized coordinate, velocity, force and acceleration
  Eigen::VectorXd gc_; // input
  Eigen::VectorXd gv_; // input
  Eigen::VectorXd gf_; // input
  Eigen::VectorXd ga_; // have to calculate

  // moving joints list
  vector<Joint*> moving_joints_;

  // num_of_moving_joints
  int num_of_moving_joints_;

  // mass matrix for composite body
  vector<vector<double>> m_comp_;
  vector<vector<Eigen::Vector3d>> r_po_w_comp_;
  vector<vector<Eigen::Matrix3d>> I_o_w_comp_;
  Eigen::MatrixXd M_com_;

  // non_linearities using RNE
  Eigen::VectorXd FT_root_;
  vector<Eigen::VectorXd> FT_;
  Eigen::VectorXd b_;

  // mass matrix & non_linearities for single rigid body
  Eigen::MatrixXd M_root_single_ = Eigen::MatrixXd::Zero(6, 6);
  Eigen::VectorXd b_root_single_ = Eigen::VectorXd::Zero(6);

  // mass matrix & non-linearities for articulated body
  Eigen::MatrixXd M_root_a_ = Eigen::MatrixXd::Zero(6, 6);
  Eigen::VectorXd b_root_a_ = Eigen::VectorXd::Zero(6);
  std::vector<Eigen::MatrixXd> M_a_;
  vector<Eigen::VectorXd> b_a_;

  // gravity
  Eigen::Vector3d gravity_ = Eigen::Vector3d(0, 0, -9.81);  

public:
  ArticulatedSystem(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf)
  {
    // input
    gf_root_w_ = gf.segment(0, 6);
    gc_ = gc.segment(6, gc.size() - 6);
    gv_ = gv.segment(6, gv.size() - 6);
    gf_ = gf.segment(6, gf.size() - 6);

    // set root
    root_joint = new Joint(0, -1, "base_joint", -1, 0, FIXED, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0});
    joints_.push_back(root_joint);
    root = new Link(0, "base", 0, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.0, (const double[]){0, 0, 0, 0, 0, 0});
    links_.push_back(root);

    root_joint->r_pc_w_ = gc.segment(0, 3);
    root_joint->r_wc_w_ = root_joint->r_pc_w_;
    root_joint->v_c_w_ = gv.segment(0, 3);
    root_joint->omega_c_w_ = gv.segment(3, 3);
    root_joint->a_c_w_ = Eigen::Vector3d::Zero();
    root_joint->alpha_c_w_ = Eigen::Vector3d::Zero();
    root_joint->v_cp_w_ = root_joint->v_c_w_;
    root_joint->omega_cp_w_ = root_joint->omega_c_w_;
    root_joint->a_cp_w_ = root_joint->a_c_w_;
    root_joint->alpha_c_w_ = Eigen::Vector3d::Zero();
    root_joint->alpha_cp_w_ = root_joint->alpha_c_w_;

    Eigen::Quaternion<double> q(gc[3], gc[4], gc[5], gc[6]);
    root_joint->R_wp_ = q.normalized().toRotationMatrix();

    // set links
    links_.push_back(new Link(links_.size(), "base_inertia", 1, (const double[]){0, 0, 0}, (const double[]){-0.018, -0.002, 0.024}, 6.222, (const double[]){0.017938806, 0.00387963, 0.001500772, 0.370887745, 6.8963e-05, 0.372497653}));
    links_.push_back(new Link(links_.size(), "face_front", 2, (const double[]){0, 0, 0}, (const double[]){0.042, -0.001, 0.004}, 0.73, (const double[]){0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938}));
    links_.push_back(new Link(links_.size(), "face_rear", 3, (const double[]){0, 0, 0}, (const double[]){0.042, -0.001, 0.004}, 0.73, (const double[]){0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938}));
    links_.push_back(new Link(links_.size(), "battery", 4, (const double[]){0, 0, 0}, (const double[]){-0.00067, -0.00023, -0.03362}, 5.53425, (const double[]){0.00749474794, 0.00016686282, 7.82763e-05, 0.0722338913, 1.42902e-06, 0.07482717535}));
    links_.push_back(new Link(links_.size(), "docking_hatch_cover", 5, (const double[]){0, 0, 0}, (const double[]){-0.003, 0.0, 0.005}, 0.065, (const double[]){0.00063283, 0.0, 3.45e-07, 0.00110971, 0.0, 0.00171883}));
    links_.push_back(new Link(links_.size(), "lidar_cage", 6, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0, 0, 0, 0}));
    links_.push_back(new Link(links_.size(), "lidar", 7, (const double[]){0, 0, 0}, (const double[]){-0.012, 0.001, -0.008}, 0.695, (const double[]){0.000846765, 6.9565e-05, 0.00027111, 0.001367583, 5.8984e-05, 0.001363673}));
    links_.push_back(new Link(links_.size(), "hatch", 8, (const double[]){0, 0, 0}, (const double[]){0.116, 0.0, 0.0758}, 0.142, (const double[]){0.001, 0.001, 0.001, 0.001, 0.001, 0.001}));
    
    links_.push_back(new Link(links_.size(), "LF_HAA", 9, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LF_HIP", 10, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LF_hip_fixed", 11, (const double[]){0, 0, 0}, (const double[]){0.048, 0.008, -0.003}, 0.74, (const double[]){0.001393106, 8.4012e-05, 2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "LF_HFE", 12, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LF_THIGH", 13, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LF_thigh_fixed", 14, (const double[]){0, 0, 0}, (const double[]){0.0, 0.018, -0.169}, 1.03, (const double[]){0.018644469, 5.2e-08, 1.0157e-05, 0.019312599, 0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "LF_KFE", 15, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LF_SHANK", 16, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LF_shank_fixed", 17, (const double[]){0, 0, 0}, (const double[]){0.03463, 0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, 2.142561e-05, 1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "LF_FOOT", 18, (const double[]){0, 0, 0}, (const double[]){0.00948, -0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, 2.63048e-06, 6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05}));

    links_.push_back(new Link(links_.size(), "RF_HAA", 19, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RF_HIP", 20, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RF_hip_fixed", 21, (const double[]){0, 0, 0}, (const double[]){0.048, -0.008, -0.0030}, 0.74, (const double[]){0.001393106, -8.4012e-05, 2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "RF_HFE", 22, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RF_THIGH", 23, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RF_thigh_fixed", 24, (const double[]){0, 0, 0}, (const double[]){0.0, -0.018, -0.169}, 1.03, (const double[]){0.018644469, -5.2e-08, 1.0157e-05, 0.019312599, -0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "RF_KFE", 25, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RF_SHANK", 26, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RF_shank_fixed", 27, (const double[]){0, 0, 0}, (const double[]){0.03463, -0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, -2.142561e-05, 1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "RF_FOOT", 28, (const double[]){0, 0, 0}, (const double[]){0.00948, 0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, -2.63048e-06, 6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05}));

    links_.push_back(new Link(links_.size(), "LH_HAA", 29, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LH_HIP", 30, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LH_hip_fixed", 31, (const double[]){0, 0, 0}, (const double[]){-0.048, 0.008, -0.003}, 0.74, (const double[]){0.001393106, -8.4012e-05, -2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "LH_HFE", 32, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LH_THIGH", 33, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LH_thigh_fixed", 34, (const double[]){0, 0, 0}, (const double[]){0.0, 0.018, -0.169}, 1.03, (const double[]){0.018644469, -5.2e-08, -1.0157e-05, 0.019312599, 0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "LH_KFE", 35, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "LH_SHANK", 36, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "LH_shank_fixed", 37, (const double[]){0, 0, 0}, (const double[]){-0.03463, 0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, -2.142561e-05, -1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "LH_FOOT", 38, (const double[]){0, 0, 0}, (const double[]){-0.00948, -0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, -2.63048e-06, -6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05}));

    links_.push_back(new Link(links_.size(), "RH_HAA", 39, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RH_HIP", 40, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RH_hip_fixed", 41, (const double[]){0, 0, 0}, (const double[]){-0.048, -0.008, -0.003}, 0.74, (const double[]){0.001393106, 8.4012e-05, -2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509}));
    links_.push_back(new Link(links_.size(), "RH_HFE", 42, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RH_THIGH", 43, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RH_thigh_fixed", 44, (const double[]){0, 0, 0}, (const double[]){0.0, -0.018, -0.169}, 1.03, (const double[]){0.018644469, 5.2e-08, -1.0157e-05, 0.019312599, -0.002520077, 0.002838361}));
    links_.push_back(new Link(links_.size(), "RH_KFE", 45, (const double[]){0, 0, 0}, (const double[]){-0.063, 7e-05, 0.00046}, 2.04, (const double[]){0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827}));
    links_.push_back(new Link(links_.size(), "RH_SHANK", 46, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0.001, (const double[]){0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001}));
    links_.push_back(new Link(links_.size(), "RH_shank_fixed", 47, (const double[]){0, 0, 0}, (const double[]){-0.03463, -0.00688, 0.00098}, 0.33742, (const double[]){0.00032748005, 2.142561e-05, -1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521}));
    links_.push_back(new Link(links_.size(), "RH_FOOT", 48, (const double[]){0, 0, 0}, (const double[]){-0.00948, 0.00948, 0.1468}, 0.25, (const double[]){0.00317174097, 2.63048e-06, -6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05}));

    // set joints
    joints_.push_back(new Joint(joints_.size(), -1, "base_to_base_inertia", 0, 1, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_face_front", 0, 2, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){0.4145, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_face_rear", 0, 3, JointType::FIXED, (const double[]){0, 0, 3.14159265359}, (const double[]){-0.4145, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_battery", 0, 4, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_to_docking_hatch_cover", 0, 5, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){0.343, 0.0, -0.07}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_to_lidar_cage", 0, 6, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){-0.364, 0.0, 0.0735}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "lidar_cage_to_lidar", 6, 7, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0.0687}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "base_hatch", 0, 8, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));

    joints_.push_back(new Joint(joints_.size(), -1, "base_LF_HAA", 0, 9, JointType::FIXED, (const double[]){2.61799387799, 0, 0}, (const double[]){0.2999, 0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LF_HAA", 9, 10, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[7], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LF_HIP_LF_hip_fixed", 10, 11, JointType::FIXED, (const double[]){-2.61799387799, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LF_hip_fixed_LF_HFE", 11, 12, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){0.0599, 0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LF_HFE", 12, 13, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[8], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LF_THIGH_LF_thigh_fixed", 13, 14, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LF_thigh_fixed_LF_KFE", 14, 15, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LF_KFE", 15, 16, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[9], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LF_shank_LF_shank_fixed", 16, 17, JointType::FIXED, (const double[]){0, 0, -1.5707963267}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LF_shank_fixed_LF_FOOT", 17, 18, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){0.08795, 0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));
    
    joints_.push_back(new Joint(joints_.size(), -1, "base_RF_HAA", 0, 19, JointType::FIXED, (const double[]){-2.61799387799, 0, 0}, (const double[]){0.2999, -0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RF_HAA", 19, 20, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[10], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RF_HIP_RF_hip_fixed", 20, 21, JointType::FIXED, (const double[]){2.61799387799, 0, 0}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RF_hip_fixed_RF_HFE", 21, 22, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){0.0599, -0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RF_HFE", 22, 23, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[11], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RF_THIGH_RF_thigh_fixed", 23, 24, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RF_thigh_fixed_RF_KFE", 24, 25, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){0, -0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RF_KFE", 25, 26, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[12], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RF_shank_RF_shank_fixed", 26, 27, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RF_shank_fixed_RF_FOOT", 27, 28, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){0.08795, -0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));
    
    joints_.push_back(new Joint(joints_.size(), -1, "base_LH_HAA", 0, 29, JointType::FIXED, (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){-0.2999, 0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LH_HAA", 29, 30, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[13], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LH_HIP_LH_hip_fixed", 30, 31, JointType::FIXED, (const double[]){-2.61799387799, 0, -3.14159265359}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LH_hip_fixed_LH_HFE", 31, 32, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){-0.0599, 0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LH_HFE", 32, 33, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[14], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LH_THIGH_LH_thigh_fixed", 33, 34, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LH_thigh_fixed_LH_KFE", 34, 35, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "LH_KFE", 35, 36, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[15], (const double[]){1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "LH_shank_LH_shank_fixed", 36, 37, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "LH_shank_fixed_LH_FOOT", 37, 38, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){-0.08795, 0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));

    joints_.push_back(new Joint(joints_.size(), -1, "base_RH_HAA", 0, 39, JointType::FIXED, (const double[]){2.61799387799, 0, -3.14159265359}, (const double[]){-0.2999, -0.104, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RH_HAA", 39, 40, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[16], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RH_HIP_RH_hip_fixed", 40, 41, JointType::FIXED, (const double[]){2.61799387799, 0, -3.14159265359}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RH_hip_fixed_RH_HFE", 41, 42, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){-0.0599, -0.08381, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RH_HFE", 42, 43, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[17], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RH_THIGH_RH_thigh_fixed", 43, 44, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RH_thigh_fixed_RH_KFE", 44, 45, JointType::FIXED, (const double[]){0, 0, -1.57079632679}, (const double[]){0, -0.1003, -0.285}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), moving_joints_.size(), "RH_KFE", 45, 46, JointType::REVOLUTE, (const double[]){0, 0, 0}, (const double[]){0, 0, 0}, gc[18], (const double[]){-1, 0, 0})); // rev
    moving_joints_.push_back(joints_.back());
    joints_.push_back(new Joint(joints_.size(), -1, "RH_shank_RH_shank_fixed", 46, 47, JointType::FIXED, (const double[]){0, 0, 1.57079632679}, (const double[]){0, 0, 0}, 0, (const double[]){0, 0, 0}));
    joints_.push_back(new Joint(joints_.size(), -1, "RH_shank_fixed_RH_FOOT", 47, 48, JointType::FIXED, (const double[]){0, 0, 0}, (const double[]){-0.08795, -0.01305, -0.33797}, 0, (const double[]){0, 0, 0}));

    num_of_moving_joints_ = moving_joints_.size();

    // initialize
    Initialize();
  }

  Eigen::Matrix3d VecToSkew(const Eigen::Vector3d& v) {
      Eigen::Matrix3d skew;
      skew <<     0,    -v.z(),     v.y(),
              v.z(),         0,    -v.x(),
             -v.y(),     v.x(),         0;
      return skew;
  }

  void Initialize() {
    for (int i = 0; i < num_of_moving_joints_; i++){
      M_a_.push_back(Eigen::MatrixXd::Zero(6, 6));
      b_a_.push_back(Eigen::VectorXd::Zero(6));
    }

    // make composite body. It means to remove fixed joints
    MakeCompositeBody();

    // set children of composite body
    SetChildren();
  }

  void MakeCompositeBody() {
    for (int i = 1; i < joints_.size(); i++){ // don't consider base
      Joint* joint = joints_[i];

      if (joint->joint_type_ != JointType::FIXED) continue;

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

  void SetChildren() {
    for (int i = 0; i < moving_joints_.size(); i++){
      Joint* joint = moving_joints_[i];
      Link* parent_link = links_[joint->parent_link_idx_];

      parent_link->children_joint_idx_.push_back(joint->idx_);
    }
  }

  void RecursiveForwardKinematics(Joint* joint) {
    Link* child_link = links_[joint->child_link_idx_];  

    // cal fk in link
    child_link->R_wp_ = joint->R_wp_ * joint->R_pc_ * joint->R_ccp_;
    child_link->r_po_w_ = child_link->R_wp_ * child_link->r_po_p_;
    child_link->r_wo_w_ = joint->r_wc_w_ + child_link->r_po_w_ ;
    child_link->inertia_w_ = (child_link->R_wp_ * child_link->R_po_) * child_link->inertia_o_ * (child_link->R_wp_ * child_link->R_po_).transpose();

    for (int i = 0; i < child_link->children_joint_idx_.size(); i++){
      Joint* child_joint = joints_[child_link->children_joint_idx_[i]];

      // cal fk in child_joint
      // joint
      child_joint->R_wp_ = child_link->R_wp_;
      child_joint->P_w_ = child_joint->R_wp_ * child_joint->R_pc_ * child_joint->P_c_;
      if (child_joint->joint_type_ == JointType::REVOLUTE){
        child_joint->S_w_.segment(3, 3) = child_joint->P_w_;
      }
      else if (child_joint->joint_type_ == JointType::PRISMATIC) {
        child_joint->S_w_.segment(0, 3) = child_joint->P_w_;
      }
      else { // fixed
        child_joint->S_w_ = Eigen::VectorXd::Zero(6);
      }
      child_joint->r_pc_w_ = child_joint->R_wp_ * child_joint->r_pc_p_;

      child_joint->r_wc_w_ = child_link->r_wo_w_ - child_link->r_po_w_ + child_joint->r_pc_w_;

      // velocity, omega
      child_joint->v_c_w_ = joint->v_cp_w_ + joint->omega_cp_w_.cross(child_joint->r_pc_w_);
      child_joint->omega_c_w_ = joint->omega_cp_w_;
      child_joint->v_cp_w_ = child_joint->v_c_w_ + (child_joint->S_w_ * gv_[child_joint->moving_idx_]).segment(0, 3);
      child_joint->omega_cp_w_ = child_joint->omega_c_w_ + (child_joint->S_w_ * gv_[child_joint->moving_idx_]).segment(3, 3);

      // S_dot
      if (child_joint->joint_type_ == JointType::REVOLUTE)
        child_joint->S_dot_w_.segment(3, 3) = child_joint->omega_c_w_.cross(child_joint->P_w_);
      else if (child_joint->joint_type_ == JointType::PRISMATIC)
        child_joint->S_dot_w_.segment(0, 3) = child_joint->omega_c_w_.cross(child_joint->P_w_);
      else
        child_joint->S_dot_w_ = Eigen::VectorXd::Zero(6);

      // acceleration, alpha
      child_joint->a_c_w_ = joint->a_cp_w_ + joint->alpha_cp_w_.cross(child_joint->r_pc_w_) + VecToSkew(joint->omega_cp_w_) * VecToSkew(joint->omega_cp_w_) * child_joint->r_pc_w_;
      child_joint->alpha_c_w_ = joint->alpha_cp_w_;
      child_joint->a_cp_w_ = child_joint->a_c_w_ + (child_joint->S_dot_w_ * gv_[child_joint->moving_idx_]).segment(0, 3);
      child_joint->alpha_cp_w_ = child_joint->alpha_c_w_ + (child_joint->S_dot_w_ * gv_[child_joint->moving_idx_]).segment(3, 3);

      RecursiveForwardKinematics(child_joint);
    }
  }

  void ForwardKinematics() {
    RecursiveForwardKinematics(root_joint);

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
    M_com_ = Eigen::MatrixXd::Zero(6 + num_of_moving_joints_, 6 + num_of_moving_joints_);
    for (int i = 0; i < num_of_moving_joints_; i++)
    {
      vector<double> m_row;
      vector<Eigen::Vector3d> r_po_w_row;
      vector<Eigen::Matrix3d> I_o_w_row;
      for (int j = 0; j < num_of_moving_joints_; j++)
      {
        m_row.push_back(0.0);
        r_po_w_row.push_back(Eigen::Vector3d::Zero());
        I_o_w_row.push_back(Eigen::Matrix3d::Zero());
      }
      m_comp_.push_back(m_row);
      r_po_w_comp_.push_back(r_po_w_row);
      I_o_w_comp_.push_back(I_o_w_row);
    }

    vector<int> leaf_moving_idx;
    for (int i = 0; i < moving_joints_.size(); i++){
      Joint* joint = moving_joints_[i];
      if (links_[joint->child_link_idx_]->children_joint_idx_.size() == 0){
        leaf_moving_idx.push_back(i);
      }
    }

    int leaf = moving_joints_.size() - 1;

    for (int j = moving_joints_.size() - 1; j >= 0; j--){
      Joint* joint_j = moving_joints_[j];
      Link* parent_link = links_[joint_j->parent_link_idx_];
      Link* child_link = links_[joint_j->child_link_idx_];

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
          m_comp_[j][leaf] = child_link->mass_ + m_comp_[k][leaf];
          r_po_w_comp_[j][leaf] = (child_link->mass_ * (child_link->r_po_w_) + m_comp_[k][leaf] * (joint_k->r_pc_w_ + r_po_w_comp_[k][leaf]))/ m_comp_[j][leaf];

          Eigen::Vector3d r_1 = child_link->r_po_w_ - r_po_w_comp_[j][leaf];
          Eigen::Vector3d r_2 = (joint_k->r_pc_w_ + r_po_w_comp_[k][leaf]) - r_po_w_comp_[j][leaf];

          I_o_w_comp_[j][leaf] = child_link->inertia_w_ + I_o_w_comp_[k][leaf] - child_link->mass_ * VecToSkew(r_1) * VecToSkew(r_1) - m_comp_[k][leaf] * VecToSkew(r_2) * VecToSkew(r_2);
        }
      }

      if (j == leaf){ // j is leaf joint
        m_comp_[j][j] = child_link->mass_;
        r_po_w_comp_[j][j] = child_link->r_po_w_;
        I_o_w_comp_[j][j] = child_link->inertia_w_;
      }

      // calculate M_single
      Eigen::MatrixXd M_single = Eigen::MatrixXd::Zero(6, 6);
      M_single.block<3,3>(0, 0) = m_comp_[j][leaf] * Eigen::MatrixXd::Identity(3, 3);
      M_single.block<3,3>(0, 3) = -m_comp_[j][leaf] * VecToSkew(r_po_w_comp_[j][leaf]);
      M_single.block<3,3>(3, 0) = m_comp_[j][leaf] * VecToSkew(r_po_w_comp_[j][leaf]);
      M_single.block<3,3>(3, 3) = I_o_w_comp_[j][leaf] - m_comp_[j][leaf] * VecToSkew(r_po_w_comp_[j][leaf]) * VecToSkew(r_po_w_comp_[j][leaf]); 

      M_com_.block<1,1>(6 + j, 6 + j) = joint_j->S_w_.transpose() * M_single * joint_j->S_w_;   

      Eigen::MatrixXd r = Eigen::MatrixXd::Identity(6, 6);
      while (parent_link->parent_joint_idx_ != -1 & parent_link->parent_joint_idx_ != root_joint->idx_){
        Joint* joint_i = joints_[parent_link->parent_joint_idx_];
        parent_link = links_[joint_i->parent_link_idx_];
        child_link = links_[joint_i->child_link_idx_];

        r.block<3,3>(0, 3) = -VecToSkew(joint_j->r_wc_w_ - joint_i->r_wc_w_ );

        M_com_.block<1,1>(6 + j, 6 + joint_i->moving_idx_) = joint_j->S_w_.transpose() * M_single * r * joint_i->S_w_;
        M_com_.block<1,1>(6 + joint_i->moving_idx_, 6 + j) = M_com_.block<1,1>(6 + j, 6 + joint_i->moving_idx_);
      }    

      // calculate M_base_leg coupling      
      r.block<3,3>(0, 3) = -VecToSkew(joint_j->r_wc_w_ - root_joint->r_wc_w_);
      M_com_.block<1,6>(6 + j, 0) = joint_j->S_w_.transpose() * M_single * r;
      M_com_.block<6,1>(0, 6 + j) = M_com_.block<1,6>(6 + j, 0);
    }

    // calculate base mass matrix
    double m_base = root->mass_;
    Eigen::Vector3d r_po_w_base = root->r_po_w_; 
    Eigen::Matrix3d I_o_w_base = root->inertia_w_;

    for (auto child : root->children_joint_idx_){
      Joint* joint = joints_[child];
      int start = joint->moving_idx_;
      int leaf = joint->moving_idx_;
      while (1)
      {
        if (links_[moving_joints_[leaf]->child_link_idx_]->children_joint_idx_.size() == 0){
          break;
        }
        leaf = joints_[links_[moving_joints_[leaf]->child_link_idx_]->children_joint_idx_[0]]->moving_idx_; // TODO: consider multiple children
      }
      Eigen::Vector3d new_r_po_w_base = (m_base * r_po_w_base + m_comp_[start][leaf] * (moving_joints_[start]->r_pc_w_ + r_po_w_comp_[start][leaf])) / (m_base + m_comp_[start][leaf]);
      double new_m_base = m_base + m_comp_[start][leaf];

      Eigen::Vector3d r_1 = r_po_w_base - new_r_po_w_base;
      Eigen::Vector3d r_2 = (moving_joints_[start]->r_pc_w_ + r_po_w_comp_[start][leaf]) - new_r_po_w_base;

      I_o_w_base = I_o_w_base + I_o_w_comp_[start][leaf] - m_base * VecToSkew(r_1) * VecToSkew(r_1) - m_comp_[start][leaf] * VecToSkew(r_2) * VecToSkew(r_2);
      m_base = new_m_base;
      r_po_w_base = new_r_po_w_base;
    }

    Eigen::MatrixXd M_single = Eigen::MatrixXd::Identity(6, 6);
    M_single.block<3,3>(0, 0) = m_base * Eigen::MatrixXd::Identity(3, 3);
    M_single.block<3,3>(0, 3) = -m_base * VecToSkew(r_po_w_base);
    M_single.block<3,3>(3, 0) = m_base * VecToSkew(r_po_w_base);
    M_single.block<3,3>(3, 3) = I_o_w_base - m_base * VecToSkew(r_po_w_base) * VecToSkew(r_po_w_base); 
    M_com_.block<6,6>(0, 0) = M_single;
  }

  void RecursiveNonlinearities(Joint* joint) {
    Link* child_link = links_[joint->child_link_idx_];

    // calculate b_single
    child_link->b_single_.segment(0, 3) = child_link->mass_ * VecToSkew(joint->omega_cp_w_) * VecToSkew(joint->omega_cp_w_) * child_link->r_po_w_;
    child_link->b_single_.segment(3, 3) = VecToSkew(joint->omega_cp_w_) * (child_link->inertia_w_ - child_link->mass_ * VecToSkew(child_link->r_po_w_) * VecToSkew(child_link->r_po_w_)) * joint->omega_cp_w_;
    Eigen::Matrix<double, 6, 1> FT_by_gravity;
    FT_by_gravity << child_link->mass_ * gravity_, VecToSkew(child_link->r_po_w_) * (child_link->mass_ * gravity_);

    child_link->b_single_ = child_link->b_single_ - FT_by_gravity;

    // calculate M_single
    child_link->M_single_.block<3,3>(0, 0) = child_link->mass_ * Eigen::MatrixXd::Identity(3, 3);
    child_link->M_single_.block<3,3>(0, 3) = -child_link->mass_ * VecToSkew(child_link->r_po_w_);
    child_link->M_single_.block<3,3>(3, 0) = child_link->mass_ * VecToSkew(child_link->r_po_w_);
    child_link->M_single_.block<3,3>(3, 3) = child_link->inertia_w_ - child_link->mass_ * VecToSkew(child_link->r_po_w_) * VecToSkew(child_link->r_po_w_);

    // make twist
    Eigen::Matrix<double, 6, 1> twist_dot;
      twist_dot << joint->a_cp_w_ , joint->alpha_cp_w_;

    if (joint->moving_idx_ != -1)
      FT_[joint->moving_idx_] = child_link->M_single_ * twist_dot + child_link->b_single_;
    else
      FT_root_ = child_link->M_single_ * twist_dot + child_link->b_single_;

    // cal inertia and non-linear term by children
    for (auto child_joint_idx : child_link->children_joint_idx_){
      Joint* child_joint = joints_[child_joint_idx];
      RecursiveNonlinearities(child_joint);

      Eigen::Matrix<double, 6, 1> FT_by_child;
      FT_by_child << FT_[child_joint->moving_idx_].segment(0, 3), FT_[child_joint->moving_idx_].segment(3, 3) + VecToSkew(child_joint->r_pc_w_) * (FT_[child_joint->moving_idx_].segment(0, 3));
      
      if (joint->moving_idx_ != -1)
        FT_[joint->moving_idx_] = FT_[joint->moving_idx_] + FT_by_child;
      else
        FT_root_ = FT_root_ + FT_by_child;
    }

    if (joint->moving_idx_ != -1)
      b_.segment(6 + joint->moving_idx_, 1) = joint->S_w_.transpose() * FT_[joint->moving_idx_];
    else
      b_.segment(0, 6) = FT_root_;
  }

  void CalNonlinearities() {
    b_ = Eigen::VectorXd::Zero(6 + num_of_moving_joints_);
    FT_root_ = Eigen::VectorXd::Zero(6);
    for (int i = 0; i < num_of_moving_joints_; i++){
      FT_.push_back(Eigen::VectorXd::Zero(6));
    }

    RecursiveNonlinearities(root_joint);
  }
  
  void RecursiveArticulatedBodyInertia(Joint* joint) {
    Link* child_link = links_[joint->child_link_idx_];

    // calculate b_single
    child_link->b_single_.segment(0, 3) = child_link->mass_ * VecToSkew(joint->omega_cp_w_) * VecToSkew(joint->omega_cp_w_) * child_link->r_po_w_;
    child_link->b_single_.segment(3, 3) = VecToSkew(joint->omega_cp_w_) * (child_link->inertia_w_ - child_link->mass_ * VecToSkew(child_link->r_po_w_) * VecToSkew(child_link->r_po_w_)) * joint->omega_cp_w_;
    Eigen::Matrix<double, 6, 1> FT_by_gravity;
    FT_by_gravity << child_link->mass_ * gravity_, VecToSkew(child_link->r_po_w_) * (child_link->mass_ * gravity_);

    child_link->b_single_ = child_link->b_single_ - FT_by_gravity;

    // calculate M_single
    child_link->M_single_.block<3,3>(0, 0) = child_link->mass_ * Eigen::MatrixXd::Identity(3, 3);
    child_link->M_single_.block<3,3>(0, 3) = -child_link->mass_ * VecToSkew(child_link->r_po_w_);
    child_link->M_single_.block<3,3>(3, 0) = child_link->mass_ * VecToSkew(child_link->r_po_w_);
    child_link->M_single_.block<3,3>(3, 3) = child_link->inertia_w_ - child_link->mass_ * VecToSkew(child_link->r_po_w_) * VecToSkew(child_link->r_po_w_);
    // set M_a and b_a about current link
    if (joint->moving_idx_ != -1)
    {
      M_a_[joint->moving_idx_] = child_link->M_single_;
      b_a_[joint->moving_idx_] = child_link->b_single_;
    }
    else
    {
      M_root_a_ = child_link->M_single_;
      b_root_a_ = child_link->b_single_;      
    }

    // add inertia and non-linear term by children
    for (auto child_joint_idx : child_link->children_joint_idx_){
      Joint* child_joint = joints_[child_joint_idx];
      RecursiveArticulatedBodyInertia(child_joint);

      // make twist
      Eigen::Matrix<double, 6, 1> twist;
      twist << joint->v_cp_w_ , joint->omega_cp_w_;
      Eigen::MatrixXd X = Eigen::MatrixXd::Identity(6, 6);
      Eigen::MatrixXd X_dot = Eigen::MatrixXd::Zero(6, 6);
      Eigen::VectorXd S = Eigen::VectorXd::Zero(6);
      Eigen::VectorXd S_dot = Eigen::VectorXd::Zero(6);
      X.block<3,3>(3, 0) = VecToSkew(child_joint->r_pc_w_);
      X_dot.block<3,3>(3, 0) = VecToSkew(child_joint->v_cp_w_ - joint->v_cp_w_);
      S = child_joint->S_w_;
      S_dot = child_joint->S_dot_w_;

      if (joint->moving_idx_ != -1)
      {
        M_a_[joint->moving_idx_] = M_a_[joint->moving_idx_] + X * M_a_[child_joint->moving_idx_] * (-S / (S.transpose() * M_a_[child_joint->moving_idx_] * S)* 
                  (S.transpose() * M_a_[child_joint->moving_idx_] * X.transpose()) + X.transpose());
        b_a_[joint->moving_idx_] = b_a_[joint->moving_idx_] + X * (M_a_[child_joint->moving_idx_] * (S / (S.transpose() * M_a_[child_joint->moving_idx_] * S) * 
                  (gf_[child_joint->moving_idx_] - S.transpose() * M_a_[child_joint->moving_idx_] * (S_dot * gv_[child_joint->moving_idx_] + X_dot.transpose() * twist) 
                  - S.transpose() * b_a_[child_joint->moving_idx_]) + S_dot * gv_[child_joint->moving_idx_] + X_dot.transpose() * twist) + b_a_[child_joint->moving_idx_]);
      }
      else
      {
        M_root_a_ = M_root_a_ + X * M_a_[child_joint->moving_idx_] * (-S / (S.transpose() * M_a_[child_joint->moving_idx_] * S)* 
                  (S.transpose() * M_a_[child_joint->moving_idx_] * X.transpose()) + X.transpose());
        b_root_a_ = b_root_a_ + X * (M_a_[child_joint->moving_idx_] * (S / (S.transpose() * M_a_[child_joint->moving_idx_] * S) * 
                  (gf_[child_joint->moving_idx_] - S.transpose() * M_a_[child_joint->moving_idx_] * (S_dot * gv_[child_joint->moving_idx_] + X_dot.transpose() * twist) 
                  - S.transpose() * b_a_[child_joint->moving_idx_]) + S_dot * gv_[child_joint->moving_idx_] + X_dot.transpose() * twist) + b_a_[child_joint->moving_idx_]);
      }

    }
  }

  void CalArticulatedBodyInertia() {    
    RecursiveArticulatedBodyInertia(root_joint);

    // for (int i = 0; i < num_of_moving_joints_; i++){
    //   cout << "joint " << moving_joints_[i]->name_ << endl;
    //   cout << M_a_[i] << endl;
    //   cout << b_a_[i] << endl;
    // }
    // cout << "root"
    // cout << M_root_a_ << endl;
    // cout << b_root_a_ << endl;
  }

  void RecursiveGeneralizedAcceleration(Joint* joint) {
    Link* child_link = links_[joint->child_link_idx_];
    Joint* parent_joint = joints_[links_[joint->parent_link_idx_]->parent_joint_idx_];

    Eigen::MatrixXd X = Eigen::MatrixXd::Identity(6, 6);
    Eigen::MatrixXd X_dot = Eigen::MatrixXd::Zero(6, 6);
    Eigen::VectorXd S = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd S_dot = Eigen::VectorXd::Zero(6);
    Eigen::Matrix<double, 6, 1> twist_p, twist_dot_p, twist_dot_b;
    
    X.block<3,3>(3, 0) = VecToSkew(joint->r_pc_w_);
    X_dot.block<3,3>(3, 0) = VecToSkew(joint->v_cp_w_ - parent_joint->v_cp_w_);
    twist_p << parent_joint->v_cp_w_, parent_joint->omega_cp_w_;
    twist_dot_p << parent_joint->a_cp_w_, parent_joint->alpha_cp_w_;
    S = joint->S_w_;
    S_dot = joint->S_dot_w_;

    // calculate generalized acceleration
    ga_[joint->moving_idx_] = 1 / (S.transpose() * M_a_[joint->moving_idx_] * S) * (gf_[joint->moving_idx_] - S.transpose() * M_a_[joint->moving_idx_] * (S_dot * gv_[joint->moving_idx_] + X_dot.transpose() * twist_p + X.transpose() * twist_dot_p) 
              - S.transpose() * b_a_[joint->moving_idx_]);

    // calculate twist_dot_b
    twist_dot_b = X.transpose() * twist_dot_p + X_dot.transpose() * twist_p + S * ga_[joint->moving_idx_] + S_dot * gv_[joint->moving_idx_];
    joint->a_cp_w_ = twist_dot_b.segment(0, 3);
    joint->a_c_w_ = joint->a_cp_w_ - (joint->S_dot_w_ * gv_[joint->moving_idx_]).segment(0, 3);
    joint->alpha_cp_w_ = twist_dot_b.segment(3, 3);
    joint->alpha_c_w_ = joint->alpha_cp_w_ - (joint->S_dot_w_ * gv_[joint->moving_idx_]).segment(3, 3);

    for (int i = 0; i < child_link->children_joint_idx_.size(); i++){
      Joint* child_joint = joints_[child_link->children_joint_idx_[i]];
      RecursiveGeneralizedAcceleration(child_joint);
    }
  }

  void CalGeneralizedAcceleration() {
    // cal body acceleration and angular acceleration
    root_joint->a_c_w_ = (M_root_a_.inverse() * (gf_root_w_ - b_root_a_)).segment(0, 3);
    root_joint->a_cp_w_ = root_joint->a_c_w_;
    root_joint->alpha_c_w_ = (M_root_a_.inverse() * (gf_root_w_ - b_root_a_)).segment(3, 3);
    root_joint->alpha_cp_w_ = root_joint->alpha_c_w_;

    // cal moving joint acceleration and angular acceleration
    ga_ = Eigen::VectorXd::Zero(num_of_moving_joints_);
    for (int i = 0; i < root->children_joint_idx_.size(); i++){
      Joint* child_joint = joints_[root->children_joint_idx_[i]];
      RecursiveGeneralizedAcceleration(child_joint);
    }
  }

  // get
  Eigen::MatrixXd GetMassMatrix(){
    CalMassMatrix();
    return M_com_;
  }

  Eigen::MatrixXd GetNonlinearities(){
    CalNonlinearities();
    return b_;
  } 

  Eigen::MatrixXd GetGeneralizedAcceration(){
    CalArticulatedBodyInertia();
    CalGeneralizedAcceleration();

    Eigen::MatrixXd ga = Eigen::MatrixXd::Zero(6 + num_of_moving_joints_, 1);
    ga.block<3, 1>(0, 0) = root_joint->a_cp_w_;
    ga.block<3, 1>(3, 0) = root_joint->alpha_cp_w_;
    ga.block(6, 0, num_of_moving_joints_, 1) = ga_;

    return ga;
  }

};

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  Eigen::VectorXd gv = Eigen::VectorXd::Zero(gc.size());
  Eigen::VectorXd gf = Eigen::VectorXd::Zero(gc.size());
  ArticulatedSystem anymal(gc, gv, gf);
  anymal.ForwardKinematics();

  return anymal.GetMassMatrix(); // Eigen::MatrixXd::Ones(18,18); // 
}

inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  Eigen::VectorXd gf = Eigen::VectorXd::Zero(gc.size());
  ArticulatedSystem anymal(gc, gv, gf);
  anymal.ForwardKinematics();

  return anymal.GetNonlinearities();
}

inline Eigen::VectorXd computeGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  ArticulatedSystem anymal(gc, gv, gf);
  anymal.ForwardKinematics();

  return anymal.GetGeneralizedAcceration();
}