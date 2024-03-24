//
// Created by Jemin Hwangbo on 2022/03/17.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "exercise1_20234191.hpp"
#include "raisim/RaisimServer.hpp"
#include <cmath>

int main(int argc, char* argv[]) {
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();

  // anymal
  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");
  anymal->setName("anymal");
  server.focusOn(anymal);

  // anymal configuration
  Eigen::VectorXd jointNominalConfig(anymal->getGeneralizedCoordinateDim());
  jointNominalConfig << 0, 1, 0.54, 1.0, 0.0, 0.2, 0.3, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.8, -0.1, 0.2, -0.03, -0.4, 0.8;
  anymal->setGeneralizedCoordinate(jointNominalConfig);
  anymal->updateKinematics();

  // debug sphere
  auto debugSphere = server.addVisualSphere("debug_sphere", 0.02);
  debugSphere->setColor(1,0,0,1);

  Eigen::Vector3d position_e = getEndEffectorPosition(jointNominalConfig);
  debugSphere->setPosition(position_e);

  // solution sphere
  auto answerSphere = server.addVisualSphere("answer_sphere", 0.04);
  answerSphere->setColor(0,1,0,1);
  raisim::Vec<3> pos;
  anymal->getFramePosition("LH_shank_fixed_LH_FOOT", pos);
  answerSphere->setPosition(pos.e());

  cout << "mse : " << std::sqrt((position_e[0] - pos.e()[0])*(position_e[0] - pos.e()[0]) + (position_e[0] - pos.e()[0])*(position_e[0] - pos.e()[0]) + (position_e[0] - pos.e()[0])*(position_e[0] - pos.e()[0])) << endl;

  // visualization
  server.launchServer();
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}
