#include <mc_control/MCController.h>
#include <mc_observers/ObserverMacros.h>
#include "../include/mc_contact_estimation/ContactEstimation.h"
#include <mc_state_observation/gui_helpers.h>
#include <RBDyn/ID.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/FA.h>
#include <RBDyn/Coriolis.h>
#include <RBDyn/Jacobian.h>
#include <Eigen/QR> 

namespace mc_state_observation
{

// namespace so = stateObservation;

ContactEstimation::ContactEstimation(const std::string & type, double dt): mc_observers::Observer(type, dt)
{

}

void ContactEstimation::configure(const mc_control::MCController & ctl, const mc_rtc::Configuration & config)
{
  config_.load(config);
  forceSensorsFrame_ = config_("force_sensors");
  config_("gainInt",gainInt_);
  config_("gainExt",gainExt_);
  dt_ = ctl.timeStep;
  residualsInt_ = Eigen::VectorXd::Zero(ctl.robot().mb().nrDof() - 6);
  residualsExt_ = Eigen::VectorXd::Zero(6);
}

void ContactEstimation::reset(const mc_control::MCController & ctl)
{
  const auto & c = config_;
  residualsInt_ = Eigen::VectorXd::Zero(ctl.robot().mb().nrDof() - 6);
  residualsExt_ = Eigen::VectorXd::Zero(6);
  integralsInt_ = Eigen::VectorXd::Zero(ctl.robot().mb().nrDof() - 6);
}

bool ContactEstimation::run(const mc_control::MCController & ctl)
{
  bool ret = true;
  auto & robot = ctl.realRobot();
  rbd::MultiBodyConfig mbc = robot.mbc();
  
  sva::ForceVecd contactsSum_ = sva::ForceVecd::Zero();

  for(auto & c : ctl.solver().contacts())
  {
    const auto f = ctl.solver().desiredContactForce(c);
    const int indx = c.r1Surface()->bodyIndex(robot);
    mbc.force[indx] = c.r1Surface()->X_0_s(robot).inv().dualMul(f);
    contactsSum_ += mbc.force[indx];
  }

  const int n = robot.mb().nrDof();
  rbd::ForwardDynamics fd(robot.mb());
  rbd::InverseDynamics id(robot.mb());
  rbd::Coriolis coriolis(robot.mb());

  id.inverseDynamics(robot.mb(),mbc);
  fd.computeH(robot.mb(),robot.mbc());
  fd.computeC(robot.mb(),robot.mbc());
  auto C_mat = coriolis.coriolis(robot.mb(),robot.mbc());

  const auto & I_0c = fd.H().block(0,0,6,6);
  const auto & p_0c = fd.C().segment(0,6);
  const auto I_0c_inv =  I_0c.inverse();
  const auto & F = fd.H().block(0,6,6,n - 6);
  FT_Iinv_ = F.transpose() * I_0c_inv;
  const auto & H_fb = fd.H().block(6,6,n - 6,n - 6) - FT_Iinv_ * F;
  const auto & C_fb = fd.C().segment(6,n - 6) - FT_Iinv_ * p_0c;
  const auto v0 = mbc.bodyVelW[0].vector();

  const Eigen::MatrixXd Hdot_wb = C_mat + C_mat.transpose();
  const Eigen::MatrixXd Hdot = Hdot_wb.block(6,6,n - 6,n - 6);
  const Eigen::MatrixXd Fdot = Hdot_wb.block(0,6,6,n - 6);
  const Eigen::MatrixXd I_0c_dot = Hdot_wb.block(0,0,6,6);

  const Eigen::VectorXd qdot = flatten(mbc.alpha);

  const Eigen::MatrixXd H_fb_dot = Hdot - Fdot.transpose() * I_0c_inv * F - FT_Iinv_ * Fdot - Fdot.transpose() * (-I_0c_inv * I_0c_dot * I_0c_inv) * F;

  auto tau_tot = flatten(mbc.jointTorque);  

  Eigen::VectorXd tau_c = Eigen::VectorXd::Zero(n);
  measuredContactsSum_ = sva::ForceVecd::Zero();
  for(auto & frame : forceSensorsFrame_)
  {
    tau_c += measuredContactTorque(frame,robot);
    measuredContactsSum_ += robot.frame(frame).position().inv().dualMul(robot.frame(frame).wrench());
  }
  // tau_c = contactTorque(robot.mb(),mbc);
  // measuredContactsSum_ = contactsSum_;

  if((tau_tot + tau_c).segment(6,n-6).norm() < 1e7)
  {
    integralsInt_ += (residualsInt_ + (tau_tot + tau_c).segment(6,n-6) + H_fb_dot * qdot.segment(6,n-6) - C_fb) * dt_;
    integralsExt_ += (residualsExt_ + I_0c_dot * v0 + Fdot * qdot.segment(6,n-6) - p_0c + measuredContactsSum_.vector())*dt_ ;

    residualsInt_ = gainInt_ * (H_fb * qdot.segment(6,n-6) - integralsInt_);
    residualsExt_ = gainExt_ * (I_0c * v0 + F * qdot.segment(6,n-6) - integralsExt_);
  }


  return ret;
}

void ContactEstimation::update(mc_control::MCController & ctl)
{
  if(!ctl.datastore().has("ContactEstimation::ContactWrench"))
  {
    ctl.datastore().make_call("ContactEstimation::ContactWrench",
      [this](const mc_rbdyn::Robot & robot,const std::vector<std::string> & frame, const std::vector<sva::PTransformd> & offsets) -> std::vector<sva::ForceVecd>
      {
        return contactWrench(robot,frame,offsets);
      }
      );
  }
}

Eigen::VectorXd ContactEstimation::flatten(const std::vector<std::vector<double>> & vec)
{
  std::vector<double> flatten_vec;
  size_t indx = 0;
  for (size_t i = 0 ; i < vec.size() ; i++)
  {
    for (size_t j = 0 ; j < vec[i].size(); j++)
    {
      flatten_vec.push_back(vec[i][j]);
    }
  }
  return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(flatten_vec.data(), flatten_vec.size());
}

std::vector<std::vector<double>> ContactEstimation::unFlatten(const Eigen::VectorXd & vec , const rbd::MultiBody & mb)
{
  std::vector<std::vector<double>> output;
  size_t count = 0;
  for (size_t i = 0 ; i < mb.nrBodies() ; i++ )
  {
    std::vector<double> joint;
    for (size_t j = 0 ; j < mb.joint(i).dof() ; j++ )
    {
      joint.push_back(vec(count));
      count +=1;
    }
    output.push_back(joint);
  }
  return output;

}

Eigen::VectorXd ContactEstimation::gravityTorque(const rbd::MultiBody & mb, rbd::MultiBodyConfig mbc)
{
  rbd::InverseDynamics id(mb);

  mbc.alpha = unFlatten(Eigen::VectorXd::Zero(mb.nrDof()),mb);
  mbc.alphaD = unFlatten(Eigen::VectorXd::Zero(mb.nrDof()),mb);

  rbd::forwardVelocity(mb,mbc);
  rbd::forwardAcceleration(mb,mbc);

  id.inverseDynamics(mb,mbc);
  return flatten(mbc.jointTorque);

}

Eigen::VectorXd ContactEstimation::contactTorque(const rbd::MultiBody & mb, const rbd::MultiBodyConfig & mbc)
{

  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6,mb.nrDof());
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(mb.nrDof());
  const int ndof = mb.nrDof();
  for (size_t i = 0 ; i < mb.nrBodies() ; i++)
  {
    auto jac = rbd::Jacobian(mb,mb.bodies()[i].name());
    auto Jlocal = jac.jacobian(mb,mbc);
    const sva::PTransformd X_0_b = mbc.bodyPosW[i];
    jac.fullJacobian(mb,Jlocal,J);
    Eigen::MatrixXd JT = J.transpose();
    JT.block(6,0,ndof - 6,6) -= FT_Iinv_;
    tau += JT * (sva::PTransformd(Eigen::Matrix3d::Identity(),X_0_b.translation()).dualMul(mbc.force[i])).vector();
  }

  return tau;

}

Eigen::VectorXd ContactEstimation::measuredContactTorque(const std::string & frame,const mc_rbdyn::Robot & robot)
{

  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6,robot.mb().nrDof());
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(robot.mb().nrDof());
  const int ndof = robot.mb().nrDof();
  const Eigen::Matrix3d R_frame_0 = robot.frame(frame).position().rotation().transpose();
  const std::string & body_name = robot.frame(frame).body();
  
  const sva::ForceVecd force = sva::PTransformd(R_frame_0,Eigen::Vector3d::Zero()).dualMul(robot.frame(frame).wrench());
  auto jac = rbd::Jacobian(robot.mb(),body_name);
  auto Jlocal = jac.jacobian(robot.mb(),robot.mbc());
  jac.fullJacobian(robot.mb(),Jlocal,J);
  Eigen::MatrixXd JT = J.transpose();
  JT.block(6,0,ndof - 6,6) -= FT_Iinv_;
  tau += JT * force.vector();

  return tau;

}

std::vector<sva::ForceVecd> ContactEstimation::contactWrench(const mc_rbdyn::Robot & robot,const std::vector<std::string> & frame, const std::vector<sva::PTransformd> & offsets)
{
  const int ndof = robot.mb().nrDof(); 
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(ndof,6 * frame.size());
  size_t count = 0;
  for (auto & fr : frame )
  {
    const std::string body_name = robot.frame(fr).body();
    const sva::PTransformd & X_frame_p = offsets[count];
    auto jac = rbd::Jacobian(robot.mb(),body_name,X_frame_p.translation());
    auto Jlocal = jac.jacobian(robot.mb(),robot.mbc());
    Eigen::MatrixXd Jfull = Eigen::MatrixXd::Zero(6,ndof);
    jac.fullJacobian(robot.mb(),Jlocal,Jfull);
    Eigen::MatrixXd JT = Jfull.transpose();
    JT.block(0,0,6,6) = Eigen::Matrix6d::Identity();
    JT.block(6,0,ndof - 6,6) -= FT_Iinv_;
    J.block(0,count * 6,ndof,6) = JT;
    count +=1;
  }
  const Eigen::MatrixXd Jinv = J.completeOrthogonalDecomposition().pseudoInverse();
  Eigen::VectorXd r = Eigen::VectorXd(ndof);
  r << residualsExt_, residualsInt_;
  const Eigen::VectorXd F = Jinv * r;
  std::vector<sva::ForceVecd> output;
  for (size_t i = 0 ; i < frame.size() ; i++ )
  {
    const sva::PTransformd & X_0_frame = robot.frame(frame[i]).position();
    output.push_back( sva::PTransformd( (offsets[i] * X_0_frame).rotation() , Eigen::Vector3d::Zero() ).dualMul( sva::ForceVecd(F.segment(6*i,6)) ));
  }
  return output;
}

void ContactEstimation::addToLogger(const mc_control::MCController &,
                                   mc_rtc::Logger & logger,
                                   const std::string & category)
{
  
}

void ContactEstimation::removeFromLogger(mc_rtc::Logger & logger, const std::string & category)
{
 
}

void ContactEstimation::addToGUI(const mc_control::MCController & ctl,
                                mc_rtc::gui::StateBuilder & gui,
                                const std::vector<std::string> & category)
{
  gui.addElement(category,
                  mc_rtc::gui::NumberInput("Gain",[this]() -> const double {return gainInt_;}, [this](const double i){gainInt_ = i;}),
                  mc_rtc::gui::NumberInput("Gain Ext",[this]() -> const double {return gainExt_;}, [this](const double i){gainExt_ = i;}),
                  mc_rtc::gui::Label("Residuals Int norm",[this]() -> const std::string {return std::to_string(residualsInt_.norm());}),     
                  mc_rtc::gui::Label("Residuals Ext norm",[this]() -> const std::string {return std::to_string(residualsExt_.norm());})      
                  );
}


} // namespace mc_state_observation
EXPORT_OBSERVER_MODULE("ContactEstimation", mc_state_observation::ContactEstimation)
