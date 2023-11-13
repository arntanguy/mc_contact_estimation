#include <mc_control/MCController.h>
#include <mc_observers/ObserverMacros.h>
#include "../include/mc_contact_estimation/ContactEstimation.h"
#include <mc_state_observation/gui_helpers.h>
#include "eigen-quadprog/QuadProg.h"
#include "eigen-quadprog/eigen_quadprog_api.h"
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
}

void ContactEstimation::reset(const mc_control::MCController & ctl)
{
  const auto & c = config_;
  residualsInt_ = Eigen::VectorXd::Zero(ctl.robot().mb().nrDof() - 6);
  integralsInt_ = Eigen::VectorXd::Zero(ctl.robot().mb().nrDof() - 6);
}

bool ContactEstimation::run(const mc_control::MCController & ctl)
{
  bool ret = true;
  auto & robot = ctl.realRobot();
  // auto & robot = ctl.robot();

  rbd::MultiBodyConfig mbc = robot.mbc();
  const sva::PTransformd & X_0_fb = robot.posW();

  sva::ForceVecd contactsSum_ = sva::ForceVecd::Zero();


  const int n = robot.mb().nrDof();
  rbd::ForwardDynamics fd(robot.mb());
  rbd::InverseDynamics id(robot.mb());
  rbd::Coriolis coriolis(robot.mb());

  fd.computeH(robot.mb(),robot.mbc());
  mbc.force[0] = sva::ForceVecd::Zero();
  fd.computeC(robot.mb(),mbc);
  const auto C = fd.C();
  const auto H = fd.H();
  const auto I_0c = H.block(0,0,6,6);
  const auto I_0c_inv =  I_0c.inverse();
  const auto F = H.block(0,6,6,n - 6);
  FT_Iinv_ = F.transpose() * I_0c_inv;
  mat = Eigen::MatrixXd::Zero(n,n);
  mat.block(6,0,n-6,6) = - FT_Iinv_;
  mat.block(6,6,n-6,n-6) = Eigen::MatrixXd::Identity(n-6,n-6);
  mat.block(0,0,6,6) = Eigen::Matrix6d::Identity();

  computeKnownContacts(robot,H,C,ctl.solver().contacts(),mbc);

  for(auto & c : ctl.solver().contacts())
  {
    const int indx = c.r1Surface()->bodyIndex(robot);
    const auto f = ctl.solver().desiredContactForce(c);
    // mc_rtc::log::info("err contact f\n{}",mbc.force[indx] - c.r1Surface()->X_0_s(robot).inv().dualMul(f));
    // mbc.force[indx] = c.r1Surface()->X_0_s(robot).inv().dualMul(f);
    contactsSum_ += X_0_fb.dualMul(mbc.force[indx]);
  }

  id.inverseDynamics(robot.mb(),mbc);
  auto C_mat = coriolis.coriolis(robot.mb(),mbc);

  
  

  const Eigen::MatrixXd H_fb = (mat * H).block(6,6,n-6,n-6);
  const Eigen::VectorXd C_fb = (mat * C).segment(6,n-6);
  const Eigen::Vector6d p_0c = C.segment(0,6);
  
  sva::MotionVecd a_0(Eigen::Vector3d::Zero(), mbc.gravity);
  const sva::PTransformd & X_p_i = mbc.parentToSon[0];

  const auto v0 = (X_0_fb * mbc.bodyVelW[0]).vector();

  const Eigen::MatrixXd Hdot_wb = C_mat + C_mat.transpose();
  const Eigen::MatrixXd Hdot = Hdot_wb.block(6,6,n - 6,n - 6);
  const Eigen::MatrixXd Fdot = Hdot_wb.block(0,6,6,n - 6);
  const Eigen::MatrixXd I_0c_dot = Hdot_wb.block(0,0,6,6);

  const Eigen::VectorXd qdot = flatten(mbc.alpha);
  const Eigen::VectorXd qddot = flatten(mbc.alphaD);

  const Eigen::MatrixXd H_fb_dot = Hdot - Fdot.transpose() * I_0c_inv * F - FT_Iinv_ * Fdot - Fdot.transpose() * (-I_0c_inv * I_0c_dot * I_0c_inv) * F;

  const auto tau_tot =flatten(mbc.jointTorque);  


  Eigen::VectorXd tau_c = Eigen::VectorXd::Zero(n);
  measuredContactsSum_ = sva::ForceVecd::Zero();
  for(auto & frame : forceSensorsFrame_)
  {
    tau_c += measuredContactTorque(frame,robot);
    const sva::PTransformd X_0_f = robot.frame(frame).position();
    measuredContactsSum_ += (X_0_fb * X_0_f.inv()).dualMul(robot.frame(frame).wrench());
  }
  // auto tau_g = gravityTorque(robot.mb(),mbc);
  // measuredContactsSum_ += sva::ForceVecd(tau_g.segment(0,6));

  // tau_c = contactTorque(robot.mb(),mbc);
  // measuredContactsSum_ = contactsSum_;
  // fd.computeC(robot.mb(),mbc);
  // const Eigen::VectorXd tau_fb = H * qddot + fd.C();
  // mc_rtc::log::info(tau_fb - tau_tot);

  if((tau_tot + tau_c).segment(6,n-6).norm() < 1e7)
  {
    integralsInt_ += (residualsInt_ + (mat*(tau_tot + tau_c)).segment(6,n-6) + H_fb_dot * qdot.segment(6,n-6) - C_fb) * dt_;
    integralsExt_ += (residualsExt_ + measuredContactsSum_.vector() +  (Hdot_wb.block(0,0,6,n) * qdot )-  p_0c )*dt_ ;

    residualsInt_ = gainInt_ * (H_fb * qdot.segment(6,n-6) - integralsInt_);
    residualsExt_ = gainExt_ * (H.block(0,0,6,n) * qdot - integralsExt_);
  }

  return ret;
}

void ContactEstimation::update(mc_control::MCController & ctl)
{
  if(!ctl.datastore().has("ContactEstimation::contactWrench"))
  {
    ctl.datastore().make_call("ContactEstimation::contactWrench",
      [this](const mc_rbdyn::Robot & robot,const std::vector<std::string> & frame, const std::vector<sva::PTransformd> & offsets) -> std::vector<sva::ForceVecd>
      {
        return contactWrench(robot,frame,offsets);
      }
      );
    ctl.datastore().make_call("ContactEstimation::contactWrenchSum",
      [this](const mc_rbdyn::Robot & robot,const std::string & frame, const sva::PTransformd & offsets) -> sva::ForceVecd
      {
        return contactWrenchSum(robot,frame,offsets);
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
  for (int i = 0 ; i < mbc.force.size() ; i++)
  {
    mbc.force[i] = sva::ForceVecd::Zero();
  }

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
    auto Jlocal = jac.bodyJacobian(mb,mbc);
    const sva::PTransformd X_0_b = mbc.bodyPosW[i];
    jac.fullJacobian(mb,Jlocal,J);
    Eigen::MatrixXd JT = J.transpose();
    tau += JT * (X_0_b.dualMul(mbc.force[i])).vector();
  }

  return tau;

}

Eigen::VectorXd ContactEstimation::measuredContactTorque(const std::string & frame,const mc_rbdyn::Robot & robot)
{

  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6,robot.mb().nrDof());

  const int ndof = robot.mb().nrDof();
  const sva::PTransformd & X_0_f = robot.frame(frame).position();
  const std::string & body_name = robot.frame(frame).body();
  const sva::PTransformd X_0_b = robot.mbc().bodyPosW[robot.mb().bodyIndexByName(body_name)];
  
  const sva::ForceVecd force = robot.frame(frame).wrench();
  auto jac = rbd::Jacobian(robot.mb(),body_name);
  auto Jlocal = jac.bodyJacobian(robot.mb(),robot.mbc());
  jac.fullJacobian(robot.mb(),Jlocal,J);
  
  return J.transpose() * ((X_0_b * X_0_f.inv()).dualMul(force)).vector();


}

std::vector<sva::ForceVecd> ContactEstimation::contactWrench(const mc_rbdyn::Robot & robot,const std::vector<std::string> & frame, const std::vector<sva::PTransformd> & offsets)
{
  const int ndof = robot.mb().nrDof(); 
  const sva::PTransformd X_0_fb = robot.posW();

  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(ndof,6 * frame.size());
  size_t count = 0;
  for (auto & fr : frame )
  {
    const std::string & body_name = robot.frame(fr).body();
    const sva::PTransformd & X_f_p = offsets[count];
    const sva::PTransformd & X_0_f = robot.frame(fr).position();
    const sva::PTransformd & X_0_b = robot.mbc().bodyPosW[ robot.mb().bodyIndexByName(body_name)];
    const sva::PTransformd X_b_p = X_0_b * (X_f_p * X_0_f).inv();
    auto jac = rbd::Jacobian(robot.mb(),body_name,X_b_p.translation());
    auto Jlocal = jac.bodyJacobian(robot.mb(),robot.mbc());
    Eigen::MatrixXd Jfull = Eigen::MatrixXd::Zero(6,ndof);
    jac.fullJacobian(robot.mb(),Jlocal,Jfull);
    J.block(0,count * 6,ndof,6) = mat * Jfull.transpose();
    count +=1;
  }
  const Eigen::MatrixXd Jinv = J.completeOrthogonalDecomposition().pseudoInverse();

  Eigen::VectorXd r = Eigen::VectorXd(ndof);
  r << residualsExt_, residualsInt_ ;
  const Eigen::VectorXd F = Jinv * r;
  std::vector<sva::ForceVecd> output;
  for (size_t i = 0 ; i < frame.size() ; i++ )
  {
    const std::string & body_name = robot.frame(frame[i]).body();
    const sva::PTransformd X_0_b = robot.mbc().bodyPosW[ robot.mb().bodyIndexByName(body_name)];
    const sva::PTransformd & X_0_f = robot.frame(frame[i]).position();
    const sva::PTransformd & X_b_p = offsets[i];

    output.push_back( (X_0_f * (( X_b_p * X_0_b ).inv()) ).dualMul(sva::ForceVecd(F.segment(6*i,6))) );


  }
  return output;
}

sva::ForceVecd ContactEstimation::contactWrenchSum(const mc_rbdyn::Robot & robot,const std::string & frame, const sva::PTransformd & offset)
{
  const sva::PTransformd X_0_fb = robot.posW();

  const sva::PTransformd & X_0_f = robot.frame(frame).position();

  return  (offset * X_0_f * (( X_0_fb ).inv()) ).dualMul(sva::ForceVecd(residualsExt_)) ;
  
}

void ContactEstimation::computeKnownContacts(const mc_rbdyn::Robot & robot, const Eigen::MatrixXd & H , const Eigen::VectorXd & C, const std::vector<mc_rbdyn::Contact> & contacts ,rbd::MultiBodyConfig & mbc)
{
  const auto & mb  = robot.mb();
  const int ndof = mb.nrDof();
  const int nvar = contacts.size() * 6;
  const auto qddot = flatten(mbc.alphaD);
  Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(6,nvar);
  const Eigen::VectorXd beq = H.block(0,0,6,ndof) * qddot + C.segment(0,6);
  sva::PTransformd X_0_fb = robot.posW();
  Eigen::MatrixXd Aineq = Eigen::MatrixXd::Zero(0,nvar);
  Eigen::VectorXd bineq =Eigen::VectorXd::Zero(Aineq.rows());

  size_t count = 0;

  std::vector<double> body_indx;

  for(auto & c : contacts)
  {
   
    const sva::PTransformd X_0_s = c.r1Surface()->X_0_s(robot);

    Aeq.block(0,6 * count,6,6) = (X_0_fb * X_0_s.inv()).dualMatrix();

    count +=1;

  }

  Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(nvar,nvar);
  Eigen::VectorXd p = Eigen::VectorXd::Zero(nvar);

  // QP Problem
  Eigen::QuadProgDense QP;

  QP.problem(nvar, Aeq.rows(), Aineq.rows());
  bool QPsuccess = QP.solve(Q, p, Aeq, beq, Aineq, bineq);
  count = 0;
  for(auto & c : contacts)
  {
    const int indx = c.r1Surface()->bodyIndex(robot);
    const sva::PTransformd X_0_s = c.r1Surface()->X_0_s(robot);
    mbc.force[indx] = X_0_s.inv().dualMul(sva::ForceVecd(QP.result().segment(6*count , 6)));
    count+=1;
  }


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
                  mc_rtc::gui::NumberInput("Gain Int",[this]() -> const double {return gainInt_;}, [this](const double i){gainInt_ = i;}),
                  mc_rtc::gui::NumberInput("Gain Ext",[this]() -> const double {return gainExt_;}, [this](const double i){gainExt_ = i;}),
                  mc_rtc::gui::Label("Residuals Int norm",[this]() -> const std::string {return std::to_string(residualsInt_.norm());}),     
                  mc_rtc::gui::ArrayInput("Residuals Ext",{"cx","cy","cz","fx","fy","fz"},[this]() -> Eigen::Vector6d
                                                                            {
                                                                              return residualsExt_;
                                                                            },
                                                                            [this](const Eigen::Vector6d & t) {}
                                                                            ));
}


} // namespace mc_state_observation
EXPORT_OBSERVER_MODULE("ContactEstimation", mc_state_observation::ContactEstimation)
