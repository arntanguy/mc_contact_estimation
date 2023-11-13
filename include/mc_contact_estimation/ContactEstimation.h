#pragma once

#include <mc_observers/Observer.h>
#include <mc_rbdyn/Robot.h>

namespace mc_state_observation
{

struct ContactEstimation : public mc_observers::Observer
{

 public:
  ContactEstimation(const std::string & type, double dt);

  void configure(const mc_control::MCController & ctl, const mc_rtc::Configuration &) override;

  void reset(const mc_control::MCController & ctl) override;

  bool run(const mc_control::MCController & ctl) override;

  void update(mc_control::MCController & ctl) override;


protected:

  void computeKnownContacts(const mc_rbdyn::Robot & robot, const Eigen::MatrixXd & H , const Eigen::VectorXd & C, const std::vector<mc_rbdyn::Contact> & contacts ,rbd::MultiBodyConfig & mbc);

  Eigen::VectorXd flatten(const std::vector<std::vector<double>> & vec);

  Eigen::VectorXd gravityTorque(const rbd::MultiBody & mb, rbd::MultiBodyConfig mbc);

  Eigen::VectorXd contactTorque(const rbd::MultiBody & mb, const rbd::MultiBodyConfig & mbc);

  Eigen::VectorXd measuredContactTorque(const std::string & frame,const mc_rbdyn::Robot & robot);

  std::vector<sva::ForceVecd> contactWrench(const mc_rbdyn::Robot & robot, const std::vector<std::string> & frame, const std::vector<sva::PTransformd> & offsets);

  sva::ForceVecd contactWrenchSum(const mc_rbdyn::Robot & robot,const std::string & frame, const sva::PTransformd & offsets);

  std::vector<std::vector<double>> unFlatten(const Eigen::VectorXd & vec , const rbd::MultiBody & mb);

  /*! \brief Add observer from logger
   *
   * @param category Category in which to log this observer
   */
  void addToLogger(const mc_control::MCController &, mc_rtc::Logger &, const std::string & category) override;

  /*! \brief Remove observer from logger
   *
   * @param category Category in which this observer entries are logged
   */
  void removeFromLogger(mc_rtc::Logger &, const std::string & category) override;

  /*! \brief Add observer information the GUI.
   *
   * @param category Category in which to add this observer
   */
  void addToGUI(const mc_control::MCController &,
                mc_rtc::gui::StateBuilder &,
                const std::vector<std::string> & /* category */) override;

protected:
  /// @{
  std::string robot_ = ""; ///< Name of robot to which the IMU sensor belongs
  
  std::string datastoreName_ = ""; ///< Name on the datastore (default name())

  mc_rtc::Configuration config_;

  Eigen::MatrixXd FT_Iinv_;

  double gainInt_ = 1.;
  double gainExt_ = 1.;

  Eigen::VectorXd integralsInt_;
  Eigen::VectorXd residualsInt_ ;

  Eigen::MatrixXd mat;

  sva::ForceVecd measuredContactsSum_= sva::ForceVecd::Zero();

  Eigen::Vector6d integralsExt_ = Eigen::Vector6d::Zero();
  Eigen::Vector6d residualsExt_ = Eigen::Vector6d::Zero();

  std::vector<std::string> forceSensorsFrame_;

  double dt_ = 5e-3;


};
} // namespace mc_state_observation


