#include "Action_XtalB.h"
#include "CpptrajStdio.h"

// Action_XtalB::Help()
void Action_XtalB::Help() const {

}

// Action_XtalB::Init()
Action::RetType Action_XtalB::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  printf("Hello world!  Intializing!\n");
  return Action::OK;
}

// Action_XtalB::Setup()
Action::RetType Action_XtalB::Setup(ActionSetup& setup)
{
  printf("Hello world!  Setting up!\n");
  return Action::OK;
}

// Action_XtalB::DoAction()
Action::RetType Action_XtalB::DoAction(int frameNum, ActionFrame& frm)
{
  printf("Hello world!  Working!\n");
  return Action::OK;
}
