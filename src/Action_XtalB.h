#ifndef INC_ACTION_XTALB_H
#define INC_ACTION_XTALB_H
#include "Action.h"
/// <Enter description of Action_XtalB here>
class Action_XtalB : public Action {
  public:
    Action_XtalB() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_XtalB(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

};
#endif
