#ifndef INC_ACTION_XTALSYMM_H
#define INC_ACTION_XTALSYMM_H
#include "Action.h"
#include "Matrix_3x3.h"
#include "ReferenceAction.h"
#include "Vec3.h"

// XtalDock: stores the results of crystal symmetry operations on solitary subunits,
//           attempting to align them back to the original asymmetric unit.  The goal
//           is to cluster the origins of multiple operations such that one cluster
//           has all of its required origins in about the same spot and gives
//           consistently low atom positional RMSDs.
class XtalDock {
  public:
    int subunit;  // The subunit that this is operating on
    int opID;     // The operation in use, indexing the lists R and T in Action_XtalSymm
    double rmsd;  // The un-fitted rmsd between the original and superimposed subunits
    Vec3 displc;  // The optimal displacement beteen the two subunits' centers of mass,
                  //   scaled to simulation cell fractional coordinates
    Vec3 origin;  // The origin that got the best rmsd

    /// Assignment
    XtalDock& operator=(const XtalDock& rhs) {
      if (this == &rhs) return *this;
      subunit = rhs.subunit;
      opID    = rhs.opID;
      rmsd    = rhs.rmsd;
      displc  = rhs.displc;
      origin  = rhs.origin;
      return *this;
    }
};

// XtalSymm: an action to superimpose symmetry-related parts of a simulation system using
//           crystallographic symmetry operations
class Action_XtalSymm : public Action {
  public:
    Action_XtalSymm() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_XtalSymm(); }
    void Help() const;
    ~Action_XtalSymm();
  private:

    // The space group
    int nops;              // The number of symmetry operations in the canonical space group
    int nCopyA;            // The number of unit cell replicas along each of three dimensions.
    int nCopyB;            //   This generalizes the space group into supercells based on that
    int nCopyC;            //   space group, and allows us to talk about symmetry operations
                           //   in terms of the supercell.
    std::string spaceGrp;  //

    // Reference frame
    AtomMask tgtMask_;
    ReferenceAction REF_;
    Frame RefFrame_;
    bool useFirst_;
  
    // Masks for the asymmetric units
    int nmasks;
    AtomMask*  Masks;
    std::vector<int> subunitOpID;
  
    // Rotation matrices and translation vectors
    Matrix_3x3* R;
    Vec3* T;
    Vec3* RefT;
    bool* rotIdentity;

    // Methods
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    Action::RetType LoadSpaceGroupSymOps();
    bool OperationAvailable(XtalDock* leads, int* HowToGetThere, int ncurr);
    bool OriginsAlign(XtalDock* leads, int* HowToGetThere, int ncurr);
    double BestSuperposition(int, int, XtalDock*, int&);
    Vec3 BestOrigin(Frame&, Frame*, std::vector<int>&);
    void Print() {}
    void ClearMemory();
};
#endif
