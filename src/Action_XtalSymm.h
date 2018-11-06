#ifndef INC_ACTION_XTALSYMM_H
#define INC_ACTION_XTALSYMM_H
#include "Action.h"
#include "Matrix_3x3.h"
#include "ReferenceAction.h"
#include "Vec3.h"

#define IASU_GRID_BINS   192
#define DASU_GRID_BINS   192.0

//---------------------------------------------------------------------------------------------
// XtalDock: stores the results of crystal symmetry operations on solitary subunits,
//           attempting to align them back to the original asymmetric unit.  The goal
//           is to cluster the origins of multiple operations such that one cluster
//           has all of its required origins in about the same spot and gives
//           consistently low atom positional RMSDs.
//---------------------------------------------------------------------------------------------
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

//---------------------------------------------------------------------------------------------
// TransOp: stores an initial translation (as three doubles) plus the index of an operation.
//          This information will help guide points (or entire coordinate sets) from some
//          starting position into one of the regions of space for a given asymmetric unit.
//---------------------------------------------------------------------------------------------
class TransOp {
  public:
    int opID;      // Operation index
    double tr_x;   // Initial X translation
    double tr_y;   // Initial Y translation
    double tr_z;   // Initial Z translation
};

//---------------------------------------------------------------------------------------------
// XtalSymm: an action to superimpose symmetry-related parts of a simulation system using
//           crystallographic symmetry operations
//---------------------------------------------------------------------------------------------
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
    int sgID;              // The space group identification number
 
    // Reference frame
    AtomMask tgtMask_;
    ReferenceAction REF_;  
    Frame RefFrame_;       // The reference (or first trajectory) frame
    bool useFirst_;        // Flag to use first trajectory frame, true if no reference frame

    // Re-imaging for solvent atoms
    int nMolecule;         // Total number of molecules in the system, taken from the topology
    bool allToFirstASU_;   // Flag to have all atoms not explicitly in a designated asymmetric
                           //   unit re-imaged to the primary ASU volume
    bool molCentToASU_;    // Flag to use molecule centroids, not individual atoms, in the
                           //   above re-imaging
    int* molLimits;        // Start and end points for each molecule (all molecules are assumed
                           //   to be contiguous within the topology, but it is not assumed
                           //   that molecule i+1 starts where molecule i ends)
    bool* molInSolvent;    // Flags to indicate whether each molecule is part of the non-ASU,
                           //   free-floating "solvent" component
  
    // Masks for the asymmetric units and solvent particles
    int nmasks;
    AtomMask*  Masks;
    AtomMask   SolventMask;
    std::vector<int> subunitOpID;
  
    // Rotation matrices and translation vectors
    Matrix_3x3* R;
    Vec3* T;
    Vec3* RefT;
    bool* rotIdentity;

    // Grid for ASU assignment of loose molecules
    TransOp* AsuGrid;
  
    // Methods
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    Action::RetType LoadSpaceGroupSymOps();
    bool OperationAvailable(XtalDock* leads, int* HowToGetThere, int ncurr);
    bool OriginsAlign(XtalDock* leads, int* HowToGetThere, int ncurr);
    double BestSuperposition(int, int, XtalDock*, int&);
    Vec3 BestOrigin(Frame&, Frame*, std::vector<int>&);
    TransOp DetectAsuResidence(double x, double y, double z);
    Action::RetType BuildAsuGrid();
    void Print() {}
    void ClearMemory();
    double dmin(double, double);
    double dmin(double, double, double);
    double dmin(double, double, double, double);
    double dmax(double, double);
    double dmax(double, double, double);
    double dmax(double, double, double, double);
    bool PointInPrimaryASU(double x, double y, double z);
};
#endif
