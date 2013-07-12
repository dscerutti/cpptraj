// Action_FixAtomOrder 
#include "Action_FixAtomOrder.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Action_FixAtomOrder::Action_FixAtomOrder() : debug_(0), newParm_(0) {} 

Action_FixAtomOrder::~Action_FixAtomOrder() {
  if (newParm_!=0) delete newParm_;
}

void Action_FixAtomOrder::Help() {
  mprintf("\t[outprefix <name>]\n"
          "\tFix atom ordering so that all atoms in molecules are sequential.\n");
}

// Action_FixAtomOrder::init()
Action::RetType Action_FixAtomOrder::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  prefix_ = actionArgs.GetStringKey("outprefix");
  mprintf("    FIXATOMORDER: Will attempt to fix atom ordering when atom numbering\n"
          "                  in molecules is non-sequential.\n");
  if (!prefix_.empty())
    mprintf("\tRe-ordered topology will be output with prefix %s\n",prefix_.c_str());
  return Action::OK;
}

void Action_FixAtomOrder::VisitAtom(int atomnum, int mol, Topology const& Parm) {
  // Return if this atom already has a molecule number
  if (molNums_[atomnum]!=-1) return;
  // Mark this atom as visited
  molNums_[atomnum] = mol;
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = Parm[atomnum].bondbegin();
                           bondedatom != Parm[atomnum].bondend(); bondedatom++)
    VisitAtom(*bondedatom, mol, Parm);
}


// Action_FixAtomOrder::setup()
Action::RetType Action_FixAtomOrder::Setup(Topology* currentParm, Topology** parmAddress) {
  // If topology already has molecule info assume no need to reorder.
  if (currentParm->Nmol() > 0) {
    mprintf("Warning: %s already has molecule information. No reordering will occur.\n",
            currentParm->c_str());
    return Action::ERR;
  }
  molNums_.resize( currentParm->Natom(), -1 );
  // Perform recursive search along bonds of each atom.
  int Nmol = 0;
  for (int atomnum = 0; atomnum < currentParm->Natom(); ++atomnum)
  {
    if (molNums_[atomnum] == -1) {
      VisitAtom( atomnum, Nmol, *currentParm );
      ++Nmol;
    }
  }
  mprintf("\tDetected %i molecules.\n", Nmol);
  if (debug_ > 0) {
    for (MapType::const_iterator mnum = molNums_.begin(); mnum != molNums_.end(); ++mnum)
      mprintf("\t\tAtom %u assigned to molecule %i\n", mnum - molNums_.begin() + 1, *mnum + 1);
  }
  // Place all atoms in molecule 0 first, molecule 1 next and so on
  atomMap_.clear();
  atomMap_.reserve( currentParm->Natom() );
  for (int mol = 0; mol < Nmol; mol++) {
    for (int atomnum = 0; atomnum < currentParm->Natom(); ++atomnum)
    {
      if (molNums_[atomnum] == mol)
        atomMap_.push_back( atomnum );
    }
  }
  if (debug_ > 0) {
    mprintf("\tNew atom mapping:\n");
    for (MapType::const_iterator atom = atomMap_.begin();
                                 atom != atomMap_.end(); ++atom)
      mprintf("\t\tNew atom %8u => old atom %8i\n", atom - atomMap_.begin() + 1, *atom + 1);
  }
  // Create new topology based on map
  if (newParm_ != 0) delete newParm_;
  newParm_ = currentParm->ModifyByMap( atomMap_ );
  if (newParm_ == 0) {
    mprinterr("Error: Could not create re-ordered topology.\n");
    return Action::ERR;
  }
  newParm_->Summary();
  // Allocate space for new frame
  newFrame_.SetupFrameM( newParm_->Atoms() );

  // If prefix given then output stripped parm
  if (!prefix_.empty()) {
    std::string newfilename = prefix_ + "." + currentParm->OriginalFilename();
    mprintf("\tWriting out amber topology file %s to %s\n",newParm_->c_str(),newfilename.c_str());
    ParmFile pfile;
    if ( pfile.Write( *newParm_, newfilename, ParmFile::AMBERPARM, 0 ) ) {
      mprinterr("Error: Could not write out stripped parm file %s\n",
                newParm_->c_str());
    }
  }

  // Set parm
  *parmAddress = newParm_;

  return Action::OK;  
}

// Action_FixAtomOrder::action()
Action::RetType Action_FixAtomOrder::DoAction(int frameNum, Frame* currentFrame, 
                                              Frame** frameAddress)
{
  // Reorder atoms in the frame
  newFrame_.SetCoordinatesByMap( *currentFrame, atomMap_ );
  *frameAddress = &newFrame_;
  return Action::OK;
}

