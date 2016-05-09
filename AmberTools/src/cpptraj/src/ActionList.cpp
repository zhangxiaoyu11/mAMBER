// ActionList
#include "ActionList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ActionList::ActionList() : debug_(0) {}

// DESTRUCTOR
ActionList::~ActionList() {
  Clear();
}

void ActionList::Clear() {
  // No need to cast back to whatever action was allocd since Action destructor is virtual
  for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act)
    delete *act;
  actionlist_.clear();
  actioncmd_.clear();
  actionstatus_.clear(); 
}

// ActionList::SetDebug()
void ActionList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("ActionList DEBUG LEVEL SET TO %i\n",debug_);
}

int ActionList::AddAction(DispatchObject::DispatchAllocatorType Alloc, ArgList& argIn,
                          TopologyList* PFL, FrameList* FL, DataSetList* DSL,
                          DataFileList* DFL)
{
  Action* act = (Action*)Alloc();
  // Attempt to initialize action
  if ( act->Init( argIn, PFL, FL, DSL, DFL, debug_ ) != Action::OK ) {
    mprinterr("Error: Could not initialize action [%s]\n", argIn.Command());
    delete act;
    return 1;
  }
  argIn.CheckForMoreArgs();
  actionlist_.push_back( act );
  actioncmd_.push_back( argIn.ArgLine() );
  actionstatus_.push_back( INIT );
  return 0;
}

// ActionList::SetupActions()
/** Attempt to set up all actions in the action list with the given parm
  * If an action cannot be set up skip it.
  */
int ActionList::SetupActions(Topology **ParmAddress) {
  if (actionlist_.empty()) return 0;
  Topology *OriginalParm = *ParmAddress;
  mprintf(".....................................................\n");
  mprintf("ACTION SETUP FOR PARM [%s] (%zu actions):\n",(*ParmAddress)->c_str(),actionlist_.size());
  unsigned int actnum = 0;
  for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    // Only attempt to set up action if active 
    if (actionstatus_[actnum] != INACTIVE) {
      mprintf("  %u: [%s]\n", actnum, actioncmd_[actnum].c_str());
      actionstatus_[actnum] = SETUP;
      Action::RetType err = (*act)->Setup(*ParmAddress, ParmAddress);
      if (err == Action::ERR) {
        mprintf("Warning: Setup failed for [%s]: Skipping\n",
                actioncmd_[actnum].c_str());
        // Reset action status to INIT (pre-setup)
        actionstatus_[actnum] = INIT;
        //return 1;
      } else if (err == Action::USEORIGINALFRAME) {
        // Return value of 2 requests return to original parm
        *ParmAddress = OriginalParm;
      }
      //fprintf(stdout,"DEBUG: After Action %i Setup parmName is %s\n",act,P->parmName);
    }
    ++actnum;
  }
  //mprintf(".....................................................\n");

  return 0;
}

// ActionList::DoActions()
/** Perform actions in the action list on the given Frame. Skip actions not 
  * initialized or not setup. 
  * \param FrameAddress Memory address of the current frame.
  * \param frameNumIn The current frame number.
  * \return true if coordinate output should be suppressed.
  * \return false if coordinate output should be performed.
  */
bool ActionList::DoActions(Frame **FrameAddress, int frameNumIn) {
  Frame *OriginalFrame = *FrameAddress;

  //fprintf(stdout,"DEBUG: Performing %i actions on frame %i.\n",Naction,frameNumIn);
  unsigned int actnum = 0;
  for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act) 
  {
    // Only do actions which were properly set up
    if (actionstatus_[actnum] == SETUP) { 
      // Perform action on frame
      Action::RetType err = (*act)->DoAction(frameNumIn, *FrameAddress, FrameAddress);
      // Check for action special conditions/errors
      if (err != Action::OK) {
        if (err == Action::USEORIGINALFRAME) {
          // Return value of 2 requests return to original frame
          *FrameAddress = OriginalFrame;
        } else if (err == Action::SUPPRESSCOORDOUTPUT) {
          // Skip the rest of the actions and suppress output. Necessary when
          // e.g. performing a running average over coords.
          return true;
        } else {
          // If here return type is ACTION_ERR.
          // Treat actions that fail as if they could not be set up
          mprintf("Warning: Action [%s] failed, frame %i.\n", actioncmd_[actnum].c_str(),
                frameNumIn);
          actionstatus_[actnum] = INIT;
        }
      }
    }
    ++actnum;
  }
  return false;
}

// ActionList::Print()
void ActionList::Print() {
  unsigned int actnum = 0;
  for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    // Skip deactivated actions
    if (actionstatus_[actnum++] != INACTIVE)
      (*act)->Print();
  }
}

void ActionList::List() const {
  mprintf("ACTIONS:\n");
  if (actionlist_.empty())
    mprintf("  No Actions.\n");
  else
    for (unsigned int actnum = 0; actnum < actionlist_.size(); ++actnum)
      mprintf("  %u: [%s]\n", actnum, actioncmd_[actnum].c_str());   
}
