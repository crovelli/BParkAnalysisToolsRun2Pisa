/// This BPark class is an auxiliary class which contains basic
/// functionality useful for any analysis of B-Parked events.
/// It derives from BParkBase.

#ifndef BPark_h
#define BPark_h

#include "./BParkBase.h"

#include <TLorentzVector.h>
#include <string>
#include <vector>
#include <map>

class BPark : public BParkBase{

public:

  BPark(TTree *tree=0);
  virtual ~BPark();

private:

protected:

};

#endif
