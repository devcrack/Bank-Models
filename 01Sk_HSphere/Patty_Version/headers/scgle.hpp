#ifndef SCGLE_HPP
#define SCGLE_HPP

#include"file_manage.hpp"
#include<cstdlib>
#include <sys/stat.h>

class scgle {
private:
  double _kMax;
  double _deltaK;
  double _phi;
  double _tmp_i;
  double _tmp_f;
  int _kas;
public:
  scgle();
  void some_processing();
  ~scgle();
};

  
#endif
