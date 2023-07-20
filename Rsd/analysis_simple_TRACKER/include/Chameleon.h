#ifndef CHAMELEON_H__
#define CHAMELEON_H__

#include <string>

class Chameleon {
public:
  Chameleon() {};
  explicit Chameleon(const std::string&);
  explicit Chameleon(double);
  explicit Chameleon(const char*);

  Chameleon(const Chameleon&);
  Chameleon& operator=(Chameleon const&);

  Chameleon& operator=(double);
  Chameleon& operator=(std::string const&);

public:
  operator std::string() const;
  operator double     () const;
private:
  std::string value_;
};

#endif
