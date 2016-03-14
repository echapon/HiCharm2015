#ifndef bin_h
#define bin_h

#include <utility>
#include <tuple>

using namespace std;

// a simple template class to store a bin and overload the equality operator

// define a few common uses of the template class
template <typename T> class bin : public pair<T,T> {
   public:
      bin(T a, T b) : pair<T,T>(a,b) {};
      T low() const {return this->first;}
      T high() const {return this->second;}
};
typedef bin<double> binD;
typedef bin<float> binF;
typedef bin<int> binI;

// associate three such bins to make an analysis bin
class anabin : public tuple<binD,binD,binI> {
   public:
      anabin(double rapmin, double rapmax, double ptmin, double ptmax, int centmin, int centmax) :
         tuple<binD,binD,binI> (binD(rapmin,rapmax), binD(ptmin, ptmax), binI(centmin, centmax)) {};
      binD rapbin() const {return get<0>(*this);};
      binD ptbin() const {return get<1>(*this);};
      binI centbin() const {return get<2>(*this);};
      void setrapbin(binD rapbin) {get<0>(*this) = rapbin;};
      void setptbin(binD ptbin) {get<1>(*this) = ptbin;};
      void setcentbin(binI centbin) {get<2>(*this) = centbin;};
};

#endif // #ifndef bin_h
