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
class anabin : public tuple<binF,binF,binI> {
   public:
      anabin(float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax) :
         tuple<binF,binF,binI> (binF(rapmin,rapmax), binF(ptmin, ptmax), binI(centmin, centmax)) {};
      binF rapbin() const {return get<0>(*this);};
      binF ptbin() const {return get<1>(*this);};
      binI centbin() const {return get<2>(*this);};
      void setrapbin(binF rapbin) {get<0>(*this) = rapbin;};
      void setptbin(binF ptbin) {get<1>(*this) = ptbin;};
      void setcentbin(binI centbin) {get<2>(*this) = centbin;};
};

#endif // #ifndef bin_h
