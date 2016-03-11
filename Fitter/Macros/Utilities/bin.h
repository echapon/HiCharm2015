// a simple template class to store a bin and overload the equality operator

template <typename Type> 
class bin {
   public:
      Type low;
      Type high;

      // constructors
      bin(Type l, Type h) :
         low(l), high(h) { 
         };
      bin(const bin<Type>& b) :
         low(b.low), high(b.high) {
         };
      bin() : 
         low(0), high(0) {
         };

      // destructor
      ~bin(){};

      // operator
      bool operator==(const bin other) const {
         return low==other.low && high==other.high;
      };
      bool operator>(const bin other) const {return true;}
      bool operator<(const bin other) const {return true;}
};

// define a few common uses of the template class
typedef bin<double> binD;
typedef bin<float> binF;
typedef bin<int> binI;

// associate three such bins to make an analysis bin
class anabin {
   public:
      binD rapbin;
      binD ptbin;
      binI centbin;

      // constructors
      anabin(double rapmin, double rapmax, double ptmin, double ptmax, int centmin, int centmax) :
         rapbin(rapmin, rapmax),
         ptbin(ptmin, ptmax),
         centbin(centmin, centmax) {
         };
      anabin(const anabin& a) :
         rapbin(a.rapbin),
         ptbin(a.ptbin),
         centbin(a.centbin) {
         };
      anabin() :
         rapbin(), ptbin(), centbin() {
         };

      // destructor
      ~anabin(){};

      // operator
      bool operator==(const anabin other) const {
         return rapbin==other.rapbin && ptbin==other.ptbin && centbin==other.centbin;
      };
      bool operator>(const anabin other) const {return true;}
      bool operator<(const anabin other) const {return true;}
};
