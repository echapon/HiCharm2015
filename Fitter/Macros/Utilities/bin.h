// a simple template class to store a bin and overload the equality operator

template <typename Type> 
class bin {
   public:
      Type low;
      Type high;

      bin(Type l, Type h) :
         low(l), high(h) { 
         };

      ~bin(){};

      bool operator==(bin other) {
         return low==other.low && high==other.high;
      };
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

      anabin(double rapmin, double rapmax, double ptmin, double ptmax, int centmin, int centmax) :
         rapbin(rapmin, rapmax),
         ptbin(ptmin, ptmax),
         centbin(centmin, centmax) {
         };

      ~anabin(){};

      bool operator==(anabin other) {
         return rapbin==other.rapbin && ptbin==other.ptbin && centbin==other.centbin;
      };
};
