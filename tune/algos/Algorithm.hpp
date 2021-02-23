
#include <vector>

class Algorithm {
    public: 
        virtual float performAnalysis(std::vector<float> data);
        virtual ~Algorithm() { };
        Algorithm() { };
};