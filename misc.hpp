#ifndef MISC_HEAD
#define MISC_HEAD

namespace misc{

    /*logspace where the list returned as a vector*/
    template<class LD>
    std::vector<LD> logspace(LD min, LD max, unsigned int length){
        std::vector<LD> X;
        for(unsigned i = 0; i<length ; ++i){
            X.push_back(std::pow( 10, min + i*(max-min)/( length-1 )));
        }
        return X;
    }

    template<class LD>
    std::vector<LD> linspace(LD min, LD max, unsigned int length){
        std::vector<LD> X;
        for(unsigned i = 0; i<length ; ++i){
            X.push_back( min + i*(max-min)/( length-1 ) );
        }
        return X;
    }

}



#endif