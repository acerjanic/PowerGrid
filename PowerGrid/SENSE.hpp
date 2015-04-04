//
//  SENSE.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_SENSE_hpp
#define PowerGrid_SENSE_hpp
// We are using two template types at the moment. One for the type of data to be processed (ie Col<cx_double>) and one for the type of G object (ie Gfft<Col<cx_double>>
template <typename T1,typename Tobj>
class SENSE {
public:
    SENSE();
    
    //Class variables go here
    
    //Class constructor
    SENSE(Tobj G) {

    }
    
    //Overloaded operators go here
    
    //Forward transformation is *
    // d is the vector of data of type T1, note it is const, so we don't modify it directly rather return another vector of type T1
    T1 operator*(const T1& d) {
        //Refer to any class variables with this->
        // ie this->variable1 to get to a class variable variable1
    }
    
    T1 operator/(const T1& d) {
        
    }
};

#endif
