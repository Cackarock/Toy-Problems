//File name: pi-mc-test.cpp
//Created on: 04-September-2020
//Author: ravi

extern "C"{
  #include "../src/pi-mc/pi-mc.h"
}
#include "gtest/gtest.h"
#include <fstream>
#include <numeric>

class MonteCarloTest: public testing::Test{
    public:

    protected:
        //Function to set up space used by all tests
        virtual void SetUp(){
        }
};

/****************************************************************************
 *
 ****************************************************************************/
TEST_F(MonteCarloTest, NoException){
  ASSERT_NO_THROW(s_monte_carlo(10));
}

/****************************************************************************
 *
 ****************************************************************************/
TEST_F(MonteCarloTest, FalseValues){
  ASSERT_NE(10, s_monte_carlo(100));
}


/****************************************************************************
 *
 ****************************************************************************/
TEST_F(MonteCarloTest, PiAccuracy){
  long double pild = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899L; 
  EXPECT_NEAR(pild, s_monte_carlo(100), 0.2);
}
