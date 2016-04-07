// In this example, we use a more advanced feature of Google Test called
// test fixture.
//
// A test fixture is a place to hold objects and functions shared by
// all tests in a test case.  Using a test fixture avoids duplicating
// the test code necessary to initialize and cleanup those common
// objects for each test.  It is also useful for defining sub-routines
// that your tests need to invoke a lot.
//
// <TechnicalDetails>
//
// The tests share the test fixture in the sense of code sharing, not
// data sharing.  Each test is given its own fresh copy of the
// fixture.  You cannot expect the data modified by one test to be
// passed on to another test, which is a bad idea.
//
// The reason for this design is that tests should be independent and
// repeatable.  In particular, a test should not fail as the result of
// another test's failure.  If one test depends on info produced by
// another test, then the two tests should really be one big test.
//
// The macros for indicating the success/failure of a test
// (EXPECT_TRUE, FAIL, etc) need to know what the current test is
// (when Google Test prints the test result, it tells you which test
// each failure belongs to).  Technically, these macros invoke a
// member function of the Test class.  Therefore, you cannot use them
// in a global function.  That's why you should put test sub-routines
// in a test fixture.
//
// </TechnicalDetails>

#include "../util/enumtree.hh"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

// To use a test fixture, derive a class from testing::Test.
class EnumtreeTest : public testing::Test {
 protected:  // You should make the members protected s.t. they can be
             // accessed from sub-classes.
  
  // Declares the variables your tests want to use.
  hyperparams hpa;
  trees x;
  trees y;
  
  // virtual void SetUp() will be called before each test is run.  You
  // should define it if you need to initialize the varaibles.
  // Otherwise, this can be skipped.
  virtual void SetUp() {
    trees_cons(x, 4);
    hpa.MAX_SUBTYPE = 4;
    hpa.TOTAL_CN = 3;
    hpa.MAX_TREE = x.size();

    y.assign(x.size(), subtypes (hpa.MAX_SUBTYPE + 1, subtype()));
    for (int t=0; t<(int)x.size(); ++t)
      {
        copy(y[t], x[t]);
      }
  }

  // virtual void TearDown() will be called after each test is run.
  // You should define it if there is cleanup work to do.  Otherwise,
  // you don't have to provide it.
  //
  // virtual void TearDown() {
  // }

};

// When you have a test fixture, you define a test using TEST_F
// instead of TEST.

TEST_F(EnumtreeTest, trees_cons)
{
  for (int t=0; t<(int)x.size(); ++t)
    {
      EXPECT_EQ(x[t][0].index, 0);
      EXPECT_EQ(x[t][0].total_cn, 2);
      EXPECT_EQ(x[t][0].variant_cn, 0);
      EXPECT_TRUE(x[t][0].parent == NULL);
      EXPECT_TRUE(x[t][0].above == NULL);
      EXPECT_THAT(x[t][0].children, ::testing::ElementsAre(&x[t][1]));
    }

  
  EXPECT_EQ(x[0][1].index, 1);
  EXPECT_TRUE(x[0][1].parent == &x[0][0]);
  EXPECT_TRUE(x[0][1].above == &x[0][0]);
  EXPECT_THAT(x[0][1].children, ::testing::ElementsAre(&x[0][2]));

  EXPECT_EQ(x[0][2].index, 2);
  EXPECT_TRUE(x[0][2].parent == &x[0][1]);
  EXPECT_TRUE(x[0][2].above == &x[0][1]);
  EXPECT_THAT(x[0][2].children, ::testing::ElementsAre(&x[0][3]));

  EXPECT_EQ(x[0][3].index, 3);
  EXPECT_TRUE(x[0][3].parent == &x[0][2]);
  EXPECT_TRUE(x[0][3].above == &x[0][2]);
  EXPECT_THAT(x[0][3].children, ::testing::ElementsAre(&x[0][4]));

  EXPECT_EQ(x[0][4].index, 4);
  EXPECT_TRUE(x[0][4].parent == &x[0][3]);
  EXPECT_TRUE(x[0][4].above == &x[0][3]);
  EXPECT_THAT(x[0][4].children, ::testing::IsEmpty());

  
  EXPECT_EQ(x[1][1].index, 1);
  EXPECT_TRUE(x[1][1].parent == &x[1][0]);
  EXPECT_TRUE(x[1][1].above == &x[1][0]);
  EXPECT_THAT(x[1][1].children, ::testing::ElementsAre(&x[1][2]));

  EXPECT_EQ(x[1][2].index, 2);
  EXPECT_TRUE(x[1][2].parent == &x[1][1]);
  EXPECT_TRUE(x[1][2].above == &x[1][1]);
  EXPECT_THAT(x[1][2].children, ::testing::ElementsAre(&x[1][3], &x[1][4]));

  EXPECT_EQ(x[1][3].index, 3);
  EXPECT_TRUE(x[1][3].parent == &x[1][2]);
  EXPECT_TRUE(x[1][3].above == &x[1][2]);
  EXPECT_THAT(x[1][3].children, ::testing::IsEmpty());

  EXPECT_EQ(x[1][4].index, 4);
  EXPECT_TRUE(x[1][4].parent == &x[1][2]);
  EXPECT_TRUE(x[1][4].above == &x[1][3]);
  EXPECT_THAT(x[1][4].children, ::testing::IsEmpty());

  
  EXPECT_EQ(x[2][1].index, 1);
  EXPECT_TRUE(x[2][1].parent == &x[2][0]);
  EXPECT_TRUE(x[2][1].above == &x[2][0]);
  EXPECT_THAT(x[2][1].children, ::testing::ElementsAre(&x[2][2], &x[2][4]));

  EXPECT_EQ(x[2][2].index, 2);
  EXPECT_TRUE(x[2][2].parent == &x[2][1]);
  EXPECT_TRUE(x[2][2].above == &x[2][1]);
  EXPECT_THAT(x[2][2].children, ::testing::ElementsAre(&x[2][3]));

  EXPECT_EQ(x[2][3].index, 3);
  EXPECT_TRUE(x[2][3].parent == &x[2][2]);
  EXPECT_TRUE(x[2][3].above == &x[2][2]);
  EXPECT_THAT(x[2][3].children, ::testing::IsEmpty());

  EXPECT_EQ(x[2][4].index, 4);
  EXPECT_TRUE(x[2][4].parent == &x[2][1]);
  EXPECT_TRUE(x[2][4].above == &x[2][2]);
  EXPECT_THAT(x[2][4].children, ::testing::IsEmpty());

  
  EXPECT_EQ(x[3][1].index, 1);
  EXPECT_TRUE(x[3][1].parent == &x[3][0]);
  EXPECT_TRUE(x[3][1].above == &x[3][0]);
  EXPECT_THAT(x[3][1].children, ::testing::ElementsAre(&x[3][2], &x[3][3]));

  EXPECT_EQ(x[3][2].index, 2);
  EXPECT_TRUE(x[3][2].parent == &x[3][1]);
  EXPECT_TRUE(x[3][2].above == &x[3][1]);
  EXPECT_THAT(x[3][2].children, ::testing::IsEmpty());

  EXPECT_EQ(x[3][3].index, 3);
  EXPECT_TRUE(x[3][3].parent == &x[3][1]);
  EXPECT_TRUE(x[3][3].above == &x[3][2]);
  EXPECT_THAT(x[3][3].children, ::testing::ElementsAre(&x[3][4]));

  EXPECT_EQ(x[3][4].index, 4);
  EXPECT_TRUE(x[3][4].parent == &x[3][3]);
  EXPECT_TRUE(x[3][4].above == &x[3][3]);
  EXPECT_THAT(x[3][4].children, ::testing::IsEmpty());


  EXPECT_EQ(x[4][1].index, 1);
  EXPECT_TRUE(x[4][1].parent == &x[4][0]);
  EXPECT_TRUE(x[4][1].above == &x[4][0]);
  EXPECT_THAT(x[4][1].children, ::testing::ElementsAre(&x[4][2], &x[4][3], &x[4][4]));

  EXPECT_EQ(x[4][2].index, 2);
  EXPECT_TRUE(x[4][2].parent == &x[4][1]);
  EXPECT_TRUE(x[4][2].above == &x[4][1]);
  EXPECT_THAT(x[4][2].children, ::testing::IsEmpty());

  EXPECT_EQ(x[4][3].index, 3);
  EXPECT_TRUE(x[4][3].parent == &x[4][1]);
  EXPECT_TRUE(x[4][3].above == &x[4][2]);
  EXPECT_THAT(x[4][3].children, ::testing::IsEmpty());

  EXPECT_EQ(x[4][4].index, 4);
  EXPECT_TRUE(x[4][4].parent == &x[4][1]);
  EXPECT_TRUE(x[4][4].above == &x[4][3]);
  EXPECT_THAT(x[4][4].children, ::testing::IsEmpty());
}

TEST_F(EnumtreeTest, copy)
{
  for (int t=0; t<(int)y.size(); ++t)
    {
      EXPECT_EQ(y[t][0].index, 0);
      EXPECT_EQ(y[t][0].total_cn, 2);
      EXPECT_EQ(y[t][0].variant_cn, 0);
      EXPECT_TRUE(y[t][0].parent == NULL);
      EXPECT_TRUE(y[t][0].above == NULL);
      EXPECT_THAT(y[t][0].children, ::testing::ElementsAre(&y[t][1]));
    }

  
  EXPECT_EQ(y[0][1].index, 1);
  EXPECT_TRUE(y[0][1].parent == &y[0][0]);
  EXPECT_TRUE(y[0][1].above == &y[0][0]);
  EXPECT_THAT(y[0][1].children, ::testing::ElementsAre(&y[0][2]));

  EXPECT_EQ(y[0][2].index, 2);
  EXPECT_TRUE(y[0][2].parent == &y[0][1]);
  EXPECT_TRUE(y[0][2].above == &y[0][1]);
  EXPECT_THAT(y[0][2].children, ::testing::ElementsAre(&y[0][3]));

  EXPECT_EQ(y[0][3].index, 3);
  EXPECT_TRUE(y[0][3].parent == &y[0][2]);
  EXPECT_TRUE(y[0][3].above == &y[0][2]);
  EXPECT_THAT(y[0][3].children, ::testing::ElementsAre(&y[0][4]));

  EXPECT_EQ(y[0][4].index, 4);
  EXPECT_TRUE(y[0][4].parent == &y[0][3]);
  EXPECT_TRUE(y[0][4].above == &y[0][3]);
  EXPECT_THAT(y[0][4].children, ::testing::IsEmpty());

  
  EXPECT_EQ(y[1][1].index, 1);
  EXPECT_TRUE(y[1][1].parent == &y[1][0]);
  EXPECT_TRUE(y[1][1].above == &y[1][0]);
  EXPECT_THAT(y[1][1].children, ::testing::ElementsAre(&y[1][2]));

  EXPECT_EQ(y[1][2].index, 2);
  EXPECT_TRUE(y[1][2].parent == &y[1][1]);
  EXPECT_TRUE(y[1][2].above == &y[1][1]);
  EXPECT_THAT(y[1][2].children, ::testing::ElementsAre(&y[1][3], &y[1][4]));

  EXPECT_EQ(y[1][3].index, 3);
  EXPECT_TRUE(y[1][3].parent == &y[1][2]);
  EXPECT_TRUE(y[1][3].above == &y[1][2]);
  EXPECT_THAT(y[1][3].children, ::testing::IsEmpty());

  EXPECT_EQ(y[1][4].index, 4);
  EXPECT_TRUE(y[1][4].parent == &y[1][2]);
  EXPECT_TRUE(y[1][4].above == &y[1][3]);
  EXPECT_THAT(y[1][4].children, ::testing::IsEmpty());

  
  EXPECT_EQ(y[2][1].index, 1);
  EXPECT_TRUE(y[2][1].parent == &y[2][0]);
  EXPECT_TRUE(y[2][1].above == &y[2][0]);
  EXPECT_THAT(y[2][1].children, ::testing::ElementsAre(&y[2][2], &y[2][4]));

  EXPECT_EQ(y[2][2].index, 2);
  EXPECT_TRUE(y[2][2].parent == &y[2][1]);
  EXPECT_TRUE(y[2][2].above == &y[2][1]);
  EXPECT_THAT(y[2][2].children, ::testing::ElementsAre(&y[2][3]));

  EXPECT_EQ(y[2][3].index, 3);
  EXPECT_TRUE(y[2][3].parent == &y[2][2]);
  EXPECT_TRUE(y[2][3].above == &y[2][2]);
  EXPECT_THAT(y[2][3].children, ::testing::IsEmpty());

  EXPECT_EQ(y[2][4].index, 4);
  EXPECT_TRUE(y[2][4].parent == &y[2][1]);
  EXPECT_TRUE(y[2][4].above == &y[2][2]);
  EXPECT_THAT(y[2][4].children, ::testing::IsEmpty());

  
  EXPECT_EQ(y[3][1].index, 1);
  EXPECT_TRUE(y[3][1].parent == &y[3][0]);
  EXPECT_TRUE(y[3][1].above == &y[3][0]);
  EXPECT_THAT(y[3][1].children, ::testing::ElementsAre(&y[3][2], &y[3][3]));

  EXPECT_EQ(y[3][2].index, 2);
  EXPECT_TRUE(y[3][2].parent == &y[3][1]);
  EXPECT_TRUE(y[3][2].above == &y[3][1]);
  EXPECT_THAT(y[3][2].children, ::testing::IsEmpty());

  EXPECT_EQ(y[3][3].index, 3);
  EXPECT_TRUE(y[3][3].parent == &y[3][1]);
  EXPECT_TRUE(y[3][3].above == &y[3][2]);
  EXPECT_THAT(y[3][3].children, ::testing::ElementsAre(&y[3][4]));

  EXPECT_EQ(y[3][4].index, 4);
  EXPECT_TRUE(y[3][4].parent == &y[3][3]);
  EXPECT_TRUE(y[3][4].above == &y[3][3]);
  EXPECT_THAT(y[3][4].children, ::testing::IsEmpty());


  EXPECT_EQ(y[4][1].index, 1);
  EXPECT_TRUE(y[4][1].parent == &y[4][0]);
  EXPECT_TRUE(y[4][1].above == &y[4][0]);
  EXPECT_THAT(y[4][1].children, ::testing::ElementsAre(&y[4][2], &y[4][3], &y[4][4]));

  EXPECT_EQ(y[4][2].index, 2);
  EXPECT_TRUE(y[4][2].parent == &y[4][1]);
  EXPECT_TRUE(y[4][2].above == &y[4][1]);
  EXPECT_THAT(y[4][2].children, ::testing::IsEmpty());

  EXPECT_EQ(y[4][3].index, 3);
  EXPECT_TRUE(y[4][3].parent == &y[4][1]);
  EXPECT_TRUE(y[4][3].above == &y[4][2]);
  EXPECT_THAT(y[4][3].children, ::testing::IsEmpty());

  EXPECT_EQ(y[4][4].index, 4);
  EXPECT_TRUE(y[4][4].parent == &y[4][1]);
  EXPECT_TRUE(y[4][4].above == &y[4][3]);
  EXPECT_THAT(y[4][4].children, ::testing::IsEmpty());
}
