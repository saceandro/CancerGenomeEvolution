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

#include "../util/loglib.hh"
#include "gtest/gtest.h"

#define dblmax (0x7fefffffffffffffULL)

// To use a test fixture, derive a class from testing::Test.
class LogTest : public testing::Test {
 protected:  // You should make the members protected s.t. they can be
             // accessed from sub-classes.
  
  // Declares the variables your tests want to use.
  double_union dm;
  double_union dm_1;
  Log rei;
  Log two;
  Log negative_two;
  Log three;
  Log negative_three;
  Log smallest;
  Log negative_smallest;
  Log dblmin;
  Log negative_dblmin;
  Log thousand;

  // virtual void SetUp() will be called before each test is run.  You
  // should define it if you need to initialize the varaibles.
  // Otherwise, this can be skipped.
  virtual void SetUp() {
    rei = Log(0.0);
    two = Log(2.0);
    negative_two = Log(-2.0);
    three = Log(3.0);
    negative_three = Log(-3.0);
    
    dm.d = DBL_MAX;
    dm_1.c = DBL_MAX_1;
    
    smallest = Log(-dm_1.d, 1);
    negative_smallest = Log(-dm_1.d, -1);
    dblmin = Log(DBL_MIN);
    negative_dblmin = Log(-DBL_MIN);

    thousand = Log(1e300, 1);
  }

  // virtual void TearDown() will be called after each test is run.
  // You should define it if there is cleanup work to do.  Otherwise,
  // you don't have to provide it.
  //
  // virtual void TearDown() {
  // }

  // A helper function that some test uses.
  // static int Double(int n) {
  //   return 2*n;
  // }

  // A helper function for testing Queue::Map().
  // void MapTester(const Queue<int> * q) {
  //   // Creates a new queue, where each element is twice as big as the
  //   // corresponding one in q.
  //   const Queue<int> * const new_q = q->Map(Double);

  //   // Verifies that the new queue has the same size as q.
  //   ASSERT_EQ(q->Size(), new_q->Size());

  //   // Verifies the relationship between the elements of the two queues.
  //   for ( const QueueNode<int> * n1 = q->Head(), * n2 = new_q->Head();
  //         n1 != NULL; n1 = n1->next(), n2 = n2->next() ) {
  //     EXPECT_EQ(2 * n1->element(), n2->element());
  //   }

  //   delete new_q;
  // }
};

// When you have a test fixture, you define a test using TEST_F
// instead of TEST.

TEST_F(LogTest, Constructor)
{
  EXPECT_EQ(rei.get_val(), -DBL_MAX);
  EXPECT_EQ(rei.get_sign(), 0);
  EXPECT_EQ(two.get_val(), log(2.0));
  EXPECT_EQ(two.get_sign(), 1);
  EXPECT_EQ(negative_two.get_val(), log(2.0));
  EXPECT_EQ(negative_two.get_sign(), -1);
  EXPECT_EQ(three.get_val(), log(3.0));
  EXPECT_EQ(three.get_sign(), 1);
  EXPECT_EQ(negative_three.get_val(), log(3.0));
  EXPECT_EQ(negative_three.get_sign(), -1);
  EXPECT_EQ(smallest.get_val(), -dm_1.d);
  EXPECT_EQ(smallest.get_sign(), 1);
  
  EXPECT_EQ(rei.eval(), 0.0);
  EXPECT_DOUBLE_EQ(two.eval(), 2.0);
  EXPECT_DOUBLE_EQ(negative_two.eval(), -2.0);
  EXPECT_DOUBLE_EQ(three.eval(), 3.0);
  EXPECT_DOUBLE_EQ(negative_three.eval(), -3.0);
  EXPECT_EQ(smallest.eval(), 0.0);
  EXPECT_EQ(negative_smallest.eval(), 0.0);

  // Assertions below fail due to the loss of accuracy in return for the wider range of real number represantation.
  // EXPECT_DOUBLE_EQ(dblmin.eval(), DBL_MIN);
  // EXPECT_DOUBLE_EQ(negative_dblmin.eval(), -DBL_MIN);
}


TEST_F(LogTest, Add)
{
  EXPECT_DOUBLE_EQ((two + three).eval(), 5.0);
  EXPECT_DOUBLE_EQ((negative_two + negative_three).eval(), -5.0);
  EXPECT_DOUBLE_EQ((negative_two + three).eval(), 1.0);
  EXPECT_DOUBLE_EQ((two + negative_three).eval(), -1.0);
  
  EXPECT_DOUBLE_EQ((three + two).eval(), 5.0);
  EXPECT_DOUBLE_EQ((negative_three + negative_two).eval(), -5.0);
  EXPECT_DOUBLE_EQ((three + negative_two).eval(), 1.0);
  EXPECT_DOUBLE_EQ((negative_three + two).eval(), -1.0);
  
  EXPECT_EQ((two + negative_two).get_sign(), 0);
  EXPECT_EQ((two + negative_two).get_val(), -DBL_MAX);
  EXPECT_EQ((two + negative_two).eval(), 0.0);
  EXPECT_EQ((three + negative_three).get_sign(), 0);
  EXPECT_EQ((three + negative_three).get_val(), -DBL_MAX);
  EXPECT_EQ((three + negative_three).eval(), 0.0);
  EXPECT_EQ((smallest + negative_smallest).get_sign(), 0);
  EXPECT_EQ((smallest + negative_smallest).get_val(), -DBL_MAX);
  EXPECT_EQ((smallest + negative_smallest).eval(), 0.0);

  EXPECT_EQ((smallest + smallest).get_sign(), 1);
  EXPECT_EQ((smallest + smallest).get_val(), -dm_1.d + log1p(1.0));
  EXPECT_EQ((negative_smallest + negative_smallest).get_sign(), -1);
  EXPECT_EQ((negative_smallest + negative_smallest).get_val(), -dm_1.d + log1p(1.0));

  EXPECT_EQ(((two + smallest) + negative_smallest).get_sign(), 1);
  EXPECT_EQ(((two + smallest) + negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((two + smallest) + negative_smallest).eval(), 2.0);

  EXPECT_EQ(((smallest + two) + negative_smallest).get_sign(), 1);
  EXPECT_EQ(((smallest + two) + negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((smallest + two) + negative_smallest).eval(), 2.0);

  EXPECT_EQ(((two + negative_smallest) + smallest).get_sign(), 1);
  EXPECT_EQ(((two + negative_smallest) + smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((two + negative_smallest) + smallest).eval(), 2.0);

  EXPECT_EQ(((negative_smallest + two) + smallest).get_sign(), 1);
  EXPECT_EQ(((negative_smallest + two) + smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_smallest + two) + smallest).eval(), 2.0);

  EXPECT_EQ(((negative_two + smallest) + negative_smallest).get_sign(), -1);
  EXPECT_EQ(((negative_two + smallest) + negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_two + smallest) + negative_smallest).eval(), -2.0);

  EXPECT_EQ(((smallest + negative_two) + negative_smallest).get_sign(), -1);
  EXPECT_EQ(((smallest + negative_two) + negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((smallest + negative_two) + negative_smallest).eval(), -2.0);

  EXPECT_EQ(((negative_two + negative_smallest) + smallest).get_sign(), -1);
  EXPECT_EQ(((negative_two + negative_smallest) + smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_two + negative_smallest) + smallest).eval(), -2.0);

  EXPECT_EQ(((negative_smallest + negative_two) + smallest).get_sign(), -1);
  EXPECT_EQ(((negative_smallest + negative_two) + smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_smallest + negative_two) + smallest).eval(), -2.0);

  EXPECT_EQ((rei + two).get_sign(), 1);
  EXPECT_EQ((rei + two).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((rei + two).eval(), 2.0);
  
  EXPECT_EQ((two + rei).get_sign(), 1);
  EXPECT_EQ((two + rei).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((two + rei).eval(), 2.0);
  
  EXPECT_EQ((rei + negative_two).get_sign(), -1);
  EXPECT_EQ((rei + negative_two).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((rei + negative_two).eval(), -2.0);

  EXPECT_EQ((negative_two + rei).get_sign(), -1);
  EXPECT_EQ((negative_two + rei).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((negative_two + rei).eval(), -2.0);
}

TEST_F(LogTest, Sub)
{
  EXPECT_DOUBLE_EQ((two - three).eval(), -1.0);
  EXPECT_DOUBLE_EQ((negative_two - negative_three).eval(), 1.0);
  EXPECT_DOUBLE_EQ((negative_two - three).eval(), -5.0);
  EXPECT_DOUBLE_EQ((two - negative_three).eval(), 5.0);
  
  EXPECT_DOUBLE_EQ((three - two).eval(), 1.0);
  EXPECT_DOUBLE_EQ((negative_three - negative_two).eval(), -1.0);
  EXPECT_DOUBLE_EQ((three - negative_two).eval(), 5.0);
  EXPECT_DOUBLE_EQ((negative_three - two).eval(), -5.0);
  
  EXPECT_EQ((two - two).get_sign(), 0);
  EXPECT_EQ((two - two).get_val(), -DBL_MAX);
  EXPECT_EQ((two - two).eval(), 0.0);
  EXPECT_EQ((three - three).get_sign(), 0);
  EXPECT_EQ((three - three).get_val(), -DBL_MAX);
  EXPECT_EQ((three - three).eval(), 0.0);
  EXPECT_EQ((smallest - smallest).get_sign(), 0);
  EXPECT_EQ((smallest - smallest).get_val(), -DBL_MAX);
  EXPECT_EQ((smallest - smallest).eval(), 0.0);

  EXPECT_EQ((negative_two - negative_two).get_sign(), 0);
  EXPECT_EQ((negative_two - negative_two).get_val(), -DBL_MAX);
  EXPECT_EQ((negative_two - negative_two).eval(), 0.0);
  EXPECT_EQ((negative_three - negative_three).get_sign(), 0);
  EXPECT_EQ((negative_three - negative_three).get_val(), -DBL_MAX);
  EXPECT_EQ((negative_three - negative_three).eval(), 0.0);
  EXPECT_EQ((negative_smallest - negative_smallest).get_sign(), 0);
  EXPECT_EQ((negative_smallest - negative_smallest).get_val(), -DBL_MAX);
  EXPECT_EQ((negative_smallest - negative_smallest).eval(), 0.0);

  EXPECT_EQ((smallest - negative_smallest).get_sign(), 1);
  EXPECT_EQ((smallest - negative_smallest).get_val(), -dm_1.d + log1p(1.0));
  EXPECT_EQ((negative_smallest - smallest).get_sign(), -1);
  EXPECT_EQ((negative_smallest - smallest).get_val(), -dm_1.d + log1p(1.0));
  
  EXPECT_EQ(((two - smallest) - negative_smallest).get_sign(), 1);
  EXPECT_EQ(((two - smallest) - negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((two - smallest) - negative_smallest).eval(), 2.0);

  EXPECT_EQ(((smallest - two) - negative_smallest).get_sign(), -1);
  EXPECT_EQ(((smallest - two) - negative_smallest).get_val(), log(2.0) + log1p(-exp(-dm_1.d + log1p(1.0) - log(2.0)))); // not correct ! information loss ! log1p(...) = 0.0 !
  EXPECT_DOUBLE_EQ(((smallest - two) - negative_smallest).eval(), -2.0);

  EXPECT_EQ(((two - negative_smallest) - smallest).get_sign(), 1);
  EXPECT_EQ(((two - negative_smallest) - smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((two - negative_smallest) - smallest).eval(), 2.0);

  EXPECT_EQ(((negative_smallest - two) - smallest).get_sign(), -1);
  EXPECT_EQ(((negative_smallest - two) - smallest).get_val(), log(2.0)  + log1p(-exp(-dm_1.d + log1p(1.0) - log(2.0)))); // not correct ! information loss !
  EXPECT_DOUBLE_EQ(((negative_smallest - two) - smallest).eval(), -2.0);

  EXPECT_EQ(((negative_two - smallest) - negative_smallest).get_sign(), -1);
  EXPECT_EQ(((negative_two - smallest) - negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_two - smallest) - negative_smallest).eval(), -2.0);

  EXPECT_EQ(((smallest - negative_two) - negative_smallest).get_sign(), 1);
  EXPECT_EQ(((smallest - negative_two) - negative_smallest).get_val(), log(2.0)  + log1p(-exp(-dm_1.d + log1p(1.0) - log(2.0)))); // not correct ! information loss !
  EXPECT_DOUBLE_EQ(((smallest - negative_two) - negative_smallest).eval(), 2.0);

  EXPECT_EQ(((negative_two - negative_smallest) - smallest).get_sign(), -1);
  EXPECT_EQ(((negative_two - negative_smallest) - smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_two - negative_smallest) - smallest).eval(), -2.0);

  EXPECT_EQ(((negative_smallest - negative_two) - smallest).get_sign(), 1);
  EXPECT_EQ(((negative_smallest - negative_two) - smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_smallest - negative_two) - smallest).eval(), 2.0);

  
  EXPECT_EQ(((two - smallest) + smallest).get_sign(), 1);
  EXPECT_EQ(((two - smallest) + smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((two - smallest) + smallest).eval(), 2.0);

  EXPECT_EQ(((smallest - two) - smallest).get_sign(), -1);
  EXPECT_EQ(((smallest - two) - smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((smallest - two) - smallest).eval(), -2.0);

  EXPECT_EQ(((two - negative_smallest) + negative_smallest).get_sign(), 1);
  EXPECT_EQ(((two - negative_smallest) + negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((two - negative_smallest) + negative_smallest).eval(), 2.0);

  EXPECT_EQ(((negative_smallest - two) - negative_smallest).get_sign(), -1);
  EXPECT_EQ(((negative_smallest - two) - negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_smallest - two) - negative_smallest).eval(), -2.0);

  EXPECT_EQ(((negative_two - smallest) + smallest).get_sign(), -1);
  EXPECT_EQ(((negative_two - smallest) + smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_two - smallest) + smallest).eval(), -2.0);

  EXPECT_EQ(((smallest - negative_two) - smallest).get_sign(), 1);
  EXPECT_EQ(((smallest - negative_two) - smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((smallest - negative_two) - smallest).eval(), 2.0);

  EXPECT_EQ(((negative_two - negative_smallest) + negative_smallest).get_sign(), -1);
  EXPECT_EQ(((negative_two - negative_smallest) + negative_smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_two - negative_smallest) + negative_smallest).eval(), -2.0);

  EXPECT_EQ(((negative_smallest - negative_two) + smallest).get_sign(), 1);
  EXPECT_EQ(((negative_smallest - negative_two) + smallest).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ(((negative_smallest - negative_two) + smallest).eval(), 2.0);

  
  EXPECT_EQ((rei - two).get_sign(), -1);
  EXPECT_EQ((rei - two).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((rei - two).eval(), -2.0);
  
  EXPECT_EQ((two - rei).get_sign(), 1);
  EXPECT_EQ((two - rei).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((two - rei).eval(), 2.0);
  
  EXPECT_EQ((rei - negative_two).get_sign(), 1);
  EXPECT_EQ((rei - negative_two).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((rei - negative_two).eval(), 2.0);

  EXPECT_EQ((negative_two - rei).get_sign(), -1);
  EXPECT_EQ((negative_two - rei).get_val(), log(2.0));
  EXPECT_DOUBLE_EQ((negative_two - rei).eval(), -2.0);
}

TEST_F(LogTest, Multiply)
{
  EXPECT_DOUBLE_EQ((two * three).eval(), 6.0);
  EXPECT_DOUBLE_EQ((negative_two * negative_three).eval(), 6.0);
  EXPECT_DOUBLE_EQ((negative_two * three).eval(), -6.0);
  EXPECT_DOUBLE_EQ((two * negative_three).eval(), -6.0);
  
  EXPECT_DOUBLE_EQ((three * two).eval(), 6.0);
  EXPECT_DOUBLE_EQ((negative_three * negative_two).eval(), 6.0);
  EXPECT_DOUBLE_EQ((three * negative_two).eval(), -6.0);
  EXPECT_DOUBLE_EQ((negative_three * two).eval(), -6.0);
  
  EXPECT_DOUBLE_EQ((two * negative_two).eval(), -4.0);
  EXPECT_DOUBLE_EQ((three * negative_three).eval(), -9.0);
  EXPECT_EQ((smallest * negative_smallest).get_sign(), 0);
  EXPECT_EQ((smallest * negative_smallest).get_val(), -DBL_MAX);
  EXPECT_EQ((smallest * negative_smallest).eval(), 0.0);

  EXPECT_EQ((smallest * smallest).get_sign(), 0);
  EXPECT_EQ((smallest * smallest).get_val(), -DBL_MAX);
  EXPECT_EQ((smallest * smallest).eval(), 0);
  EXPECT_EQ((negative_smallest * negative_smallest).get_sign(), 0);
  EXPECT_EQ((negative_smallest * negative_smallest).get_val(), -DBL_MAX);
  EXPECT_EQ((negative_smallest * negative_smallest).eval(), 0);

  EXPECT_EQ((rei * two).get_sign(), 0);
  EXPECT_EQ((rei * two).get_val(), -DBL_MAX);
  EXPECT_EQ((rei * two).eval(), 0.0);
  
  EXPECT_EQ((two * rei).get_sign(), 0);
  EXPECT_EQ((two * rei).get_val(), -DBL_MAX);
  EXPECT_EQ((two * rei).eval(), 0.0);
  
  EXPECT_EQ((rei * negative_two).get_sign(), 0);
  EXPECT_EQ((rei * negative_two).get_val(), -DBL_MAX);
  EXPECT_EQ((rei * negative_two).eval(), 0.0);

  EXPECT_EQ((negative_two * rei).get_sign(), 0);
  EXPECT_EQ((negative_two * rei).get_val(), -DBL_MAX);
  EXPECT_EQ((negative_two * rei).eval(), 0.0);
}

TEST_F(LogTest, DIVIDE)
{
  EXPECT_DOUBLE_EQ((two / three).eval(), 2.0/3.0);
  EXPECT_DOUBLE_EQ((negative_two / negative_three).eval(), 2.0/3.0);
  EXPECT_DOUBLE_EQ((negative_two / three).eval(), -2.0/3.0);
  EXPECT_DOUBLE_EQ((two / negative_three).eval(), -2.0/3.0);
  
  EXPECT_DOUBLE_EQ((three / two).eval(), 1.5);
  EXPECT_DOUBLE_EQ((negative_three / negative_two).eval(), 1.5);
  EXPECT_DOUBLE_EQ((three / negative_two).eval(), -1.5);
  EXPECT_DOUBLE_EQ((negative_three / two).eval(), -1.5);
  
  EXPECT_DOUBLE_EQ((two / negative_two).eval(), -1.0);
  EXPECT_DOUBLE_EQ((three / negative_three).eval(), -1.0);

  EXPECT_EQ((smallest / smallest).get_sign(), 1);
  EXPECT_EQ((smallest / smallest).get_val(), 0.0);
  EXPECT_EQ((smallest / smallest).eval(), 1.0);

  EXPECT_EQ((smallest / thousand).get_sign(), 0);
  EXPECT_EQ((smallest / thousand).get_val(), -DBL_MAX);
  EXPECT_EQ((smallest / thousand).eval(), 0.0);

  // EXPECT_EQ(isinf((smallest / three).get_val()), -1);

  EXPECT_EQ((negative_smallest / negative_smallest).get_sign(), 1);
  EXPECT_EQ((negative_smallest / negative_smallest).get_val(), 0.0);
  EXPECT_EQ((negative_smallest / negative_smallest).eval(), 1.0);

  EXPECT_EQ((rei / two).get_sign(), 0);
  EXPECT_EQ((rei / two).get_val(), -DBL_MAX);
  EXPECT_EQ((rei / two).eval(), 0.0);
  
  EXPECT_EQ((rei / negative_two).get_sign(), 0);
  EXPECT_EQ((rei / negative_two).get_val(), -DBL_MAX);
  EXPECT_EQ((rei / negative_two).eval(), 0.0);
}

TEST_F(LogTest, floatingpoint)
{
  EXPECT_EQ(dm.d, DBL_MAX);
  EXPECT_LE(-dm_1.d - 1e300, -dm.d);
  // std::cout << std::scientific << std::setprecision(20);
  // std::cout << dm.d << std::endl;
  // std::cout << dm_1.d << std::endl;
}


// TEST_F(LogTest, HALF_DBOULE)
// {
//   EXPECT_LT(((dblmin / two - dblmin) * two).eval(), 1e-10);
//   EXPECT_LT(((dblmin / two) + (dblmin / two) - dblmin).eval(), 1e-10);
// }

// // Tests the default c'tor.
// TEST_F(QueueTest, DefaultConstructor) {
//   // You can access data in the test fixture here.
//   EXPECT_EQ(0u, q0_.Size());
// }

// // Tests Dequeue().
// TEST_F(QueueTest, Dequeue) {
//   int * n = q0_.Dequeue();
//   EXPECT_TRUE(n == NULL);

//   n = q1_.Dequeue();
//   ASSERT_TRUE(n != NULL);
//   EXPECT_EQ(1, *n);
//   EXPECT_EQ(0u, q1_.Size());
//   delete n;

//   n = q2_.Dequeue();
//   ASSERT_TRUE(n != NULL);
//   EXPECT_EQ(2, *n);
//   EXPECT_EQ(1u, q2_.Size());
//   delete n;
// }

// // Tests the Queue::Map() function.
// TEST_F(QueueTest, Map) {
//   MapTester(&q0_);
//   MapTester(&q1_);
//   MapTester(&q2_);
// }
