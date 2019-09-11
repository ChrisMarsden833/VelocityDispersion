#define BOOST_TEST_MODULE My Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "vector"
#include "utillity.h"
#include "stdlib.h"
#include "stdio.h"

#include "testutil.h"
#include "testinteg.h"


BOOST_AUTO_TEST_CASE(Utility)
{
    test_util();
}

BOOST_AUTO_TEST_CASE(Integration)
{
    test_integ();
}