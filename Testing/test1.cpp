#define BOOST_TEST_MODULE My Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(test) {
    BOOST_CHECK_EQUAL(3, 3);
}
