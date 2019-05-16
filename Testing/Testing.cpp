#define BOOST_TEST_MODULE My Test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(utillity)

    BOOST_AUTO_TEST_CASE( test_case3 )
    {
        BOOST_CHECK( true );
    }

    BOOST_AUTO_TEST_CASE( test_case4 )
    {
        BOOST_CHECK( false );
    }

BOOST_AUTO_TEST_SUITE_END()