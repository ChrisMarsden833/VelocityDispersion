#ifndef VELOCITYDISPERSION_TESTING_INTEGRATION_H
#define VELOCITYDISPERSION_TESTING_INTEGRATION_H

#include "../integration.h"
#include "math.h"
#include "vector"
#include <stdio.h>

#define PI 3.14159265

float testing_sin(float x, std::vector<float>);

bool test_ARE_sin(void);

bool test_RE_exception(void);

#endif //VELOCITYDISPERSION_TESTING_INTEGRATION_H
