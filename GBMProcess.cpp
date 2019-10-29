#include "GBMProcess.h"
#include "TermStructure.h"
#include <stdexcept>

double GBMProcess::getDiv(Date d) {
	return divTermStructure_.value(d);
}

double GBMProcess::getRf(Date d) {
	return yieldTermStructure_.value(d);
}

double GBMProcess::getVol(Date d) {
	return volTermStructure_.value(d);
}
