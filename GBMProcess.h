#pragma once

#include "TermStructure.h"
#include "date.h"

class GBMProcess
{
public:
	GBMProcess() {}
	GBMProcess(double spot, TermStructure div, YieldTermStructure yield, VolatilityTermStructure vol)
		: spot_(spot), divTermStructure_(div), yieldTermStructure_(yield), volTermStructure_(vol) {}
	~GBMProcess() {};

	double getSpot() { return spot_; }
	void setSpot(double s) { spot_ = s; }
	double getRf(Date d);
	double getDiv(Date d);
	double getVol(Date d);

	double getDiscount(Date d) { return yieldTermStructure_.discount(d); }
	double getForwardRate(Date d1, Date d2) { return yieldTermStructure_.forwardRate(d1, d2); }
	double getVriance(Date d) { return volTermStructure_.variance(d); }

	void setDivTermStructure(TermStructure div) { divTermStructure_ = div; }
	void setYieldTermStructure(YieldTermStructure r) { yieldTermStructure_ = r; }
	void setVolTermStructure(VolatilityTermStructure vol) { volTermStructure_ = vol; }
private:
	double spot_;
	TermStructure divTermStructure_;
	YieldTermStructure yieldTermStructure_;
	VolatilityTermStructure volTermStructure_;
};
