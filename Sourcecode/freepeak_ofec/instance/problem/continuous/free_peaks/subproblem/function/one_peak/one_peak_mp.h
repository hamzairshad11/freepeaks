#ifndef OFEC_FREE_PEAKS_ONE_PEAK_MP_H
#define OFEC_FREE_PEAKS_ONE_PEAK_MP_H

#include "one_peak_base.h"

namespace ofec::free_peaks {
    class OnePeakMP final : public OnePeakBase {
        OFEC_CONCRETE_INSTANCE(OnePeakMP)
    private:
        Real m_radius = 25;
    public:
        void addInputParameters();
        Real evaluate_(Real dummy, size_t var_size) override;
    };
}

#endif
