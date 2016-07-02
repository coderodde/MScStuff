#ifndef RODDE_CURRENT_TIME_H
#define	RODDE_CURRENT_TIME_H

#include <chrono>
#include <cstdint>

namespace rodde {

    namespace current_time {

        inline uint64_t milliseconds() {
            static std::chrono::system_clock clock;
            return std::chrono::duration_cast<std::chrono::milliseconds>
                    (clock.now().time_since_epoch()).count();
        }
    }
}

#endif	/* RODDE_CURRENT_TIME_H */
