#include "direction.hpp"
#include "topologies.hpp"

namespace Feynumeric
{
    const Topology s_channel = {{
                    {"i1", "v1"},
                    {"v1", "v2"},
                    {"i2", "v1"},
                    {"v2", "o1"},
                    {"v2", "o2"}
    }};

    const Topology t_channel = {{
                                        {"i1", "v1"},
                                        {"v1", "v2"},
                                        {"i2", "v1"},
                                        {"v1", "o1"},
                                        {"v2", "o2"}
                                }};

    const Topology u_channel = {{
                                        {"i1", "v1"},
                                        {"v1", "v2"},
                                        {"i2", "v1"},
                                        {"v1", "o2"},
                                        {"v2", "o1"}
                                }};

    const Topology Decay_1_to_2 = {{
                                           {"i1", "v1"},
                                           {"v1", "o1"},
                                           {"v1", "o2"}
    }};

    const Topology Decay_1_to_M2_1 = {{
                                           {"i1", "v1"},
                                           {"v1", "o3"},
                                           {"v1", "v2"},
                                           {"v2", "o1"},
                                   }};
}
