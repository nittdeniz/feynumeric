#include "direction.hpp"
#include "topologies.hpp"

namespace Feynumeric
{
    const Topology s_channel = {{
                                        {"i0", "v0"},
                                        {"v0", "v1"},
                                        {"i1", "v0"},
                                        {"v1", "o0"},
                                        {"v1", "o1"}
    }};

    const Topology t_channel = {{
                                        {"i0", "v0"},
                                        {"v0", "v1"},
                                        {"i1", "v1"},
                                        {"v0", "o0"},
                                        {"v1", "o1"}
                                }};

    const Topology u_channel = {{
                                        {"i0", "v0"},
                                        {"v0", "v1"},
                                        {"i1", "v1"},
                                        {"v0", "o1"},
                                        {"v1", "o0"}
                                }};

    const Topology u2_channel = {{
                                        {"i0", "v0"},
                                        {"v0", "v1"},
                                        {"v1", "v2"},
                                        {"i1", "v2"},
                                        {"v0", "o1"},
                                        {"v2", "o0"}
                                }};

    const Topology Decay_1_to_2 = {{
                                           {"i0", "v0"},
                                           {"v0", "o0"},
                                           {"v0", "o1"}
    }};

    const Topology Decay_1_to_M2_1 = {{
                                           {"i0", "v0"},
                                           {"v0", "o2"},
                                           {"v0", "v1"},
                                           {"v1", "o0"},
                                           {"v1", "o1"}
                                   }};

    const Topology Decay_1_to_M2_1_cross = {{
                                              {"i0", "v0"},
                                              {"v0", "o1"},
                                              {"v0", "v1"},
                                              {"v1", "o0"},
                                              {"v1", "o2"}
                                      }};

    const Topology Decay_1_to_1_M2 = {{
		                                      {"i0", "v0"},
		                                      {"v0", "o0"},
		                                      {"v0", "v1"},
		                                      {"v1", "o1"},
		                                      {"v1", "o2"}
    }};

	const Topology Decay_1_to_1_M2_cross = {{
			                                  {"i0", "v0"},
			                                  {"v0", "o1"},
			                                  {"v0", "v1"},
			                                  {"v1", "o0"},
			                                  {"v1", "o2"}
	                                  }};
}
