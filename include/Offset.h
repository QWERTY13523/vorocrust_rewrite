#pragma once

#include <string>

namespace Offset {

struct PowerDiagramParams {
    double r_value = 0.001; // Constant R; gradient is 0 < 1. Change in offset.cpp if needed.
    bool r_signed_by_region = true; // Outside: +R, inside: -R.
    double vertex_merge_epsilon = 1e-10;
};

// Reads a seed CSV (with face info), builds base points per seed pair, and exports
// two power diagrams: outside and inside. Returns false on failure.
bool build_power_diagram_from_seed_csv(
    const std::string& seed_csv,
    const std::string& outside_obj,
    const std::string& inside_obj,
    const std::string& interface_obj = "",
    const std::string& basepoints_obj = "power_basepoints.obj",
    const std::string& offsets_obj = "power_offsets.obj",
    const PowerDiagramParams& params = PowerDiagramParams(),
    std::string* error = nullptr
);

} // namespace Offset
