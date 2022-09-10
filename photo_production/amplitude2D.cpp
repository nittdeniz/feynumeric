#include <feynumeric/command_line_manager.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/feynman_diagram.hpp>
#include <feynumeric/feynman_process.hpp>
#include <feynumeric/format.hpp>
#include <feynumeric/particle.hpp>
#include <feynumeric/particle_manager.hpp>
#include <feynumeric/timer.hpp>
#include <feynumeric/topologies.hpp>

#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <fstream>
#include <iostream>

std::string const CMD_CHANNEL_S                  = "s";
std::string const CMD_CHANNEL_T                  = "t";
std::string const CMD_CHANNEL_U                  = "u";
std::string const CMD_CHANNEL_C                  = "c";

std::string const CMD_PROCESS_DECAY = "decay";
std::string const CMD_PROCESS_SCATTERING_PIPLUS_PROTON = "piplus_proton";

int main(int argc, char** argv) {
    using namespace Feynumeric;

    Command_Line_Manager cmd(argc, argv);

    cmd.register_command("particle_file", true);
    cmd.register_command("coupling_file", true);
    cmd.register_command("particle", true);
    cmd.register_command("channel", CMD_CHANNEL_S, "which channel to use [s, t, u, c] or any combination.");
    cmd.register_command("process", true);
    cmd.register_command("start_s", std::string("1.08"));
    cmd.register_command("start_m", std::string("1.08"));
    cmd.register_command("end_s", std::string("2.08"));
    cmd.register_command("end_m", std::string("2.08"));
    cmd.register_command("steps_m", std::string("100"));
    cmd.register_command("steps_s", std::string("100"));
    cmd.register_command("steps_c", std::string("100"));
    cmd.register_command("help", false, "list all command line parameters");

    cmd.crash_on_missing_mandatory_command();

    const double start_s = cmd.as_double("start_s");
    const double end_s = cmd.as_double("end_s");
    const double start_m = cmd.as_double("start_m");
    const double end_m = cmd.as_double("end_m");
    std::size_t steps_m = static_cast<std::size_t>(cmd.as_int("steps_m")) + 1;
    double stepsize_m = (end_m-start_m)/(steps_m-1);
    std::size_t steps_s = static_cast<std::size_t>(cmd.as_int("steps_s")) + 1;
    double stepsize_s = (end_s-start_s)/(steps_s-1);
    std::size_t steps_c = static_cast<std::size_t>(cmd.as_int("steps_c")) + 1;
    double stepsize_c = (0.999+0.999)/(steps_c-1);

    auto const &channel = cmd.as_string("channel");
    bool const s_channel_enabled = channel == "s";
    bool const t_channel_enabled = channel == "t";
    bool const u_channel_enabled = channel == "u";
    bool const c_channel_enabled = channel == "c";

    Particle_Manager P(cmd.as_string("particle_file"));
    Particle_Ptr const &Proton = P.get("proton");
    Particle_Ptr const &Neutron = P.get("neutron");
    Particle_Ptr const &Pi_Plus = P.get("pi+");
    Particle_Ptr const &Pi_Minus = P.get("pi-");
    Particle_Ptr const &Pi_Zero = P.get("pi0");

    init_vertices(P, cmd.as_string("coupling_file"));

    Particle_Ptr particle_ptr = P.get(cmd.as_string("particle"));
    particle_ptr->feynman_virtual = [](Feynumeric::Vertex::Edge_Ptr e, Kinematics const& kin){
        return Projector(e, kin);
    };

    particle_ptr->user_data("form_factor", Form_Factor::identity);
    couplings.set(coupling_string(remove_charge_ending(particle_ptr->name()), "N", "Pion"), 1.);

    std::vector<Particle_Ptr> particles(steps_m);
    std::vector<Feynman_Diagram_Ptr> diagrams1(steps_m);
    std::vector<Feynman_Diagram_Ptr> diagrams2(steps_m);
    std::vector<Feynman_Process> processes1;
    std::vector<Feynman_Process> processes2;
    processes1.reserve(steps_m);
    processes2.reserve(steps_m);

    for( std::size_t i = 0; i < steps_m; ++i ){
        particles[i] = std::make_shared<Particle>(*particle_ptr);
        particles[i]->mass(start_m + i * stepsize_m);
    }

    auto const process = cmd.as_string("process");

    if( process  == CMD_PROCESS_DECAY ){
        for( std::size_t i = 0; i < steps_m; ++i ){
            diagrams1[i] = create_diagram("Decay", Decay_1_to_2, VMP,
                                          {particles[i]},
                                          {},
                                          {Proton, Pi_Zero}
            );
            diagrams2[i] = create_diagram("Decay", Decay_1_to_2, VMP,
                                          {particles[i]},
                                          {},
                                          {Neutron, Pi_Plus}
            );
            processes1.push_back(Feynman_Process({diagrams1[i]}));
            processes2.push_back(Feynman_Process({diagrams2[i]}));
        }
        std::map<double, std::map<double, double>> result;

        Timer stopwatch;
        stopwatch.start();

        #pragma omp parallel for
        for( std::size_t i_m = 0; i_m < steps_m; ++i_m ){
            std::cout << "#" << std::flush;
            double const mass = start_m + i_m * stepsize_m;
            double const width = (processes1[i_m].decay_width(mass) + processes2[i_m].decay_width(mass));
            for( std::size_t i_s = 0; i_s < steps_s; ++i_s ){
                double const sqrt_s = start_s + i_s  * stepsize_s;
                result[mass][sqrt_s] = (processes1[i_m].decay_width(sqrt_s) + processes2[i_m].decay_width(sqrt_s))/width;
            }
        }
        stopwatch.stop();
        std::cout << "\nDone: " << stopwatch.time<std::chrono::milliseconds>() / 1000. << "\n";

        std::ofstream out(FORMAT("data/width_{}.json", particle_ptr->name()));
        if( !out ){
            critical_error("Error opening file.");
        }

        out << "[";
        bool inner_first = true;
        for( auto const& [mass, row] : result ){
            for( auto const& [sqrt_s, width] : row ){
                if( !inner_first ){
                    out << ",";
                }
                out << "[" << mass << "," << sqrt_s << "," << width << "]";
                inner_first = false;
            }
        }
        out << "]";
    }else{
        for( std::size_t i = 0; i < steps_m; ++i ){
            if( process == CMD_PROCESS_SCATTERING_PIPLUS_PROTON ) {
                if( s_channel_enabled ) {
                    diagrams1[i] = create_diagram("Scattering", s_channel, VMP,
                                                  {Proton, Pi_Plus},
                                                  {particles[i]},
                                                  {Proton, Pi_Plus}
                    );
                    processes1.push_back(Feynman_Process({diagrams1[i]}));
                }
            }
        }

        struct Result_t{
            double sqrts;
            double cos;
            std::array<Complex, 4> spins;
        };

        std::map<double, std::vector<Result_t>> result;

        for( std::size_t i_m = 0; i_m < steps_m; ++i_m ){
            double const mass = start_m + i_m * stepsize_m;
            result[mass].reserve(steps_c * steps_s);
        }

        Timer stopwatch;
        stopwatch.start();

        #pragma omp parallel for
        for( std::size_t i_m = 0; i_m < steps_m; ++i_m ){
            std::cerr << "#";
            double const mass = start_m + i_m * stepsize_m;
            for( std::size_t i_s = 0; i_s < steps_s; ++i_s ){
                double const sqrt_s = start_s + i_s  * stepsize_s;
                for( std::size_t i_c = 0; i_c < steps_c; ++i_c ){
                    double cos_theta = -0.999 + i_c * stepsize_c;
                    if( cos_theta > 0.999 ) cos_theta = 0.999;
                    auto vec = processes1[i_m].M_costheta(sqrt_s, cos_theta);
                    Result_t temp;
                    temp.sqrts = sqrt_s;
                    temp.cos = cos_theta;
                    std::copy(vec.begin(), vec.end(), temp.spins.begin());
                    result[mass].push_back(temp);
                }
            }
        }
        std::cerr << "\n";

        stopwatch.stop();
        std::cerr << "\nDone: " << stopwatch.time<std::chrono::milliseconds>() / 1000. << "\n";

        std::ofstream out(FORMAT("data/scattering_{}.json", particle_ptr->name()));

        if( !out ){
            critical_error("Error opening file.");
        }

        out << "[";
        bool outer_first = true;
        for( std::size_t i_spin = 0; i_spin < 4; ++i_spin ) {
            if( !outer_first ) out << ",";
            out << "{";
            bool middle_first = true;
            for (auto const& [mass, results]: result) {
                if( !middle_first ) out << ",";
                out << '"' << mass << "\":[";
                bool inner_first = true;
                for (auto result : results ){
                    if( !inner_first ) out << ',';
                    out << "[" << result.sqrts << "," <<result.cos << "," << result.spins[i_spin].real() << "," << result.spins[i_spin].imag() << "]";
                    inner_first = false;
                }
                out << "]";
                middle_first = false;
            }
            out << "}";
            outer_first = false;
        }
        out << "]";
    }

}