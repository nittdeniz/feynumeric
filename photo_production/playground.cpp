#include <feynumeric/command_line_manager.hpp>
#include <feynumeric/core.hpp>
#include <feynumeric/table.hpp>
#include <feynumeric/messages.hpp>
#include <feynumeric/timer.hpp>


#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"


#include <iostream>
#include <iomanip>

std::string const CMD_PROCESS_ELASTIC_SCATTERING = "elasticscattering";
std::string const CMD_PROCESS_PHOTO_PRODUCTION   = "photoproduction";
std::string const CMD_CHANNEL_S                  = "s";
std::string const CMD_CHANNEL_T                  = "t";
std::string const CMD_CHANNEL_U                  = "u";
std::string const CMD_CHANNEL_C                  = "c";
std::string const CMD_CROSS_SECTION_TOTAL        = "total";
std::string const CMD_CROSS_SECTION_DIFFERENTIAL = "differential";

int main(int argc, char** argv)
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    Command_Line_Manager cmd(argc, argv);

    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("coupling_constants", true, "file with coupling constants");
    cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE,
                         FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING,
                                Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ,
                                Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));
    cmd.register_command("channel", CMD_CHANNEL_S, "which channel to use [s, t, u, c] or any combination.");
    cmd.register_command("process", CMD_PROCESS_PHOTO_PRODUCTION,
                         FORMAT("which process to use: {} {}", CMD_PROCESS_PHOTO_PRODUCTION, CMD_PROCESS_ELASTIC_SCATTERING));
    cmd.register_command("sqrt_s", false, "the energy in the center of mass frame in GeV");
    cmd.register_command("help", false, "list all command line parameters");
    cmd.register_command("cross_section", CMD_CROSS_SECTION_TOTAL,
                         FORMAT("get {} or {} cross section", CMD_CROSS_SECTION_DIFFERENTIAL, CMD_CROSS_SECTION_TOTAL));

    if( cmd.is_enabled("help"))
    {
        return EXIT_SUCCESS;
    }

    cmd.crash_on_missing_mandatory_command();

    auto const &channel = cmd.as_string("channel");
    bool const s_channel_enabled = channel.find('s') != std::string::npos;
    bool const t_channel_enabled = channel.find('t') != std::string::npos;
    bool const u_channel_enabled = channel.find('u') != std::string::npos;
    bool const c_channel_enabled = channel.find('c') != std::string::npos;

    Particle_Manager P(cmd.as_string("particle_file"));
    Particle_Ptr const &Proton = P.get("proton");
    Particle_Ptr const &Neutron = P.get("neutron");
    Particle_Ptr const &Pi_Plus = P.get("pi+");
    Particle_Ptr const &Pi_Minus = P.get("pi-");
    Particle_Ptr const &Pi_Zero = P.get("pi0");

    init_vertices(P, cmd.as_string("coupling_constants"));

    std::string const &form_factor = cmd.as_string("form_factor");

    FORM_FACTOR_FUNCTION ff;
    if( Form_Factor::ff_dict.contains(form_factor))
    {
        ff = Form_Factor::ff_dict[form_factor];
    }else
    {
        critical_error("Unknown form factor");
    }

    status(FORMAT("Form factor: {}", form_factor));

    std::vector<std::string> const nucleon_resonances = {"N1440", "N1520", "N1535", "N1650", "N1675", "N1680", "N1700", "N1710", "N1720", "N1875", "N1880", "N1895", "N1900"};
    std::vector<std::string> const delta_resonances = {"D1232", "D1600", "D1620", "D1700", "D1750", "D1900", "D1905", "D1910", "D1920", "D1930", "D1940", "D1950"};

    std::vector<std::string> resonances;
    resonances.insert(resonances.end(), nucleon_resonances.cbegin(), nucleon_resonances.cend());
    resonances.insert(resonances.end(), delta_resonances.cbegin(), delta_resonances.cend());

    auto pp_string = [](std::string const &p)
    { return FORMAT("{}pp", p); };
    auto p_string = [](std::string const &p)
    { return FORMAT("{}p", p); };
    auto n_string = [](std::string const &p)
    { return FORMAT("{}n", p); };
    auto m_string = [](std::string const &p)
    { return FORMAT("{}m", p); };

    std::vector<Feynman_Diagram_Ptr> pip_proton_elastic_diagrams;
    std::vector<Feynman_Diagram_Ptr> pim_proton_elastic_diagrams;
    std::vector<Feynman_Diagram_Ptr> pim_proton_charge_ex_diagrams;

    for( auto const &resonance : resonances )
    {
        std::vector<std::string> strs = {pp_string(resonance), p_string(resonance), n_string(resonance), m_string(resonance)};

        for( auto const &str: strs )
        {
            if( P.exists(str))
            {
                P[str]->user_data("form_factor", ff);
                auto const &file_name = FORMAT("./data/dyson_factors/{}_{}.txt", resonance, cmd.as_string("form_factor"));
                std::ifstream ifs(file_name);
                if( ifs )
                {
                    Table dyson({});
                    ifs >> dyson;
                    P[str]->user_data("dyson_factors", dyson);
                    P[str]->width([&, str](double p2)
                                  {
                                      auto result = P[str]->user_data<Table>("dyson_factors").interpolate(std::sqrt(p2));
                                      return result;
                                  });
                }else
                {
                    P[str]->width([&, str](double p2)
                                  {
                                      return P[str]->width();
                                  });
                }

                if( s_channel_enabled ){
                    if( P[str]->charge() == 2 ){
                        auto temp = create_diagram(FORMAT("pi_plus proton elastic {} s", P[str]->name()), s_channel, VMP,
                                                   {Proton, Pi_Plus},
                                                   {P[str]},
                                                   {Proton, Pi_Plus}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }else if( P[str]->charge() == 0 ){
                        auto temp = create_diagram(FORMAT("pi_minus proton elastic {} s", P[str]->name()), s_channel, VMP,
                                                   {Proton, Pi_Minus},
                                                   {P[str]},
                                                   {Proton, Pi_Minus}
                        );
                        pim_proton_elastic_diagrams.push_back(temp);
                        temp = create_diagram(FORMAT("pi_minus proton charge ex {} s", P[str]->name()), s_channel, VMP,
                                              {Proton, Pi_Minus},
                                              {P[str]},
                                              {Neutron, Pi_Zero}
                        );
                        pim_proton_charge_ex_diagrams.push_back(temp);
                    }
                }
                if( u_channel_enabled){
                    if( P[str]->charge() == 0 ){
                        auto temp = create_diagram(FORMAT("pi_plus proton elastic {} u", P[str]->name()), u_channel, VMP,
                                               {Proton, Pi_Plus},
                                               {P[str]},
                                               {Proton, Pi_Plus}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }else if( P[str]->charge() == 1 ){
                        auto temp = create_diagram(FORMAT("pi_minus proton charge_ex {} u", P[str]->name()), u_channel, VMP,
                                                   {Proton, Pi_Minus},
                                                   {P[str]},
                                                   {Neutron, Pi_Zero}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }
                    else if( P[str]->charge() == 2 ){
                        auto temp = create_diagram(FORMAT("pi_minus proton elastic {} u", P[str]->name()), u_channel, VMP,
                                                   {Proton, Pi_Minus},
                                                   {P[str]},
                                                   {Proton, Pi_Minus}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }
                }
                /*
                if( t_channel_enabled){
                    if( P[str]->charge() == 0 ){
                        auto temp = create_diagram(FORMAT("pi_plus proton elastic {} u", P[str]->name()), u_channel, VMP,
                                                   {Proton, Pi_Plus},
                                                   {P[str]},
                                                   {Proton, Pi_Plus}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }else if( P[str]->charge() == 1 ){
                        auto temp = create_diagram(FORMAT("pi_minus proton charge_ex {} u", P[str]->name()), u_channel, VMP,
                                                   {Proton, Pi_Plus},
                                                   {P[str]},
                                                   {Neutron, Pi_Zero}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }
                    else if( P[str]->charge() == 2 ){
                        auto temp = create_diagram(FORMAT("pi_minus proton elastic {} u", P[str]->name()), u_channel, VMP,
                                                   {Proton, Pi_Minus},
                                                   {P[str]},
                                                   {Proton, Pi_Minus}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }
                }
                */
            }
        }
    }

    Timer stopwatch;
    stopwatch.start();

    Feynman_Process scattering_pip_proton_elastic(pip_proton_elastic_diagrams);
    Feynman_Process scattering_pim_proton_elastic(pim_proton_elastic_diagrams);
    Feynman_Process scattering_pim_proton_charge_ex(pim_proton_charge_ex_diagrams);
    scattering_pip_proton_elastic.conversion_factor(1._2mbarn);
    scattering_pim_proton_elastic.conversion_factor(1._2mbarn);
    scattering_pim_proton_charge_ex.conversion_factor(1._2mbarn);

    double start = cmd.exists("start") ? cmd.as_double("start") : 1.1;
    double end = cmd.exists("end") ? cmd.as_double("end") : 2.0;
    std::size_t steps = cmd.exists("steps") ? static_cast<std::size_t>(cmd.as_int("steps")) : 100ULL;
    std::ofstream file_out("cross_section.txt");
    file_out << "pipprotonelastic=";
    scattering_pip_proton_elastic.print_sigma_table(file_out, start, end, steps);
    file_out << ";\npimprotonelastic=";
    scattering_pim_proton_elastic.print_sigma_table(file_out, start, end, steps);
    file_out << ";\npimprotonchargeex=";
    scattering_pim_proton_charge_ex.print_sigma_table(file_out, start, end, steps);
    file_out << ";\n";
    stopwatch.stop();
    std::cout << "Time: " << std::setw(5) <<  stopwatch.time<std::chrono::milliseconds>()/1000. << "s\n";

    return EXIT_SUCCESS;
}