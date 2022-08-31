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
    cmd.register_command("fit", std::string("0"), "perform a fit for the given cross section.");

    if( cmd.is_enabled("help"))
    {
        return EXIT_SUCCESS;
    }

    cmd.crash_on_missing_mandatory_command();

    double start = cmd.exists("start") ? cmd.as_double("start") : 1.1;
    double end = cmd.exists("end") ? cmd.as_double("end") : 2.0;
    std::size_t steps = cmd.exists("steps") ? static_cast<std::size_t>(cmd.as_int("steps")) : 100ULL;

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

    std::vector<std::string> const nucleon_resonances = {"N1440", "N1520", "N1535"};//, "N1650", "N1675", "N1680", "N1700", "N1710", "N1720", "N1875", "N1880", "N1895", "N1900"};
    std::vector<std::string> const delta_resonances = {"D1232", "D1600", "D1620", "D1700", "D1750"};//, "D1900", "D1905", "D1910", "D1920", "D1930", "D1940"};//, "D1950"};
    std::vector<std::string> const mesons = {"rho0", "rho+"};//, "f0_500"};

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
                Polynomial poly;
                poly.load(FORMAT("widths/{}_width_N_Pi.poly", resonance));
                P[str]->width([&, str, poly](double p2){
                    static const auto threshold = (Proton->mass() + Pi_Plus->mass()) * (Proton->mass() + Pi_Plus->mass());
                    if( p2 < threshold ){
                        return 0.;
                    }
                    auto const sqrt = std::sqrt(p2);
                    auto const ff = P[str]->user_data<FORM_FACTOR_FUNCTION>("form_factor")(P[str], Proton, Pi_Plus, sqrt);
                    return poly(sqrt).real() * ff * ff;
                });

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
                        pim_proton_charge_ex_diagrams.push_back(temp);
                    }
                    else if( P[str]->charge() == 2 ){
                        auto temp = create_diagram(FORMAT("pi_minus proton elastic {} u", P[str]->name()), u_channel, VMP,
                                                   {Proton, Pi_Minus},
                                                   {P[str]},
                                                   {Proton, Pi_Minus}
                        );
                        pim_proton_elastic_diagrams.push_back(temp);
                    }
                }
            }
        }
    }

    if( t_channel_enabled ){
        for( auto const& meson_name : mesons ){
            auto particle = P[meson_name];
            if( particle->charge() == 0 ){
                auto temp = create_diagram(FORMAT("pi_plus proton elastic {} t", particle->name()), t_channel, VMP,
                                           {Proton, Pi_Plus},
                                           {particle},
                                           {Proton, Pi_Plus}
                );
                pip_proton_elastic_diagrams.push_back(temp);
                temp = create_diagram(FORMAT("pi_minus proton elastic {} t", particle->name()), t_channel, VMP,
                                           {Proton, Pi_Minus},
                                           {particle},
                                           {Proton, Pi_Minus}
                );
                pim_proton_elastic_diagrams.push_back(temp);
            }
            else if( particle->charge() == 1 ){
                auto temp = create_diagram(FORMAT("pi_minus proton charge_ex {} u", particle->name()), t_channel, VMP,
                                           {Proton, Pi_Minus},
                                           {particle},
                                           {Neutron, Pi_Zero}
                );
                pim_proton_charge_ex_diagrams.push_back(temp);
            }
        }
    }

    Timer stopwatch;
    std::cout << "config done\n" << std::flush;
    stopwatch.start();

    Feynman_Process scattering_pip_proton_elastic(pip_proton_elastic_diagrams);
    Feynman_Process scattering_pim_proton_elastic(pim_proton_elastic_diagrams);
    Feynman_Process scattering_pim_proton_charge_ex(pim_proton_charge_ex_diagrams);
    scattering_pip_proton_elastic.conversion_factor(1._2mbarn);
    scattering_pim_proton_elastic.conversion_factor(1._2mbarn);
    scattering_pim_proton_charge_ex.conversion_factor(1._2mbarn);

//    couplings.set(coupling_string("N", "Pion", "D1600"), couplings.get(coupling_string("N", "Pion", "D1600")) * 0.6445);
//    couplings.set(coupling_string("N", "Pion", "D1620"), couplings.get(coupling_string("N", "Pion", "D1620")) * 0.8114);
//    couplings.set(coupling_string("N", "Pion", "D1700"), couplings.get(coupling_string("N", "Pion", "D1700")) * 0.0330);
//    couplings.set(coupling_string("N", "Pion", "D1750"), couplings.get(coupling_string("N", "Pion", "D1750")) * 0.18227);

//    couplings.set(coupling_string("N", "Pion", "D1600"), couplings.get(coupling_string("N", "Pion", "D1600")) * 1);
//    couplings.set(coupling_string("N", "Pion", "D1620"), couplings.get(coupling_string("N", "Pion", "D1620")) * 1);
//    couplings.set(coupling_string("N", "Pion", "D1700"), couplings.get(coupling_string("N", "Pion", "D1700")) * 0.07);
//    couplings.set(coupling_string("N", "Pion", "D1750"), couplings.get(coupling_string("N", "Pion", "D1750")) * 0.6688);

    if( cmd.as_string("cross_section") == CMD_CROSS_SECTION_TOTAL )
    {
        std::ofstream file_out("cross_section_total.mat");
        file_out << "pipprotonelasticTotal=";
        scattering_pip_proton_elastic.print_sigma_table(file_out, start, end, steps);
        file_out << ";\npimprotonelasticTotal=";
        scattering_pim_proton_elastic.print_sigma_table(file_out, start, end, steps);
        file_out << ";\npimprotonchargeexTotal=";
        scattering_pim_proton_charge_ex.print_sigma_table(file_out, start, end, steps);
        file_out << ";\n";
    }else
    {
        double sqrt_s = cmd.exists("sqrt_s") ? cmd.as_double("sqrt_s") : 1.5_GeV;
        std::ofstream file_out(FORMAT("cross_section_differential_{}.mat", sqrt_s));
        file_out << "pipprotonelasticDiff=";
        scattering_pip_proton_elastic.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\npimprotonelasticDiff=";
        scattering_pim_proton_elastic.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\npimprotonchargeexDiff=";
        scattering_pim_proton_charge_ex.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\n";

        auto result = scattering_pip_proton_elastic.dsigma_dcos_table_trace(sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));

        std::ofstream file_out2(FORMAT("cross_section_differential_amplitudes{}.json", sqrt_s));

        file_out2 << "{";
        bool outer_first = true;
        for( std::size_t i = 0; i < resonances.size(); ++i )
        {
            if( !s_channel_enabled || !P.exists(resonances[i] + "pp") || P[resonances[i] + "pp"]->charge() != 2 ) continue;
            if( !outer_first ) file_out2 << ",";
            file_out2 << '"' << resonances[i] << '"' << ": [";
            for( std::size_t n_spins = 0; n_spins < 4; ++n_spins )
            {
                if( n_spins > 0 ) file_out2 << ",";
                file_out2 << "[";
                bool first = true;
                for( auto const &[cos, values]: result[n_spins] )
                {
                    if( !first ) file_out2 << ",";
                    file_out2 << "[" << cos << "," << values.first[i].real() << "," << values.first[i].imag()<< "]";
                    first = false;
                }
                file_out2 << "]";
            }
            file_out2 << "]";
            outer_first = false;
        }
        file_out2 << "}";

        std::ofstream file_out3(FORMAT("cross_section_differential_squared_amplitudes{}.json", sqrt_s));

        file_out3 << "{";
        for( std::size_t i = 0; i < resonances.size(); ++i )
        {
            if( i > 0 ) file_out3 << ",";
            file_out3 << '"' << resonances[i] << '"' << ":[";
            for( std::size_t n_spins = 0; n_spins < 4; ++n_spins )
            {
                if( n_spins > 0 ) file_out3 << ",";
                file_out3 << "[";
                bool first = true;
                for( auto const &[cos, values]: result[n_spins] )
                {
                    if( !first ) file_out3 << ",";
                    file_out3 << "[" << cos << "," << values.second[i] << "]";
                    first = false;
                }
                file_out3 << "]";
            }
            file_out3 << "]";
        }
        file_out3 << "}";
        file_out3 << "\n\n";
    }

    stopwatch.stop();
    std::cout << "Time: " << std::setw(5) << stopwatch.time<std::chrono::milliseconds>() / 1000. << "s\n";
    return EXIT_SUCCESS;
}