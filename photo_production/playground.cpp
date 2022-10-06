#include <feynumeric/command_line_manager.hpp>
#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>
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


    std::vector<std::string> const nucleon_resonances = {"N1440", "N1520", "N1535", "N1650", "N1675", "N1680", "N1700", "N1710", "N1720"};//, "N1875", "N1880", "N1895", "N1900"};
    std::vector<std::string> const delta_resonances = {"D1232", "D1600", "D1620", "D1700", "D1750"};//, "D1900", "D1905", "D1910", "D1920", "D1930", "D1940"};//, "D1950"};
    std::vector<std::string> const mesons = {"rho0", "rho+"};//, "f0_500"};

    std::vector<std::string> resonances;// = {"Fictional12+_32", "Fictional12-_32", "Fictional32+_32", "Fictional32-_32", "Fictional52+_32", "Fictional52-_32"};
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

    std::vector<Feynman_Diagram_Ptr> photoprod_Pi0P_diagrams;
    std::vector<Feynman_Diagram_Ptr> photoprod_Pi0N_diagrams;
    std::vector<Feynman_Diagram_Ptr> photoprod_PipN_diagrams;
    std::vector<Feynman_Diagram_Ptr> photoprod_PimP_diagrams;

    for( auto const &resonance : resonances ){
        std::vector<std::string> strs = {pp_string(resonance), p_string(resonance), n_string(resonance), m_string(resonance)};
        for( auto const &str: strs ){
            if( P.exists(str)){
                P[str]->user_data("form_factor", ff);
                auto const g = couplings.get(coupling_string(remove_charge_ending(str), "N", "Pion"));
                auto const mR = P[str]->mass();
                auto const mN = Proton->mass();
                auto const mpi = Pi_Plus->mass();
                auto const iso = P[str]->isospin().j() == 0.5? 3 : 1;
                auto const br = P[str]->user_data<double>("branching_N_pi");
                auto form_factor = [&P, str, &Proton, &Pi_Plus, ff](double s){
                    return ff(P[str], Proton, Pi_Plus, std::sqrt(s));
                };

                if( P[str]->spin().j() == 0.5 ){
                    if( P[str]->parity() == 1 ){
                        P[str]->width(
                                [g, mR, mN, mpi, iso, form_factor, br](double s)
                                {
                                    auto q = Feynumeric::momentum(std::sqrt(s), mN, mpi);
                                    auto x = std::sqrt((mN*mN+q*q)*s);

                                    return form_factor(s)/br * iso * g*g*q*(-mN * mR * mpi * mpi + mpi*mpi * x + 2 * q * q * (x + std::sqrt((mpi*mpi+q*q)*s)) ) / (4 * mpi*mpi * M_PI * s);
                                });
                    }else{
                        P[str]->width(
                                [g, mR, mN, mpi, iso, form_factor, br](double s)
                                {
                                    auto q = Feynumeric::momentum(std::sqrt(s), mN, mpi);
                                    auto x = std::sqrt((mN*mN+q*q)*s);
                                    return form_factor(s)/br * iso * g*g*q*(mN * mR * mpi * mpi + mpi*mpi * x + 2 * q * q * (x + std::sqrt((mpi*mpi+q*q)*s)) ) / (4 * mpi*mpi * M_PI * s);
                                });
                    }
                }else if( P[str]->spin().j() == 1.5 ){
                    if( P[str]->parity() == 1 ){
                        P[str]->width(
                                [g, mR, mN, mpi, iso, form_factor, br](double s)
                                {
                                    auto q = Feynumeric::momentum(std::sqrt(s), mN, mpi);
                                    return form_factor(s)/br * iso * g*g * q*q*q * (mN * mR*mR*mR + std::sqrt(s*(mN*mN+q*q)) * s) / (12 * mR * mR * mpi*mpi*mpi*mpi * M_PI);
                                }
                        );
                    }else{
                        P[str]->width(
                                [g, mR, mN, mpi, iso, form_factor, br](double s)
                                {
                                    auto q = Feynumeric::momentum(std::sqrt(s), mN, mpi);
                                    return form_factor(s)/br * iso * g*g * q*q*q * (-mN * mR*mR*mR + std::sqrt(s*(mN*mN+q*q)) * s) / (12 * mR * mR * mpi*mpi*mpi*mpi * M_PI);
                                }
                        );
                    }
                }else if( P[str]->spin().j() == 2.5 ){
                    if( P[str]->parity() == 1 ){
                        P[str]->width(
                                [g, mR, mN, mpi, iso, form_factor, br](double s)
                                {
                                    auto q = Feynumeric::momentum(std::sqrt(s), mN, mpi);
                                    return g*g * form_factor(s)/br * iso * g * g * std::pow(q, 5) * s * (-mN * mR + std::sqrt(s*(mN*mN+q*q)) * s) / (30 * std::pow(mpi, 6) * M_PI);
                                }
                            );
                    }else{
                        P[str]->width(
                                [g, mR, mN, mpi, iso, form_factor, br](double s)
                                {
                                    auto q = Feynumeric::momentum(std::sqrt(s), mN, mpi);
                                    return g*g * form_factor(s)/br * iso * g * g * std::pow(q, 5) * s * (mN * mR + std::sqrt(s*(mN*mN+q*q)) * s) / (30 * std::pow(mpi, 6) * M_PI);
                                }
                        );
                    }
                }else if( P[str]->spin().j() == 3.5 ){
                    if( P[str]->parity() == 1 ){
                        critical_error("width 7/2 not implemented");
                    }else{
                        critical_error("width 7/2 not implemented");
                    }
                }else{
                    critical_error("width not implemented");
                }
                std::cout << FORMAT("Particle: {} Width: {}", str, P[str]->width(mR*mR)) << "\n";
                if( s_channel_enabled ){
                    if( P[str]->charge() == 2 ){
                        auto temp = create_diagram(FORMAT("pi_plus proton elastic {} s", P[str]->name()), s_channel, VMP,
                                                   {Proton, Pi_Plus},
                                                   {P[str]},
                                                   {Proton, Pi_Plus}
                        );
                        pip_proton_elastic_diagrams.push_back(temp);
                    }else if( P[str]->charge() == 1 ){
                        /// Pion Photo Production
                        auto temp = create_diagram(FORMAT("gamma + proton -> pi0 p ({}) s", P[str]->name()), s_channel, VMP,
                                                   {Proton, QED::Photon},
                                                   {P[str]},
                                                   {Proton, Pi_Zero}
                        );
                        photoprod_Pi0P_diagrams.push_back(temp);
                        temp = create_diagram(FORMAT("gamma + proton -> pi+ n ({}) s", P[str]->name()), s_channel, VMP,
                                              {Proton, QED::Photon},
                                              {P[str]},
                                              {Neutron, Pi_Plus}
                        );
                        photoprod_PipN_diagrams.push_back(temp);
                    }
                    else if( P[str]->charge() == 0 )
                    {
                        /// Pion Nucleon Scattering
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

                        /// Pion Photo Production
                        temp = create_diagram(FORMAT("gamma + neutron -> pi0 n ({}) s", P[str]->name()), s_channel, VMP,
                                              {Neutron, QED::Photon},
                                              {P[str]},
                                              {Neutron, Pi_Zero}
                        );
                        photoprod_Pi0N_diagrams.push_back(temp);

                        temp = create_diagram(FORMAT("gamma + neutron -> pi- p ({}) s", P[str]->name()), s_channel, VMP,
                                              {Neutron, QED::Photon},
                                              {P[str]},
                                              {Proton, Pi_Minus}
                        );
                        photoprod_PimP_diagrams.push_back(temp);
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

    if( c_channel_enabled ){
        {
            auto temp = create_diagram("pi_plus proton elastic {} c", c_channel2_2, VMP,
                                       {Proton, Pi_Plus},
                                       {},
                                       {Proton, Pi_Plus}
            );
            pip_proton_elastic_diagrams.push_back(temp);
        }
        {
            auto temp = create_diagram("pi_minus proton elastic {} c", c_channel2_2, VMP,
                                       {Proton, Pi_Minus},
                                       {},
                                       {Proton, Pi_Minus}
            );
            pim_proton_elastic_diagrams.push_back(temp);
        }
        {
            auto temp = create_diagram("pi_minusCE proton elastic {} c", c_channel2_2, VMP,
                                       {Proton, Pi_Minus},
                                       {},
                                       {Neutron, Pi_Zero}
            );
            pim_proton_charge_ex_diagrams.push_back(temp);
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

    Feynman_Process scattering_photoprod_pip_n(photoprod_PipN_diagrams);
    Feynman_Process scattering_photoprod_pim_p(photoprod_PimP_diagrams);
    Feynman_Process scattering_photoprod_pi0_p(photoprod_Pi0P_diagrams);
    Feynman_Process scattering_photoprod_pi0_n(photoprod_Pi0N_diagrams);

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

        file_out << "photoprodPipN=";
        scattering_photoprod_pip_n.print_sigma_table(file_out, start, end, steps);
        file_out << "\nphotoprodPimP=";
        scattering_photoprod_pim_p.print_sigma_table(file_out, start, end, steps);
        file_out << "\nphotoprodPi0P=";
        scattering_photoprod_pi0_p.print_sigma_table(file_out, start, end, steps);
        file_out << "\nphotoprodPi0N=";
        scattering_photoprod_pi0_n.print_sigma_table(file_out, start, end, steps);
    }else
    {
        double sqrt_s = cmd.exists("sqrt_s") ? cmd.as_double("sqrt_s") : 1.5_GeV;
        std::ofstream file_out(FORMAT("mathematica/cross_section_differential_{}.mat", sqrt_s));
        file_out << "pipprotonelasticDiff=";
        scattering_pip_proton_elastic.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\npimprotonelasticDiff=";
        scattering_pim_proton_elastic.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\npimprotonchargeexDiff=";
        scattering_pim_proton_charge_ex.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\n";

        file_out << "pipneutronDiff=";
        scattering_pip_proton_elastic.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\npimprotonelasticDiff=";
        scattering_pim_proton_elastic.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\npimprotonchargeexDiff=";
        scattering_pim_proton_charge_ex.print_dsigma_dcos_table(file_out, sqrt_s, steps);
        file_out << ";\n";

        file_out << "photoprodPipN=";
        scattering_photoprod_pip_n.print_sigma_table(file_out, start, end, steps);
        file_out << "\nphotoprodPimP=";
        scattering_photoprod_pim_p.print_sigma_table(file_out, start, end, steps);
        file_out << "\nphotoprodPi0P=";
        scattering_photoprod_pi0_p.print_sigma_table(file_out, start, end, steps);
        file_out << "\nphotoprodPi0N=";
        scattering_photoprod_pi0_n.print_sigma_table(file_out, start, end, steps);

        std::ofstream file_out2(FORMAT("mathematica/cross_section_pipprotonelastic_differential_amplitudes{}.json", sqrt_s));
        std::ofstream file_out3(FORMAT("mathematica/cross_section_pipprotonelastic_differential_squared_amplitudes{}.json", sqrt_s));
        scattering_pip_proton_elastic.print_dsigma_dcos_table_trace(file_out2, file_out3, sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));
        std::ofstream file_out4(FORMAT("mathematica/cross_section_pimprotonelastic_differential_amplitudes{}.json", sqrt_s));
        std::ofstream file_out5(FORMAT("mathematica/cross_section_pimprotonelastic_differential_squared_amplitudes{}.json", sqrt_s));
        scattering_pim_proton_elastic.print_dsigma_dcos_table_trace(file_out4, file_out5, sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));
        std::ofstream file_out6(FORMAT("mathematica/cross_section_pimprotonCE_differential_amplitudes{}.json", sqrt_s));
        std::ofstream file_out7(FORMAT("mathematica/cross_section_pimprotonCE_differential_squared_amplitudes{}.json", sqrt_s));
        scattering_pim_proton_charge_ex.print_dsigma_dcos_table_trace(file_out6, file_out7, sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));


        std::ofstream file_out10(FORMAT("mathematica/cross_section_photoprod_pipN_differential_amplitudes{}.json", sqrt_s));
        std::ofstream file_out11(FORMAT("mathematica/cross_section_photoprod_pipN_differential_squared_amplitudes{}.json", sqrt_s));
        scattering_photoprod_pip_n.print_dsigma_dcos_table_trace(file_out10, file_out11, sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));

        std::ofstream file_out12(FORMAT("mathematica/cross_section_photoprod_pimP_differential_amplitudes{}.json", sqrt_s));
        std::ofstream file_out13(FORMAT("mathematica/cross_section_photoprod_pimP_differential_squared_amplitudes{}.json", sqrt_s));
        scattering_photoprod_pim_p.print_dsigma_dcos_table_trace(file_out12, file_out13, sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));

        std::ofstream file_out14(FORMAT("mathematica/cross_section_photoprod_pi0p_differential_amplitudes{}.json", sqrt_s));
        std::ofstream file_out15(FORMAT("mathematica/cross_section_photoprod_pi0p_differential_squared_amplitudes{}.json", sqrt_s));
        scattering_photoprod_pi0_p.print_dsigma_dcos_table_trace(file_out14, file_out15, sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));

        std::ofstream file_out16(FORMAT("mathematica/cross_section_photoprod_pi0n_differential_amplitudes{}.json", sqrt_s));
        std::ofstream file_out17(FORMAT("mathematica/cross_section_photoprod_pi0n_differential_squared_amplitudes{}.json", sqrt_s));
        scattering_photoprod_pi0_n.print_dsigma_dcos_table_trace(file_out16, file_out17, sqrt_s, Feynumeric::lin_space(-0.999, 0.999, steps));


    }

    stopwatch.stop();
    std::cout << "Time: " << std::setw(5) << stopwatch.time<std::chrono::milliseconds>() / 1000. << "s\n";
    return EXIT_SUCCESS;
}