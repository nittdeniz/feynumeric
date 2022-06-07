#include <feynumeric/feynumeric.hpp>
#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include <omp.h>

#define D1232 0
#define N1440 0
#define N1520 0
#define N1535 1
#define D1600 1

int main(int argc, char** argv)
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using namespace std::placeholders;

    Command_Line_Manager cmd(argc, argv);
    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("start", std::string("1.0"), "starting point");
    cmd.register_command("end", std::string("3.0"), "end value");
    cmd.register_command("steps", std::string("200"), "steps");
    cmd.register_command("form_factor", CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", CMD_FORM_FACTOR_NONE, CMD_FORM_FACTOR_CASSING, CMD_FORM_FACTOR_CUTKOSKY, CMD_FORM_FACTOR_MANLEY, CMD_FORM_FACTOR_MONIZ, CMD_FORM_FACTOR_BREIT_WIGNER));

    cmd.crash_on_missing_mandatory_command();

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Neutron = P.get("neutron");
    auto const& Pi_Zero = P.get("pi0");
    auto const& Pi_Minus = P.get("pi-");
    auto const& Pi_Plus = P.get("pi+");

    auto breit_wigner_modified = [](Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, double E, double M){
        double const l = .7;
        double const l4 = std::pow(l, 4);
        return l4/(std::pow(E*E-M*M, 2)+l4);
    };

    auto gaussian_modified = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E, double M)
    {
        return std::exp(-std::pow(E*E - M * M, 2)/(std::pow(0.7, 4)));
    };

    auto rayleigh_modified = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E, double M)
    {
        auto l = 0.61;
        auto l4 = std::pow(l, 4);
        auto kin_limit =0.95*( N->mass() + pi->mass());
        return (M - kin_limit)/(E-kin_limit) *  std::exp(-std::pow(E*E - M * M, 2)/l4);
    };

    auto inverse_gaussian_modified = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E, double M)
    {
        auto l = 10.;
        return std::pow(E/M, -1.5) * std::exp(-l * std::pow(E-M, 2) / (M * E));
    };

    std::string const ff_str = cmd.as_string("form_factor");
    FORM_FACTOR_FUNCTION ff;
    if( ff_dict.contains(ff_str) ){
        ff = ff_dict[ff_str];
    }
    else{
        critical_error("Unknown form factor");
    }

    status(FORMAT("Form factor: {}", ff_str));

    for( auto& [key, particle] : P ){
        auto k = key;
        auto p = particle;
        particle->user_data("form_factor", ff);
    }

    init_vertices(P);

    double const start = cmd.as_double("start");
    double const end   = cmd.as_double("end");
    int    const steps = cmd.as_int("steps");

    std::vector<double> values;
    values.reserve(steps+1);
    for( int i = 0; i <= steps; ++i ){
        values.push_back(start + (end - start) / steps * i);
    }

    std::cout << FORMAT("start: {} end: {} steps: {}\n", start, end, steps);
    #if D1232
    {/// D1232 -> Npi
        std::cout << "D1232\n";
        auto dummypp = std::make_shared<Particle>(*P.get("D1232pp"));
//        FORM_FACTOR_FUNCTION breit_wigner_d1232 = std::bind(inverse_gaussian_modified, _1, _2, _3, _4, P.get("D1232pp")->mass());
//        dummypp->user_data("form_factor", identity);
        auto decay_1 = create_diagram(FORMAT("decay {} to proton pi+", dummypp->name()), Decay_1_to_2, VMP,
                                      {dummypp},
                                      {},
                                      {Proton, Pi_Plus}
        );
        Feynman_Process decay({decay_1});
        std::map<double, double> dyson_factor;
        values.back() = dummypp->mass();
        for( auto const& value : values ){
            dummypp->mass(value);
            auto temp = decay.decay_width();// /Gamma0;
            dyson_factor[value] = std::isnan(temp)? 0 : temp;
        }
		Table(dyson_factor).write(FORMAT("data/dyson_factors/D1232_{}.txt", ff_str));
    }
    #endif
    #if N1440
    {   /// N1440 -> Npi and N1440 -> Npipi
        std::cout << "N1440\n";
        auto dummy = std::make_shared<Particle>(*P.get("N1440p"));
        dummy->user_data("form_factor", identity);

        auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_2, VMP,
                                      {dummy},
                                      {},
                                      {Proton, Pi_Zero}
        );
        auto diagram_pi2 = create_diagram(FORMAT("decay {} to neutron pi+", dummy->name()), Decay_1_to_2, VMP,
                                          {dummy},
                                          {},
                                          {Neutron, Pi_Plus}
        );

        auto diagram_pipi_1 = create_diagram(FORMAT("decay {} to proton pi+ pi-", dummy->name()), Decay_1_to_M2_1, VMP,
                                        {dummy},
                                        {P.get("D1232pp")},
                                        {Proton, Pi_Plus, Pi_Minus}
        );
        auto diagram_pipi_2 = create_diagram(FORMAT("decay {} to proton pi+ pi-", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                             {dummy},
                                             {P.get("D1232n")},
                                             {Proton, Pi_Plus, Pi_Minus}
        );
        auto diagram_pipi_3 = create_diagram(FORMAT("decay {} to proton f0_500", dummy->name()), Decay_1_to_1_M2, VMP,
                                             {dummy},
                                             {P.get("f0_500")},
                                             {Proton, Pi_Plus, Pi_Minus}
        );
        auto diagram_pipi_4 = create_diagram(FORMAT("decay {} to proton pi+ pi-", dummy->name()), Decay_1_to_M2_1, VMP,
                                             {dummy},
                                             {P.get("D1232p")},
                                             {Proton, Pi_Zero, Pi_Zero}
        );
        auto diagram_pipi_5 = create_diagram(FORMAT("decay {} to proton pi+ pi-", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                             {dummy},
                                             {P.get("D1232p")},
                                             {Proton, Pi_Zero, Pi_Zero}
        );
        auto diagram_pipi_6 = create_diagram(FORMAT("decay {} to proton f0_500", dummy->name()), Decay_1_to_1_M2, VMP,
                                             {dummy},
                                             {P.get("f0_500")},
                                             {Proton, Pi_Zero, Pi_Zero}
        );

        auto diagram_pipi_7 = create_diagram(FORMAT("decay {} to proton pi+ pi0", dummy->name()), Decay_1_to_M2_1, VMP,
                                             {dummy},
                                             {P.get("D1232p")},
                                             {Neutron, Pi_Plus, Pi_Zero}
        );
        auto diagram_pipi_8 = create_diagram(FORMAT("decay {} to proton pi+ pi0", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                             {dummy},
                                             {P.get("D1232n")},
                                             {Neutron, Pi_Plus, Pi_Zero}
        );

        Feynman_Process decay_pi1({diagram_pi1});
        Feynman_Process decay_pi2({diagram_pi2});
        Feynman_Process decay_pipi1({diagram_pipi_1, diagram_pipi_2, diagram_pipi_3});
        Feynman_Process decay_pipi2({diagram_pipi_4, diagram_pipi_5, diagram_pipi_6});
        Feynman_Process decay_pipi3({diagram_pipi_7, diagram_pipi_8});

        std::map<double, double> dyson_factor;
        values.back() = dummy->mass();
        std::size_t i{0};
        for( auto const& value : values ){
            dummy->mass(value);
            auto temp1 = decay_pi1.decay_width();
            auto temp2 = decay_pi2.decay_width();
            temp1 = std::isnan(temp1) ? 0  : temp1;
            temp2 = std::isnan(temp2) ? 0  : temp2;

            auto temp3 = decay_pipi1.decay_width();
            auto temp4 = decay_pipi2.decay_width();
            auto temp5 = decay_pipi3.decay_width();
            temp3 = std::isnan(temp3) ? 0  : temp3;
            temp4 = std::isnan(temp4) ? 0  : temp4;
            temp5 = std::isnan(temp5) ? 0  : temp5;

            dyson_factor[value] = temp1 + temp2 + temp3 + temp4 + temp5;
            std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << FORMAT("{} {} {} / {}: {} - {}\n", std::ctime(&time), dummy->name(), i++, steps+1, value, dyson_factor[value]) << std::flush;
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/N1440_{}.txt", ff_str));
    }
    #endif
    #if N1520
	{/// N1520 -> Npi and N1520 -> Npipi
	    std::cout << "N1520\n";
		auto dummy = std::make_shared<Particle>(*P.get("N1520p"));
		dummy->user_data("form_factor", identity);
		auto diagram_pi = create_diagram(FORMAT("decay {} to proton pi0 pi0", dummy->name()), Decay_1_to_2, VMP,
		                                 {dummy},
		                                 {},
		                                 {Proton, Pi_Zero}
		);
		auto diagram_pipi_1 = create_diagram(FORMAT("decay {} to proton pi+ pi-", dummy->name()), Decay_1_to_M2_1, VMP,
		                                     {dummy},
		                                     {P.get("D1232pp")},
		                                     {Proton, Pi_Plus, Pi_Minus}
		);
		auto diagram_pipi_2 = create_diagram(FORMAT("decay {} to proton pi+ pi-", dummy->name()), Decay_1_to_M2_1_cross, VMP,
		                                     {dummy},
		                                     {P.get("D1232n")},
		                                     {Proton, Pi_Plus, Pi_Minus}
		);
		auto diagram_pipi_3 = create_diagram(FORMAT("decay {} to proton rho", dummy->name()), Decay_1_to_1_M2, VMP,
		                                     {dummy},
		                                     {P.get("rho0")},
		                                     {Proton, Pi_Plus, Pi_Minus}
		);


		Feynman_Process decay_pi({diagram_pi});
		Feynman_Process decay_pipi({diagram_pipi_1, diagram_pipi_2, diagram_pipi_3});

        std::vector<Feynman_Process> copies_pi;
        std::vector<Feynman_Process> copies_pipi;
        copies_pi.reserve(steps);
        copies_pipi.reserve(steps);
        for( std::size_t i = 0; i <= steps; ++i ){
            copies_pi.emplace_back(decay_pi);
            copies_pipi.emplace_back(decay_pipi);
        }

        std::map<double, double> dyson_factor;

        values.back() = dummy->mass();
        std::size_t i{0};
        for( auto const& value : values ){
            dummy->mass(value);
            auto temp1 = copies_pi[i].decay_width();// / Gamma0_1;
            temp1 = std::isnan(temp1) ? 0  : temp1;
            auto temp2 = copies_pipi[i].decay_width();// / Gamma0_2;
            temp2 = std::isnan(temp2) ? 0  : temp2;
            dyson_factor[value] = temp1 + temp2;
            std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << FORMAT("{} {} {} / {}: {} - {}\n", std::ctime(&time), dummy->name(), i++, steps+1, value, dyson_factor[value]) << std::flush;
        }
		Table(dyson_factor).write(FORMAT("data/dyson_factors/N1520_{}.txt", ff_str));
	}
    #endif
    #if N1535
    {
        std::cout << "N1535\n";
        auto dummy = std::make_shared<Particle>(*P.get("N1535p"));
        dummy->user_data("form_factor", identity);
        /** N1535 -> Npi **/
        auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_2, VMP,
                                          {dummy},
                                          {},
                                          {Proton, Pi_Zero}
        );
        auto diagram_pi2 = create_diagram(FORMAT("decay {} to neutron pi+", dummy->name()), Decay_1_to_2, VMP,
                                          {dummy},
                                          {},
                                          {Neutron, Pi_Plus}
        );
        /** N1535 -> Neta **/
        auto diagram_eta = create_diagram(FORMAT("decay {} to proton eta", dummy->name()), Decay_1_to_2, VMP,
                                          {dummy},
                                          {},
                                          {Proton, P.get("eta")}
        );
        /** N1535 -> Npipi **/
        auto diagram_pip_pim_1 = create_diagram(FORMAT("decay {} -> D1232 -> p pi+pi-", dummy->name()), Decay_1_to_M2_1, VMP,
                                                {dummy},
                                                {P.get("D1232pp")},
                                                {Proton, Pi_Plus, Pi_Minus}
        );
        auto diagram_pip_pim_2 = create_diagram(FORMAT("decay {} -> D1232 -> p pi+pi-", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("D1232n")},
                                                {Proton, Pi_Plus, Pi_Minus}
        );
        auto diagram_pip_pim_3 = create_diagram(FORMAT("decay {} -> sigma -> p pi+pi-", dummy->name()), Decay_1_to_1_M2, VMP,
                                                {dummy},
                                                {P.get("f0_500")},
                                                {Proton, Pi_Plus, Pi_Minus}
        );
        auto diagram_pip_pim_4 = create_diagram(FORMAT("decay {} -> D1232 -> p pi+pi-", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("N1440n")},
                                                {Proton, Pi_Plus, Pi_Minus}
        );



        auto diagram_pip_pi0_1 = create_diagram(FORMAT("decay {} -> D1232p -> p pi+pi-", dummy->name()), Decay_1_to_M2_1, VMP,
                                                {dummy},
                                                {P.get("D1232p")},
                                                {Neutron, Pi_Plus, Pi_Zero}
        );
        auto diagram_pip_pi0_2 = create_diagram(FORMAT("decay {} -> D1232n -> p pi+pi-", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("D1232n")},
                                                {Neutron, Pi_Plus, Pi_Zero}
        );
        auto diagram_pip_pi0_3 = create_diagram(FORMAT("decay {} -> N1440p -> n pi+pi0", dummy->name()), Decay_1_to_M2_1, VMP,
                                                {dummy},
                                                {P.get("N1440p")},
                                                {Neutron, Pi_Plus, Pi_Zero}
        );
        auto diagram_pip_pi0_4 = create_diagram(FORMAT("decay {} -> N1440n -> p pi+pi0", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("N1440n")},
                                                {Neutron, Pi_Plus, Pi_Zero}
        );



        auto diagram_pi0_pi0_1 = create_diagram(FORMAT("decay {} -> D1232 -> p pi0pi0", dummy->name()), Decay_1_to_M2_1, VMP,
                                                {dummy},
                                                {P.get("D1232p")},
                                                {Proton, Pi_Zero, Pi_Zero}
        );
        auto diagram_pi0_pi0_2 = create_diagram(FORMAT("decay {} -> D1232 -> p pi0pi0", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("D1232p")},
                                                {Proton, Pi_Zero, Pi_Zero}
        );
        auto diagram_pi0_pi0_3 = create_diagram(FORMAT("decay {} -> N1440 -> p pi0pi0", dummy->name()), Decay_1_to_M2_1, VMP,
                                                {dummy},
                                                {P.get("N1440p")},
                                                {Proton, Pi_Zero, Pi_Zero}
        );
        auto diagram_pi0_pi0_4 = create_diagram(FORMAT("decay {} -> N1440 -> p pi0pi0", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("N1440p")},
                                                {Proton, Pi_Zero, Pi_Zero}
        );
        auto diagram_pi0_pi0_5 = create_diagram(FORMAT("decay {} -> sigma -> p pi0pi0", dummy->name()), Decay_1_to_1_M2, VMP,
                                                {dummy},
                                                {P.get("f0_500")},
                                                {Proton, Pi_Zero, Pi_Zero}
        );


        std::vector<Feynman_Process> processes;
        processes.push_back(Feynman_Process({diagram_pi1}));
        processes.push_back(Feynman_Process({diagram_pi2}));
        processes.push_back(Feynman_Process({diagram_eta}));
        processes.push_back(Feynman_Process({diagram_pip_pim_1, diagram_pip_pim_2, diagram_pip_pim_3, diagram_pip_pim_4}));
        processes.push_back(Feynman_Process({diagram_pi0_pi0_1, diagram_pi0_pi0_2, diagram_pi0_pi0_3, diagram_pi0_pi0_4, diagram_pi0_pi0_5}));
        processes.push_back(Feynman_Process({diagram_pip_pi0_1, diagram_pip_pi0_2, diagram_pip_pi0_3, diagram_pip_pi0_4}));

        std::map<double, double> dyson_factor;
        values.back() = dummy->mass();
        std::size_t i{0};
        for( auto const& value : values ){
            dummy->mass(value);
            for( auto& process : processes ){
                auto temp = process.decay_width();
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << FORMAT("{} {} {} / {}: {} - {}\n", std::ctime(&time), dummy->name(), i++, steps+1, value, dyson_factor[value]) << std::flush;
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/N1535_{}.txt", ff_str));
    }
#if D1600
    {
        std::cout << "D1600\n";
        auto dummy = std::make_shared<Particle>(*P.get("D1600pp"));
        dummy->user_data("form_factor", identity);
        /** D1600 -> Npi **/
        auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                         {dummy},
                                         {},
                                         {Proton, Pi_Plus}
        );
        /** D1600 -> Npipi **/
        auto diagram_pip_pin_1 = create_diagram(FORMAT("decay {} -> D1232 -> p pi+pi0 1", dummy->name()), Decay_1_to_M2_1, VMP,
                                            {dummy},
                                            {P.get("D1232pp")},
                                            {Proton, Pi_Plus, Pi_Zero}
        );
        auto diagram_pip_pin_2 = create_diagram(FORMAT("decay {} -> D1232 -> p pi+pi0 2", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("D1232p")},
                                                {Proton, Pi_Plus, Pi_Zero}
        );
        auto diagram_pip_pin_3 = create_diagram(FORMAT("decay {} -> N1440 -> p pi+pi0 3", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("N1440p")},
                                                {Proton, Pi_Plus, Pi_Zero}
        );

        auto diagram_pip_pip_1 = create_diagram(FORMAT("decay {} -> D1232 -> p pi+pi+ 1", dummy->name()), Decay_1_to_M2_1, VMP,
                                                {dummy},
                                                {P.get("D1232p")},
                                                {Neutron, Pi_Plus, Pi_Plus}
        );
        auto diagram_pip_pip_2 = create_diagram(FORMAT("decay {} -> D1232 -> p pi+pi+ 2", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("D1232p")},
                                                {Neutron, Pi_Plus, Pi_Plus}
        );
        auto diagram_pip_pip_3 = create_diagram(FORMAT("decay {} -> N1440 -> p pi+pi+ 3", dummy->name()), Decay_1_to_M2_1, VMP,
                                                {dummy},
                                                {P.get("N1440p")},
                                                {Neutron, Pi_Plus, Pi_Plus}
        );
        auto diagram_pip_pip_4 = create_diagram(FORMAT("decay {} -> N1440 -> p pi+pi+ 4", dummy->name()), Decay_1_to_M2_1_cross, VMP,
                                                {dummy},
                                                {P.get("N1440p")},
                                                {Neutron, Pi_Plus, Pi_Plus}
        );



        std::vector<Feynman_Process> processes;
        processes.push_back(Feynman_Process({diagram_pi1}));
        processes.push_back(Feynman_Process({diagram_pip_pin_1, diagram_pip_pin_2, diagram_pip_pin_3}));
        processes.push_back(Feynman_Process({diagram_pip_pip_1, diagram_pip_pip_2, diagram_pip_pip_3, diagram_pip_pip_4}));

        std::map<double, double> dyson_factor;
        values.back() = dummy->mass();
        std::size_t i{0};
        for( auto const& value : values ){
            dummy->mass(value);
            for( auto& process : processes ){
                auto temp = process.decay_width();
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << FORMAT("{} {} {} / {}: {} - {}\n", std::ctime(&time), dummy->name(), i++, steps+1, value, dyson_factor[value]) << std::flush;
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/D1600_{}.txt", ff_str));
    }

#endif
    #endif
}