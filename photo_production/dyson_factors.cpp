#include <feynumeric/feynumeric.hpp>
#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include <omp.h>

#define D1232 0
#define N1440 1
#define N1520 0
#define N1535 0

int main(int argc, char** argv)
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    Command_Line_Manager cmd(argc, argv);
    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("start", std::string("1.0"), "starting point");
    cmd.register_command("end", std::string("3.0"), "end value");
    cmd.register_command("steps", std::string("200"), "steps");

    cmd.crash_on_missing_mandatory_command();

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Neutron = P.get("neutron");
    auto const& Pi_Zero = P.get("pi0");
    auto const& Pi_Minus = P.get("pi-");
    auto const& Pi_Plus = P.get("pi+");

    for( auto& [key, particle] : P ){
        auto k = key;
        auto p = particle;
        particle->user_data("form_factor", identity);
    }

    init_vertices(P);

    double const start = cmd.as_double("start");
    double const end   = cmd.as_double("end");
    int    const steps = cmd.as_int("steps");
    std::cout << FORMAT("start: {} end: {} steps: {}\n", start, end, steps);
    #if D1232
    {/// D1232 -> Npi
        auto dummypp = std::make_shared<Particle>(*P.get("D1232pp"));
        dummypp->user_data("form_factor", identity);
        auto decay_1 = create_diagram(FORMAT("decay {} to proton pi+", dummypp->name()), Decay_1_to_2, VMP,
                                      {dummypp},
                                      {},
                                      {Proton, Pi_Plus}
        );
        Feynman_Process decay({decay_1});
        double Gamma0 = decay.decay_width();
        std::cout << "Gamma0: " << Gamma0 << "\n";
        std::map<double, double> dyson_factor;
        for( int i = 0; i <= steps; ++i ){
            double value = start + (end - start) / steps * i;
            dummypp->mass(value);
            auto temp = decay.decay_width();// /Gamma0;
            dyson_factor[value] = std::isnan(temp)? 0 : temp;
        }
		Table(dyson_factor).write("data/dyson_factors/D1232.txt");
    }
    #endif
    #if N1440
    {   /// N1440 -> Npi and N1440 -> Npipi
        std::cout << FORMAT("g: {}\n", P.get("N1440p")->user_data<double>("gD1232N1440pi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232pp")->user_data<double>("gD1232N1440pi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232p")->user_data<double>("gD1232N1440pi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232n")->user_data<double>("gD1232N1440pi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232m")->user_data<double>("gD1232N1440pi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232pp")->user_data<double>("gRNpi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232p")->user_data<double>("gRNpi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232n")->user_data<double>("gRNpi"));
        std::cout << FORMAT("g: {}\n", P.get("D1232m")->user_data<double>("gRNpi"));
        auto dummy = std::make_shared<Particle>(*P.get("N1440p"));
        dummy->user_data("form_factor", identity);
        P.get("D1232pp")->user_data("form_factor", identity);
        P.get("D1232p")->user_data("form_factor", identity);
        P.get("D1232n")->user_data("form_factor", identity);
        P.get("D1232m")->user_data("form_factor", identity);
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

        Feynman_Process decay_delta1({diagram_pipi_1, diagram_pipi_2});
        Feynman_Process decay_delta2({diagram_pipi_4, diagram_pipi_5});
        Feynman_Process decay_delta3({diagram_pipi_7, diagram_pipi_8});

        Feynman_Process decay_f5001({diagram_pipi_3});
        Feynman_Process decay_f5002({diagram_pipi_6});

        std::map<double, double> dyson_factor;
        for( int i = 0; i <= steps; ++i ){
            double value = start + (end - start) / steps * i;
            dummy->mass(value);
            auto temp1 = decay_pi1.decay_width();
            auto temp2 = decay_pi2.decay_width();
            temp1 = std::isnan(temp1) ? 0  : temp1;
            temp2 = std::isnan(temp2) ? 0  : temp2;
            auto delta1 = decay_delta1.decay_width();
            auto delta2 = decay_delta2.decay_width();
            auto delta3 = decay_delta3.decay_width();
//            auto temp3 = decay_pipi1.decay_width();
//            auto temp4 = decay_pipi2.decay_width();
//            auto temp5 = decay_pipi3.decay_width();
//            temp3 = std::isnan(temp3) ? 0  : temp3;
//            temp4 = std::isnan(temp4) ? 0  : temp4;
//            temp5 = std::isnan(temp5) ? 0  : temp5;
//            std::cout << FORMAT("mass: {}\ntemp1: {} temp2: {}\ntemp3: {} temp4: {} temp5: {}\n", value, temp1, temp2, temp3, temp4, temp5);
            std::cout << FORMAT("Npi: {}, d1232: {} f500: {}\n", temp1+temp2, delta1+delta2+delta3, decay_f5001.decay_width() + decay_f5002.decay_width());
            std::cout << FORMAT("w1: {}, w2: {}, w3: {}\n", delta1, delta2, delta3);
//            dyson_factor[value] = temp1 + temp2 + temp3 + temp4 + temp5;
            if( i%20 == 0 ){
                std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                std::cout << FORMAT("{} N1440 {} / {}\n", std::ctime(&time), i, steps) << std::flush;
            }
        }
        Table(dyson_factor).write("data/dyson_factors/N1440.txt");
    }
    #endif
    #if N1520
	{/// N1520 -> Npi and N1520 -> Npipi
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

//        double Gamma0_1 = decay_pi.decay_width();
//        double Gamma0_2 = decay_pipi.decay_width();
        std::map<double, double> dyson_factor;

        #pragma omp parallel for
        for( int i = 0; i <= steps; ++i ){
            double value = start + (end - start) / steps * i;
            dummy->mass(value);
            auto temp1 = copies_pi[i].decay_width();// / Gamma0_1;
            temp1 = std::isnan(temp1) ? 0  : temp1;
            auto temp2 = copies_pipi[i].decay_width();// / Gamma0_2;
            temp2 = std::isnan(temp2) ? 0  : temp2;
            dyson_factor[value] = temp1 + temp2;
            if( i%20 == 0 ){
                std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                std::cout << FORMAT("{} N1520 {} / {} core: {}\n", std::ctime(&time), i, steps, omp_get_thread_num()) << std::flush;
            }
        }
		Table(dyson_factor).write("data/dyson_factors/N1520.txt");
	}
    #endif
    #if N1535
    {/// N1535 -> Npi and N1535 -> Neta
        auto dummy = std::make_shared<Particle>(*P.get("N1535p"));
        dummy->user_data("form_factor", identity);
        auto diagram_pi = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_2, VMP,
                                         {dummy},
                                         {},
                                         {Proton, Pi_Zero}
        );
        auto diagram_eta = create_diagram(FORMAT("decay {} to proton eta", dummy->name()), Decay_1_to_2, VMP,
                                         {dummy},
                                         {},
                                         {Proton, Pi_Zero}
        );


        Feynman_Process decay_pi({diagram_pi});
        Feynman_Process decay_eta({diagram_eta});

        std::vector<Feynman_Process> copies_pi;
        std::vector<Feynman_Process> copies_eta;
        copies_pi.reserve(steps);
        copies_eta.reserve(steps);
        for( std::size_t i = 0; i <= steps; ++i ){
            copies_pi.emplace_back(decay_pi);
            copies_eta.emplace_back(decay_eta);
        }

//        double Gamma0_1 = decay_pi.decay_width();
//        double Gamma0_2 = decay_eta.decay_width();
        std::map<double, double> dyson_factor;
        #pragma omp parallel for
        for( int i = 0; i <= steps; ++i ){
            double value = start + (end - start) / steps * i;
            dummy->mass(value);
            auto temp1 = copies_pi[i].decay_width();// / Gamma0_1;
            temp1 = std::isnan(temp1) ? 0  : temp1;
            auto temp2 = copies_eta[i].decay_width();// / Gamma0_2;
            temp2 = std::isnan(temp2) ? 0  : temp2;
            dyson_factor[value] = temp1 + temp2;
            if( i%20 == 0 ){
                std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                std::cout << FORMAT("{} N1535 {} / {}\n", std::ctime(&time), i, steps) << std::flush;
            }
        }
        Table(dyson_factor).write("data/dyson_factors/N1535.txt");
    }
    #endif
}