#include <feynumeric/feynumeric.hpp>
#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include <omp.h>


int main(int argc, char** argv)
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using namespace std::placeholders;

    Command_Line_Manager cmd(argc, argv);
    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("coupling_constants", true, "file with coupling constants");
    cmd.register_command("start", std::string("1.0"), "starting point");
    cmd.register_command("end", std::string("3.0"), "end value");
    cmd.register_command("steps", std::string("200"), "steps");
    cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING, Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ, Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));

    cmd.crash_on_missing_mandatory_command();

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Neutron = P.get("neutron");
    auto const& Pi_Zero = P.get("pi0");
    auto const& Pi_Minus = P.get("pi-");
    auto const& Pi_Plus = P.get("pi+");

    auto breit_wigner_modified = [](Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, double E, double M){
        double const l = .8;
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
    if( Form_Factor::ff_dict.contains(ff_str) ){
        ff = Form_Factor::ff_dict[ff_str];
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

    init_vertices(P, cmd.as_string("coupling_constants"));

    double const start = cmd.as_double("start");
    double const end   = cmd.as_double("end");
    int    const steps = cmd.as_int("steps");

    std::vector<double> values;
    values.resize(steps+2);
    for( int i = 0; i <= steps; ++i ){
        values[i] = (start + (end - start) / steps * i);
    }
//    values = {1.232};

	bool D1232 = cmd.is_enabled("D1232");
    bool D1600 = cmd.is_enabled("D1600");
    bool D1620 = cmd.is_enabled("D1620");
    bool D1700 = cmd.is_enabled("D1700");
    bool D1750 = cmd.is_enabled("D1750");
    bool D1900 = cmd.is_enabled("D1900");
    bool D1905 = cmd.is_enabled("D1905");
    bool D1910 = cmd.is_enabled("D1910");
    bool D1920 = cmd.is_enabled("D1920");
    bool D1930 = cmd.is_enabled("D1930");
    bool D1940 = cmd.is_enabled("D1940");
    bool D1950 = cmd.is_enabled("D1950");


    std::cout << FORMAT("start: {} end: {} steps: {}\n", start, end, steps);
    if( D1232 )
    {
        std::string particle_name = "D1232";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
#pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            ++completed;
            auto value = values[i];
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, P.get(FORMAT("{}pp", particle_name))->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
            //dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    if( D1600 )
    {
        std::string particle_name = "D1600";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 1;
        #pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            ++completed;
            auto value = values[i];
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, P.get(FORMAT("{}pp", particle_name))->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    /*
    if( D1600 )
    {
        std::cout << "D1600\n";
        auto dummy = std::make_shared<Particle>(*P.get("D1600pp"));
        dummy->user_data("form_factor", Form_Factor::identity);
        /** D1600 -> Npi **//*
        auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                         {dummy},
                                         {},
                                         {Proton, Pi_Plus}
        );
        /** D1600 -> Npipi **/
        /*
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
        */
/*

        std::vector<Feynman_Process> processes;
        processes.push_back(Feynman_Process({diagram_pi1}));
//        processes.push_back(Feynman_Process({diagram_pip_pin_1, diagram_pip_pin_2, diagram_pip_pin_3}));
//        processes.push_back(Feynman_Process({diagram_pip_pip_1, diagram_pip_pip_2, diagram_pip_pip_3, diagram_pip_pip_4}));

        std::map<double, double> dyson_factor;
        values.back() = dummy->mass();
        std::size_t i{0};
        for( auto const& value : values ){
            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width();
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << FORMAT("{} {} {} / {}: {} - {}\n", std::ctime(&time), dummy->name(), i++, steps+1, value, dyson_factor[value]) << std::flush;
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/D1600_{}.txt", ff_str));
    }
    */
    if( D1620 )
    {
        std::string particle_name = "D1620";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
#pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, P.get(FORMAT("{}pp", particle_name))->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    /*
    if( D1620 )
    {
        std::cout << "D1620\n";
        auto dummy = std::make_shared<Particle>(*P.get("D1620pp"));
        dummy->user_data("form_factor", Form_Factor::identity);
        */
        /** D1620 -> Npi **/
        /*
        auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                          {dummy},
                                          {},
                                          {Proton, Pi_Plus}
        );
        /** D1620 -> Npipi **/
        /*
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

        auto diagram_pip_pin_3 = create_diagram(FORMAT("decay {} -> Rho -> p pi+pi0 4", dummy->name()), Decay_1_to_1_M2, VMP,
                                                {dummy},
                                                {P.get("rho+")},
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

        */
        /*
        std::vector<Feynman_Process> processes;
        processes.push_back(Feynman_Process({diagram_pi1}));
//        processes.push_back(Feynman_Process({diagram_pip_pin_1, diagram_pip_pin_2, diagram_pip_pin_3}));
//        processes.push_back(Feynman_Process({diagram_pip_pip_1, diagram_pip_pip_2}));

        std::map<double, double> dyson_factor;
        values.back() = dummy->mass();
        std::size_t i{0};
        for( auto const& value : values ){
            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width();
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << FORMAT("{} {} {} / {}: {} - {}\n", std::ctime(&time), dummy->name(), i++, steps+1, value, dyson_factor[value]) << std::flush;
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/D1620_{}.txt", ff_str));
    }
         */
    if( D1700 )
    {
        std::string particle_name = "D1700";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
        #pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, P.get(FORMAT("{}pp", particle_name))->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    /*
	if( D1700 )
	{
		std::cout << "D1700\n";
		auto dummy = std::make_shared<Particle>(*P.get("D1700pp"));
		dummy->user_data("form_factor", Form_Factor::identity);
		/** D1700 -> Npi **/
    /*
		auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
		                                  {dummy},
		                                  {},
		                                  {Proton, Pi_Plus}
		);
		/** D1700 -> Npipi **/
		/*
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

		auto diagram_pip_pin_3 = create_diagram(FORMAT("decay {} -> Rho -> p pi+pi0 4", dummy->name()), Decay_1_to_1_M2, VMP,
		                                        {dummy},
		                                        {P.get("rho+")},
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

		*/
/*

		std::vector<Feynman_Process> processes;
		processes.push_back(Feynman_Process({diagram_pi1}));
//		processes.push_back(Feynman_Process({diagram_pip_pin_1, diagram_pip_pin_2, diagram_pip_pin_3}));
//		processes.push_back(Feynman_Process({diagram_pip_pip_1, diagram_pip_pip_2}));

		std::map<double, double> dyson_factor;
		values.back() = dummy->mass();
		std::size_t i{0};
		for( auto const& value : values ){
			dummy->mass(value);
            dyson_factor[value] = 0;
			for( auto& process : processes ){
				auto temp = process.decay_width();
				dyson_factor[value] += std::isnan(temp)? 0 : temp;
			}
			std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
			std::cout << FORMAT("{} {} {} / {}: {} - {}\n", std::ctime(&time), dummy->name(), i++, steps+1, value, dyson_factor[value]) << std::flush;
		}
		Table(dyson_factor).write(FORMAT("data/dyson_factors/D1700_{}.txt", ff_str));
	}
    */
    if( D1900 )
    {
        std::string particle_name = "D1900";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
        #pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, P.get(FORMAT("{}pp", particle_name))->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    if( D1905 )
    {
        std::string particle_name = "D1905";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
        #pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, P.get(FORMAT("{}pp", particle_name))->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    if( D1910 )
    {
        std::string particle_name = "D1910";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
        #pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, P.get(FORMAT("{}pp", particle_name))->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    if( D1920 )
    {
        std::string particle_name = "D1920";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
#pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, particle->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    if( D1930 )
    {
        std::string particle_name = "D1930";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
#pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, particle->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    if( D1940 )
    {
        std::string particle_name = "D1940";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
#pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, particle->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
    if( D1950 )
    {
        std::string particle_name = "D1950";
        auto particle = P.get(FORMAT("{}pp", particle_name));
        std::cout << particle_name << "\n";
        std::map<double, double> dyson_factor;
        values.back() = particle->mass();
        std::size_t completed = 0;
#pragma omp parallel for
        for( std::size_t i = 0; i < values.size(); ++i ){
            auto value = values[i];
            ++completed;
            auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
            if( ff_str == Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER ){
                FORM_FACTOR_FUNCTION breit_wigner_bind = std::bind(breit_wigner_modified, _1, _2, _3, _4, particle->mass());
                dummy->user_data("form_factor", breit_wigner_bind);
            }else{
                dummy->user_data("form_factor", Form_Factor::identity);
            }
            /**  -> Npi **/
            auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                              {dummy},
                                              {},
                                              {Proton, Pi_Plus}
            );
            std::vector<Feynman_Process> processes;
            processes.push_back(Feynman_Process({diagram_pi1}));
//            dummy->mass(value);
            dyson_factor[value] = 0;
            for( auto& process : processes ){
                auto temp = process.decay_width(value);
                dyson_factor[value] += std::isnan(temp)? 0 : temp;
            }
            if( (completed / 10)%10 == 0 ){
                std::cout << "#" << std::flush;
            }
        }
        Table(dyson_factor).write(FORMAT("data/dyson_factors/{}_{}.txt", particle_name, ff_str));
    }
   
}