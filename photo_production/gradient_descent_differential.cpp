#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"


#include <feynumeric/feynumeric.hpp>
#include <feynumeric/utility.hpp>

#include <fmt/chrono.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <feynumeric/phase_space.hpp>

double sigmoid(double x){ return 1./(1. + std::exp(-x));}


Feynumeric::func_t<1> interval(double min, double max){
    return [min, max](double x){ return min + (max-min) * sigmoid(x);};
}
Feynumeric::func_t<1> inverse_interval(double min, double max){
    return [min, max](double x){ return std::log( (max-x)/(x-min) );};
}

int main(int argc, char** argv){
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    std::string const CMD_DATA_FILE_PIP = "data_file_pip";
    std::string const CMD_DATA_FILE_PIM = "data_file_pim";
    std::string const CMD_DATA_FILE_PIM_CE = "data_file_pim_ce";

    Command_Line_Manager cmd(argc, argv);
    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("coupling_constants", true, "file with coupling constants");
    cmd.register_command("n_epochs", std::string("100"), "number of epochs to train");
    cmd.register_command("rate", std::string("0.001"), "rate of descent");
    cmd.register_command(CMD_DATA_FILE_PIP, true, "file with fit data");
    cmd.register_command(CMD_DATA_FILE_PIM, true, "file with fit data");
    cmd.register_command(CMD_DATA_FILE_PIM_CE, true, "file with fit data");
    cmd.register_command("k", 1.0, "pdg-deviation rate");
    cmd.register_command("fit_params", true, "file with initial parameters and particles");
    cmd.register_command("energies", true, "comma separated energies (no spaces), e.g. 1.2,1.25,1.3 in GeV");
    cmd.register_command("channel", std::string("s"), "what channels to use, args: stuc");
    cmd.crash_on_missing_mandatory_command();

    auto const &channel = cmd.as_string("channel");
    bool const s_channel_enabled = channel.find('s') != std::string::npos;
    bool const t_channel_enabled = channel.find('t') != std::string::npos;
    bool const u_channel_enabled = channel.find('u') != std::string::npos;
    bool const c_channel_enabled = channel.find('c') != std::string::npos;

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Neutron = P.get("neutron");
    auto const& Pi_Plus = P.get("pi+");
    auto const& Pi_Minus = P.get("pi-");
    auto const& Pi_Zero = P.get("pi0");
    init_vertices(P, cmd.as_string("coupling_constants"));

    std::vector<Polynomial> s_poly;
    std::vector<Polynomial> c_poly;

    struct Particle_Fit{
        std::string name, fit_type;
        func_t<1> f;
        func_t<1> inverse_f;
        double min, max, fit_value;
        Particle_Fit(std::string const& str, std::string const& fit_type, double min, double max, std::mt19937& gen)
                : name(str)
                , fit_type(fit_type)
                , min(min)
                , max(max)
        {
            if( fit_type == "bounded"){
                f = interval(min, max);
                inverse_f = inverse_interval(min, max);
                fit_value = std::uniform_real_distribution<double>(min, max)(gen);
            }else if( fit_type == "unbounded" ){
                f = [](double x){return x;};
                inverse_f = [](double x){return x;};
                fit_value = std::uniform_real_distribution<double>(min, max)(gen);
            }else if( fit_type == "constant" ){
                f = [min](double){return min;};
                inverse_f = [min](double){return min;};
                fit_value = min;
            }else{
                critical_error(FORMAT("Unknown fit type: {}.\n", fit_type));
            }
        }
    };

    std::random_device r;
    std::mt19937 random_generator{r()};


    std::ifstream infile(cmd.as_string("fit_params"));
    if( !infile ){
        std::cerr << "Could not open fit parameter file.\n" << cmd.as_string("fit_params") << "\n";
    }
    std::string particle_name;
    std::string fit_type;
    double min, max;
    std::vector<Particle_Fit> particle_fit;
    while (infile >> particle_name >> fit_type >> min >> max)
    {
        particle_fit.emplace_back(particle_name, fit_type, min, max, random_generator);
    }


    std::uniform_real_distribution<double> dist(0, 1);

    auto const N_EPOCHS = cmd.as_int("n_epochs");
    auto rate  = cmd.as_double("rate");

    std::cout << FORMAT("n_epochs: {} rate: {}\n", N_EPOCHS, rate);

    for( auto const& pf : particle_fit ){
        std::cout << pf.name << " " << pf.fit_type << " " << pf.fit_value << "\n";
    }
    
    std::vector<Amplitude<0>> width_amplitudes;

    Timer stopwatch;
    stopwatch.start();

    for( auto const& pf : particle_fit ){
        Particle_Ptr particle;
        if( pf.name.starts_with('D') ){
            particle = P.get(FORMAT("{}pp", pf.name));
        }else{
            particle = P.get(pf.name);
        }

        auto n_spin_states = particle->spin().n_states() * Proton->spin().n_states();
        std::vector<std::vector<Polynomial>> width_polynomials(1);
        for( std::size_t i = 0; i < n_spin_states; ++i ){
            Polynomial temp;
            if( particle->name().starts_with('D') ){
                temp.load(FORMAT("data/polynomials/polynomial_decay_{}_0_{}.txt", particle->name(), i));
            }
            width_polynomials[0].push_back(temp);
        }

        Amplitude<0> M_width(width_polynomials, {particle}, {Proton, Pi_Plus});
        width_amplitudes.push_back(M_width);
    }

    struct Data{
        double angle, obs, obsx, err, dchi;
    };

    std::map<double, std::vector<Data>> data_points_pip;
    std::map<double, std::vector<Data>> data_points_pim;
    std::map<double, std::vector<Data>> data_points_pimCE;

    std::vector<std::map<double, std::vector<Data>>*> data_pointers = {&data_points_pip, &data_points_pim, &data_points_pimCE};
    std::string buffer;
    buffer.reserve(100UL);

    int i = 0;
    for( auto const& cmd_str : {CMD_DATA_FILE_PIP, CMD_DATA_FILE_PIM, CMD_DATA_FILE_PIM_CE} ){
        std::string file = cmd.as_string(cmd_str);
        std::ifstream in(file);
        if( !in ){
            critical_error(FORMAT("Could not open: {}", file));
        }
        while( std::getline(in, buffer) ){
            if( buffer.starts_with('#') ){
                continue;
            }
            std::stringstream stream(buffer);
            Data row;
            double tlab;
            stream >> tlab;
            stream >> row.angle >> row.obs >> row.obsx >> row.err >> row.dchi;
            double converted = std::sqrt((tlab/1000. + Pi_Plus->mass()) * 2. * Proton->mass() + Proton->mass() * Proton->mass() + Pi_Plus->mass() * Pi_Plus->mass());
            (*data_pointers[i])[converted].push_back(row);
        }
        i++;
    }

    std::map<double, std::vector<Data>> sample_points_pip;
    std::map<double, std::vector<Data>> sample_points_pim;
    std::map<double, std::vector<Data>> sample_points_pimCE;

    std::istringstream tokenizer(cmd.as_string("energies"));
    while( std::getline(tokenizer, buffer, ',') ){
        double sqrt_s = std::stod(buffer);
        auto ptr1_prev = data_points_pip.upper_bound(sqrt_s);
        auto ptr1_next = ptr1_prev;
        auto ptr2_prev = data_points_pim.upper_bound(sqrt_s);
        auto ptr2_next = ptr2_prev;
        auto ptr3_prev = data_points_pimCE.upper_bound(sqrt_s);
        auto ptr3_next = ptr3_prev;

        auto const threshold = 10;

        while( ptr1_prev->second.size() < threshold && ptr1_next->second.size() < threshold ){
            ptr1_prev--;
            ptr1_next++;
        }

        while( ptr2_prev->second.size() < threshold && ptr2_next->second.size() < threshold ){
            ptr2_prev--;
            ptr2_next++;
        }

        while( ptr3_prev->second.size() < threshold && ptr3_next->second.size() < threshold ){
            ptr3_prev--;
            ptr3_next++;
        }

        auto ptr1 = ptr1_prev->second.size() > ptr1_next->second.size()? ptr1_prev : ptr1_next;
        auto ptr2 = ptr2_prev->second.size() > ptr2_next->second.size()? ptr2_prev : ptr2_next;
        auto ptr3 = ptr3_prev->second.size() > ptr3_next->second.size()? ptr3_prev : ptr3_next;


        sample_points_pip[sqrt_s] = ptr1->second;
        sample_points_pim[sqrt_s] = ptr2->second;
        sample_points_pimCE[sqrt_s] = ptr3->second;

        auto to_tlab = [](double s, double mN, double mpi){ return 1000*((s*s-mN*mN-mpi*mpi)/(2*mN)-mpi);};
        auto to_tlab_bound = [&](double s){ return to_tlab(s, Proton->mass(), Pi_Plus->mass());};
        std::cout << FORMAT("sqrt_s: {}/{} {}/{} ({}) {}/{} ({}) {}/{} ({})\n",
                            sqrt_s, to_tlab_bound(sqrt_s),
                            ptr1->first, to_tlab_bound(ptr1->first), ptr1->second.size(),
                            ptr2->first, to_tlab_bound(ptr2->first), ptr2->second.size(),
                            ptr3->first, to_tlab_bound(ptr3->first), ptr3->second.size());
    }

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

    std::map<std::string, double> initial_values;

    for( auto const &resonance : particle_fit )
    {
        if( P.exists(resonance.name) && !P.get(resonance.name)->is_group() )
        {
            initial_values[resonance.name] = couplings.get(coupling_string("N", "N", resonance.name));
        }else{
            initial_values[resonance.name] = couplings.get(coupling_string("N", "Pion", resonance.name));
        }
        std::cout << "g: " << initial_values[resonance.name] << "\n";

        std::vector<std::string> strs = {pp_string(resonance.name), p_string(resonance.name), n_string(resonance.name), m_string(resonance.name)};

        for( auto const &str: strs )
        {
            if( P.exists(str))
            {
                P[str]->user_data("form_factor", Form_Factor::ff_dict[Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER]);
                Polynomial poly;
                poly.load(FORMAT("widths/{}_width_N_Pi.poly", resonance.name));
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
        if( t_channel_enabled && P.exists(resonance.name) ){
            auto particle = P[resonance.name];
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






    stopwatch.stop();
    std::cout << "config done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
//
//    stopwatch.start();
//
//    double const epsilon = 0.001;
//
//    Timer epoch_timer;
//
//    double last_loss = std::numeric_limits<double>::max();
//
//    for( std::size_t epoch = 0; epoch < N_EPOCHS; ++epoch ){
//        epoch_timer.start();
//        double loss = 0.;
//
//        for( std::size_t i = 0; i < particle_fit.size(); ++i ){
//
//        }
//
//        for( std::size_t k = 0; k < N_amplitudes; ++k ){
//            double temp = rate * derivatives[k] / sample_points.size();
//            particle_fit[k].fit_value -= temp;
//        }
//
//        epoch_timer.stop();
//        std::cout << "Epoch: " << epoch << " [" << epoch_timer.time<std::chrono::milliseconds>()/1000. << "] rate: [" << rate << "] loss: [" << loss / sample_points.size() << "]\n";
//
//        if( loss > last_loss ){
//            rate *= 0.9;
//        }
//
//        for( std::size_t i = 0; i < N_amplitudes; ++i ){
//            std::cout << FORMAT("{}: {} {}\n", particle_fit[i].name, particle_fit[i].fit_value, particle_fit[i].f(particle_fit[i].fit_value).real());
//        }
//
//        if( last_loss > loss && (last_loss - loss) / sample_points.size() < 1.e-2 ){
//            std::cout << "last_loss: " << last_loss << "\n";
//            std::cout << "loss: " << loss << "\n";
//            std::cout << "ss: " << sample_points.size() << "\n";
//            break;
//        }
//        last_loss = loss;
//    }
//
////    std::ofstream out("plot.txt");
//
//
//    scattering_amplitudes[0].scale([=](double){return particle_fit[0].f(particle_fit[0].fit_value).real();});
//    Amplitude<1> interference = scattering_amplitudes[0];
//
//    for( std::size_t i = 1; i < N_amplitudes; ++i ){
//        std::cout << particle_fit[i].name << ": " << particle_fit[i].f(particle_fit[i].fit_value).real() << "\n";
//        scattering_amplitudes[i].scale([=](double){return particle_fit[i].f(particle_fit[i].fit_value).real();});
//        interference = interference + scattering_amplitudes[i];
//    }
//
////    interference = interference + scattering_amplitudes[N_amplitudes-1];
//    auto result = interference.scattering();
//
//    std::cout << "interference={";
//    for( auto& point : data_points ){
//        std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
//    }
//    std::cout << "\b};";
//
//    for( std::size_t i = 0; i < N_amplitudes; ++i ){
//        std::cout << "\nresonance" << particle_fit[i].name << "={";
//        result = scattering_amplitudes[i].scattering();
//        for( auto& point : data_points ){
//            std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
//        }
//        std::cout << "\b};";
//    }
//
////    result = scattering_amplitudes.back().scattering();
////    std::cout << "\ntchannel={";
////    for( auto& point : data_points ){
////        std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
////    }
////    std::cout << "\b};\n";
//
//    stopwatch.stop();
//    std::cout << "Minimization done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
}