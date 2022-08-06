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

double inverse_sigmoid(double x){
    return std::log(x/(1.-x));
}

Feynumeric::func_t<1> interval(double min, double max){
    return [min, max](double x){ return min + (max-min) * sigmoid(x);};
}
Feynumeric::func_t<1> inverse_interval(double min, double max){
    return [min, max](double x){ return std::log( (max-x)/(x-min) );};
}

double interval_isigmoid(double min, double max, double x){
    return std::log( (max-x)/(x-min) );
}

double slope(double a, double b, double c){
    return std::max(0., std::abs(c - (b+a)/2.) - (b-a)/2.);
}

double out_of_interval_cost(double a, double b, double c, double k){
    auto arg = slope(a, b, c);
    return std::exp(k * arg * (c-b)) + std::exp(-k * arg * (c-a)) - 2 * std::exp(k/2. * (a-b) * arg);
}

int main(int argc, char** argv){
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    Command_Line_Manager cmd(argc, argv);
    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("coupling_constants", true, "file with coupling constants");
    cmd.register_command("n_epochs", std::string("100"), "number of epochs to train");
    cmd.register_command("rate", std::string("0.001"), "rate of descent");
    cmd.register_command("data_file", true, "file with fit data");
    cmd.register_command("k", 1.0, "pdg-deviation rate");
    cmd.register_command("fit_params", true, "file with initial parameters and particles");
//	cmd.register_command("start", std::string("1.0"), "starting point");
//	cmd.register_command("end", std::string("3.0"), "end value");
//	cmd.register_command("steps", std::string("200"), "steps");
//	cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING, Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ, Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));
    cmd.crash_on_missing_mandatory_command();

    double const k_factor = cmd.as_double("k");

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Pi_Plus = P.get("pi+");
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
//	std::default_random_engine eng{r()};
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
    std::vector<Amplitude<1>> scattering_amplitudes;
    std::vector<std::pair<Amplitude<1>,Amplitude<1>>> scattering_amplitudes_derivatives;

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

        /* scattering */
        n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
        std::vector<std::vector<Polynomial>> scattering_polynomials(2);
        for( std::size_t i = 0; i < n_spin_states; ++i ){
            Polynomial temp_s, temp_c;
            temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
            temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
            scattering_polynomials[0].push_back(temp_s);
            scattering_polynomials[1].push_back(temp_c);
        }

        func_t<1> breit_wigner = [=](double sqrt_s) mutable -> Complex {
            if( particle->is_fermion() ){
                return 1. / (sqrt_s * sqrt_s - particle->mass() * particle->mass() +
                             1.i * sqrt_s * M_width.width(sqrt_s));
            }
            else{
                return 1. / (sqrt_s * sqrt_s - particle->mass() * particle->mass());
            }
        };

        Amplitude<1> M_scattering(scattering_polynomials, {Proton, Pi_Plus}, {Proton, Pi_Plus});

        func_t<1> form_factor = [=](double sqrt_s) mutable{
            auto const lambda = 0.8;
            auto const l4 = std::pow(lambda, 4);
            auto const delta = sqrt_s * sqrt_s - particle->mass() * particle->mass();
            return l4 / ( delta * delta + l4 );
        };

        M_scattering.scale(breit_wigner);
        if( particle->is_fermion() ){
            M_scattering.scale(form_factor);
        }
        scattering_amplitudes.push_back(M_scattering);

    }

    if( scattering_amplitudes.empty() ){
        error("No particle selected.");
        return EXIT_SUCCESS;
    }
    /*
    {
        auto particle = P.get("f0_500");

        auto n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
        std::vector<std::vector<Polynomial>> scattering_polynomials(2);
        for( std::size_t i = 0; i < n_spin_states; ++i )
        {
            Polynomial temp_s, temp_c;
            temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
            temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
            scattering_polynomials[0].push_back(temp_s);
            scattering_polynomials[1].push_back(temp_c);
        }
        Amplitude<1> M_scattering(scattering_polynomials, {Proton, Pi_Plus}, {Proton, Pi_Plus});

        scattering_amplitudes.push_back(M_scattering);
    }

    {
        auto particle = P.get("rho0");

        auto n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
        std::vector<std::vector<Polynomial>> scattering_polynomials(2);
        for( std::size_t i = 0; i < n_spin_states; ++i )
        {
            Polynomial temp_s, temp_c;
            temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
            temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
            scattering_polynomials[0].push_back(temp_s);
            scattering_polynomials[1].push_back(temp_c);
        }
        Amplitude<1> M_scattering(scattering_polynomials, {Proton, Pi_Plus}, {Proton, Pi_Plus});

        scattering_amplitudes.push_back(M_scattering);
    }
    */



    stopwatch.stop();
    std::cout << "config done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
    stopwatch.start();

    std::string file = cmd.as_string("data_file");

    std::ifstream in(file);

    if( !in ){
        critical_error(FORMAT("Could not open: {}", file));
    }

    std::string buffer;
    buffer.reserve(100UL);

    struct Data{
        double srt, cross, dcros, plab;
    };

    std::vector<Data> data_points;

    while( std::getline(in, buffer) ){
        if( buffer.starts_with('#') ){
            continue;
        }
        std::stringstream stream(buffer);
        Data row;
        stream >> row.srt >> row.cross >> row.dcros >> row.plab;
        data_points.push_back(row);
    }

    std::vector<Data> sample_points(data_points.begin(), data_points.begin()+260);

    std::size_t const N_amplitudes = scattering_amplitudes.size();

    double const epsilon = 0.001;

    Timer epoch_timer;

    double last_loss = std::numeric_limits<double>::max();

    for( std::size_t epoch = 0; epoch < N_EPOCHS; ++epoch ){
        epoch_timer.start();
        double loss = 0.;

        std::vector<Amplitude<1>> scattering_left(scattering_amplitudes.begin(), scattering_amplitudes.end());
        std::vector<Amplitude<1>> scattering_mid(scattering_amplitudes.begin(), scattering_amplitudes.end());
        std::vector<Amplitude<1>> scattering_right(scattering_amplitudes.begin(), scattering_amplitudes.end());

        for( std::size_t i = 0; i < particle_fit.size(); ++i ){
            auto const& p = particle_fit[i];
            auto scale_left = p.f(p.fit_value - epsilon);
            auto scale_mid = p.f(p.fit_value);
            auto scale_right = p.f(p.fit_value + epsilon);
            scattering_left[i].scale([scale_left](double){ return scale_left; });
            scattering_mid[i].scale([scale_mid ](double){ return scale_mid;});
            scattering_right[i].scale([scale_right](double){ return scale_right; });
        }

        func_t<1> y_hat;
        std::vector<func_t<1>> y_hat0(N_amplitudes);
        std::vector<func_t<1>> y_hat1(N_amplitudes);


        Amplitude<1> interference = scattering_mid[0];
        for( std::size_t i = 1; i < N_amplitudes; ++i ){
            interference = interference + scattering_mid[i];
            Amplitude<1> interference_left = scattering_left[0];
            Amplitude<1> interference_right = scattering_right[0];
            for( std::size_t j = 1; j < N_amplitudes; ++j ){
                if( i == j ){
                    interference_left = interference_left + scattering_left[j];
                    interference_right = interference_right + scattering_right[j];
                } else{
                    interference_left = interference_left + scattering_mid[j];
                    interference_right = interference_right + scattering_mid[j];
                }
            }
//            interference_left = interference_left + scattering_mid[N_amplitudes - 1];
//            interference_right = interference_right + scattering_mid[N_amplitudes - 1];
            y_hat0[i] = interference_left.scattering();
            y_hat1[i] = interference_right.scattering();
        }
//        interference = interference + scattering_mid[N_amplitudes - 1];
        y_hat = interference.scattering();


        std::vector<double> derivatives(N_amplitudes, 0.);


        for( auto const& point : sample_points ){
            double diff = point.cross - y_hat(point.srt).real() * static_cast<double>(1._2mbarn);
            loss += diff*diff / (point.dcros * point.dcros);
            for( std::size_t k = 1; k < N_amplitudes; ++k ){
                double const loss0 = std::pow(point.cross - y_hat0[k](point.srt).real() * static_cast<double>(1._2mbarn), 2) / (point.dcros * point.dcros);
                double const loss1 = std::pow(point.cross - y_hat1[k](point.srt).real() * static_cast<double>(1._2mbarn), 2) / (point.dcros * point.dcros);
                double const d = (loss1 - loss0) / (2. * epsilon);
                derivatives[k] += d;
            }
        }

//        for( std::size_t k = 0; k < particle_fit.size(); ++k ){
//            auto const& resonance = particle_fit[k];
//            loss += out_of_interval_cost(resonance.min, resonance.max, sigmoid(sigm_branching_ratios[k]), k_factor);
//            auto y0 = out_of_interval_cost(resonance.min, resonance.max, sigmoid(sigm_branching_ratios[k] - epsilon), k_factor);
//            auto y1 = out_of_interval_cost(resonance.min, resonance.max, sigmoid(sigm_branching_ratios[k] + epsilon), k_factor);
//            derivatives[k] += (y1-y0)/(2. * epsilon);
//        }

        for( std::size_t k = 0; k < N_amplitudes; ++k ){
            double temp = rate * derivatives[k] / sample_points.size();
            particle_fit[k].fit_value -= temp;
        }

        epoch_timer.stop();
        std::cout << "Epoch: " << epoch << " [" << epoch_timer.time<std::chrono::milliseconds>()/1000. << "] rate: [" << rate << "] loss: [" << loss / sample_points.size() << "]\n";

        if( loss > last_loss ){
            rate *= 0.9;
        }

        for( std::size_t i = 0; i < N_amplitudes; ++i ){
            std::cout << FORMAT("{}: {} {}\n", particle_fit[i].name, particle_fit[i].fit_value, particle_fit[i].f(particle_fit[i].fit_value).real());
        }

        if( last_loss > loss && (last_loss - loss) / sample_points.size() < 1.e-2 ){
            std::cout << "last_loss: " << last_loss << "\n";
            std::cout << "loss: " << loss << "\n";
            std::cout << "ss: " << sample_points.size() << "\n";
            break;
        }
        last_loss = loss;
    }

//    std::ofstream out("plot.txt");


    scattering_amplitudes[0].scale([=](double){return particle_fit[0].f(particle_fit[0].fit_value).real();});
    Amplitude<1> interference = scattering_amplitudes[0];

    for( std::size_t i = 0; i < N_amplitudes; ++i ){
        std::cout << particle_fit[i].name << ": " << particle_fit[i].f(particle_fit[i].fit_value).real() << "\n";
        scattering_amplitudes[i].scale([=](double){return particle_fit[i].f(particle_fit[i].fit_value).real();});
        interference = interference + scattering_amplitudes[i];
    }

//    interference = interference + scattering_amplitudes[N_amplitudes-1];
    auto result = interference.scattering();

    std::cout << "interference={";
    for( auto& point : data_points ){
        std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
    }
    std::cout << "\b};";

    for( std::size_t i = 0; i < N_amplitudes; ++i ){
        std::cout << "\nresonance" << particle_fit[i].name << "={";
        result = scattering_amplitudes[i].scattering();
        for( auto& point : data_points ){
            std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
        }
        std::cout << "\b};";
    }

//    result = scattering_amplitudes.back().scattering();
//    std::cout << "\ntchannel={";
//    for( auto& point : data_points ){
//        std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
//    }
//    std::cout << "\b};\n";

    stopwatch.stop();
    std::cout << "Minimization done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
}