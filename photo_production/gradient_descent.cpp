#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"


#include <feynumeric/feynumeric.hpp>
#include <feynumeric/utility.hpp>

#include <fmt/chrono.h>

#include <chrono>
#include <iostream>
#include <random>
#include <feynumeric/phase_space.hpp>

double sigmoid(double x){ return 1./(1. + std::exp(-x));}

Feynumeric::func_t<1> interval(double min, double max){
    return [min, max](double x){ return min + (max-min) * sigmoid(x);};
}
double interval_isigmoid(double min, double max, double x){
    return std::log( (max-x)/(x-min) );
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
//	cmd.register_command("start", std::string("1.0"), "starting point");
//	cmd.register_command("end", std::string("3.0"), "end value");
//	cmd.register_command("steps", std::string("200"), "steps");
//	cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING, Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ, Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));
    cmd.crash_on_missing_mandatory_command();

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Pi_Plus = P.get("pi+");
    init_vertices(P, cmd.as_string("coupling_constants"));

    std::vector<Polynomial> s_poly;
    std::vector<Polynomial> c_poly;

    struct Resonance{
        std::string name;
        func_t<1> sigm;
        double start;
        Resonance(std::string const& str, double min, double max, std::mt19937& gen)
                : name(str)
                , sigm(interval(min, max))
                , start(interval_isigmoid(min, max, std::uniform_real_distribution<double>(min, max)(gen)))
        {
        }
    };

    std::random_device r;
//	std::default_random_engine eng{r()};
    std::mt19937 random_generator{r()};

    std::cout << "yes\n";

//	std::vector<std::string> resonances = {"D1232","D1600", "D1620", "D1700", "D1750", "D1900", "D1905", "D1910", "D1920", "D1930", "D1940", "D1950"};
    std::vector<Resonance> resonances{
            {"D1232", 0.99399999, 0.99400001, random_generator},
            {"D1600", 0.08, 0.24, random_generator},
            {"D1620", 0.25, 0.35, random_generator},
            {"D1700", 0.1, 0.2, random_generator},
            {"D1750", 0, 1, random_generator},
            {"D1900", 0.04, 0.12, random_generator},
            {"D1905", 0.09, 0.15, random_generator},
            {"D1910", 0.1, 0.3, random_generator},
            {"D1920", 0.05, 0.2, random_generator},
            {"D1930", 0.05, 0.15, random_generator},
            {"D1940", 0.01, 0.07, random_generator},
            {"D1950", 0.35, 0.45, random_generator},
    };

    if( cmd.exists("branching_ratios") ){
        std::ifstream infile(cmd.as_string("branching_ratios"));
        if( !infile ){
            std::cerr << "Could not open branching_ratio file.\n" << cmd.as_string("branching_ratios") << "\n";
        }
        double a, b;
        int i = 0;
        while (infile >> a >> b)
        {
            resonances[i].sigm = interval(a, b);
            resonances[i].start = interval_isigmoid(a,  b, std::uniform_real_distribution<double>(a, b)(random_generator));
            i++;
        }
    }


    std::uniform_real_distribution<double> dist(0, 1);

    auto const N_EPOCHS = cmd.as_int("n_epochs");
    auto const rate = cmd.as_double("rate");

    std::cout << FORMAT("n_epochs: {} rate: {}\n", N_EPOCHS, rate);

    auto inverse_sigmoid = [](double x){ return std::log(x/(1.-x));};

    std::vector<double> sigm_branching_ratios(resonances.size());

    for( std::size_t i = 0; i < resonances.size(); ++i ){
        sigm_branching_ratios[i] = resonances[i].start;
        std::cout << FORMAT("{}: {} {}\n", resonances[i].name, sigm_branching_ratios[i], resonances[i].sigm(sigm_branching_ratios[i]).real());
    }

    std::vector<Amplitude<0>> width_amplitudes;
    std::vector<Amplitude<1>> scattering_amplitudes;
    std::vector<std::pair<Amplitude<1>,Amplitude<1>>> scattering_amplitudes_derivatives;

    Timer stopwatch;
    stopwatch.start();

//	for( std::size_t k = 0; k < resonances.size(); ++k){
    for( auto const& resonance : resonances ){
        auto particle = P.get(FORMAT("{}pp", resonance.name));
        /* width */

        auto n_spin_states = particle->spin().n_states() * Proton->spin().n_states();
        std::vector<std::vector<Polynomial>> width_polynomials(1);
        for( std::size_t i = 0; i < n_spin_states; ++i ){
            Polynomial temp;
            temp.load(FORMAT("data/polynomials/polynomial_decay_{}_0_{}.txt", particle->name(), i));
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

        func_t<1> breit_wigner = [=](double sqrt_s) mutable {
            return 1. / ( sqrt_s * sqrt_s - particle->mass() * particle->mass() +
                          1.i * sqrt_s * M_width.width(sqrt_s) );
        };

        Amplitude<1> M_scattering(scattering_polynomials, {Proton, Pi_Plus}, {Proton, Pi_Plus});

        func_t<1> form_factor = [=](double sqrt_s) mutable{
            auto const lambda = 0.8;
            auto const l4 = std::pow(lambda, 4);
            auto const delta = sqrt_s * sqrt_s - particle->mass() * particle->mass();
            return l4 / ( delta * delta + l4 );
        };

        M_scattering.scale(breit_wigner);
        M_scattering.scale(form_factor);
        scattering_amplitudes.push_back(M_scattering);

    }

    if( scattering_amplitudes.empty() ){
        error("No particle selected.");
        return EXIT_SUCCESS;
    }

    auto particle = P.get("rho0");

    auto n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
    std::vector<std::vector<Polynomial>> scattering_polynomials(2);
    for( std::size_t i = 0; i < n_spin_states; ++i ){
        Polynomial temp_s, temp_c;
        temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
        temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
        scattering_polynomials[0].push_back(temp_s);
        scattering_polynomials[1].push_back(temp_c);
    }
    Amplitude<1> M_scattering(scattering_polynomials, {Proton, Pi_Plus}, {Proton, Pi_Plus});

    scattering_amplitudes.push_back(M_scattering);

    stopwatch.stop();
    std::cout << "config done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
    stopwatch.start();

    std::string file = cmd.as_string("data_file");

    std::ifstream in(file);

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

//    struct Point2D{
//        double x, y;
//    };
//
//    std::vector<Point2D> data_points = {{1.105, 6.15}, {1.115, 9.8}, {1.116, 12.}, {1.126, 17.6}, {1.127, 15.8}, {1.133, 20.4}, {1.137, 19.}, {1.138, 26.09}, {1.139, 25.3}, {1.144, 31.}, {1.149, 37.3}, {1.158, 50.3}, {1.159, 53.74}, {1.168, 71.4}, {1.169, 75.2}, {1.17, 77.}, {1.176, 91.}, {1.177, 96.25}, {1.178, 99.32}, {1.18, 103.7}, {1.184, 122.}, {1.189, 127.733}, {1.191, 140.395}, {1.193, 140.8}, {1.195, 150.}, {1.196, 140.5}, {1.197, 151.}, {1.198, 150.}, {1.199, 167.7}, {1.201, 173.2}, {1.202, 157.75}, {1.203, 175.}, {1.206, 179.83}, {1.207, 162.}, {1.21, 181.25}, {1.211, 198.3}, {1.213, 178.6}, {1.214, 181.5}, {1.215, 204.24}, {1.217, 202.58}, {1.219, 199.25}, {1.22, 198.35}, {1.221, 198.85}, {1.222, 198.}, {1.226, 185.303}, {1.227, 197.55}, {1.228, 192.}, {1.23, 194.2}, {1.231, 194.1}, {1.232, 188.15}, {1.235, 200.}, {1.236, 174.}, {1.237, 195.15}, {1.24, 177.9}, {1.241, 178.3}, {1.242, 174.7}, {1.244, 175.73}, {1.246, 179.}, {1.248, 148.}, {1.25, 141.}, {1.252, 154.167}, {1.254, 156.}, {1.255, 140.9}, {1.257, 148.}, {1.261, 133.603}, {1.263, 131.9}, {1.27, 127.2}, {1.272, 113.9}, {1.273, 111.9}, {1.275, 114.5}, {1.28, 100.63}, {1.282, 96.2}, {1.283, 96.6}, {1.286, 111.}, {1.287, 107.}, {1.292, 82.8}, {1.293, 84.1}, {1.299, 88.}, {1.301, 73.82}, {1.303, 73.8}, {1.312, 69.655}, {1.313, 64.7}, {1.318, 65.7}, {1.32, 60.3}, {1.322, 57.4}, {1.326, 57.39}, {1.332, 51.}, {1.338, 53.}, {1.339, 46.63}, {1.342, 46.9}, {1.345, 41.23}, {1.347, 45.22}, {1.352, 41.105}, {1.361, 40.69}, {1.362, 37.4733}, {1.366, 38.155}, {1.377, 31.15}, {1.379, 31.9}, {1.381, 34.58}, {1.39, 31.4}, {1.392, 29.02}, {1.395, 28.075}, {1.405, 26.25}, {1.408, 26.53}, {1.411, 24.36}, {1.416, 24.8}, {1.417, 26.6}, {1.418, 23.33}, {1.421, 29.44}, {1.422, 23.79}, {1.429, 23.47}, {1.43, 21.45}, {1.431, 21.3}, {1.434, 16.06}, {1.435, 22.23}, {1.443, 20.625}, {1.447, 23.8}, {1.45, 18.95}, {1.453, 24.14}, {1.454, 18.74}, {1.465, 18.635}, {1.468, 19.99}, {1.47, 18.1}, {1.477, 15.96}, {1.478, 16.1}, {1.486, 14.93}, {1.487, 15.47}, {1.488, 14.95}, {1.492, 17.26}, {1.499, 14.46}, {1.5, 16.2}, {1.503, 15.13}, {1.51, 18.6}, {1.511, 14.12}, {1.512, 15.88}, {1.515, 15.}, {1.521, 14.1}, {1.525, 16.28}, {1.529, 15.11}, {1.533, 14.35}, {1.535, 14.71}, {1.544, 16.285}, {1.547, 16.1}, {1.551, 14.805}, {1.553, 15.21}, {1.555, 15.47}, {1.557, 14.5}, {1.564, 15.75}, {1.566, 15.41}, {1.567, 15.96}, {1.571, 18.}, {1.572, 16.6}, {1.573, 17.}, {1.576, 16.42}, {1.579, 16.89}, {1.585, 17.9}, {1.587, 17.5}, {1.597, 18.6}, {1.598, 18.64}, {1.6, 18.66}, {1.603, 21.58}, {1.607, 20.7}, {1.612, 20.97}, {1.613, 20.38}, {1.615, 20.06}, {1.627, 19.5}, {1.629, 23.}, {1.633, 22.29}, {1.635, 24.98}, {1.641, 23.585}, {1.642, 21.4}, {1.643, 23.5}, {1.65, 23.38}, {1.652, 24.4}, {1.654, 22.49}, {1.663, 25.1}, {1.667, 25.6}, {1.668, 21.88}, {1.669, 24.33}, {1.672, 24.43}, {1.679, 25.35}, {1.683, 25.185}, {1.686, 25.9}, {1.692, 27.04}, {1.694, 26.005}, {1.697, 23.12}, {1.699, 25.29}, {1.701, 26.5}, {1.717, 26.98}, {1.718, 26.5}, {1.721, 25.97}, {1.724, 23.8}, {1.727, 26.28}, {1.739, 26.5}, {1.742, 28.9}, {1.743, 25.325}, {1.751, 25.85}, {1.754, 28.05}, {1.761, 27.8}, {1.764, 28.57}, {1.765, 29.69}, {1.772, 29.}, {1.774, 30.4}, {1.777, 26.78}, {1.78, 29.56}, {1.791, 30.5}, {1.793, 31.}, {1.795, 31.45}, {1.796, 31.8}, {1.804, 27.51}, {1.807, 32.56}, {1.81, 33.2}, {1.82, 35.12}, {1.821, 31.3}, {1.824, 33.7}, {1.826, 35.}, {1.829, 30.65}, {1.831, 35.85}, {1.832, 35.93}, {1.839, 36.76}, {1.844, 36.1}, {1.85, 38.33}, {1.851, 38.3}, {1.854, 35.27}, {1.865, 38.1}, {1.868, 39.53}, {1.871, 40.05}, {1.872, 38.8}, {1.874, 39.9}, {1.878, 40.5}, {1.88, 36.64}, {1.881, 39.4}, {1.883, 41.11}, {1.888, 41.25}, {1.891, 39.5}, {1.892, 37.37}, {1.894, 40.6}, {1.898, 41.46}, {1.904, 38.07}, {1.908, 41.55}, {1.911, 39.1}, {1.915, 41.56}, {1.916, 37.65}, {1.918, 40.9}, {1.919, 41.02}, {1.925, 41.4}, {1.929, 36.16}, {1.932, 42.52}, {1.933, 40.25}, {1.935, 41.4}, {1.939, 40.08}, {1.943, 39.5}, {1.952, 36.53}, {1.957, 38.32}, {1.962, 38.2}, {1.966, 39.53}, {1.968, 38.}, {1.976, 34.1}, {1.978, 38.25}, {1.981, 36.05}, {1.99, 36.1}, {1.992, 35.3}, {1.997, 34.9}, {1.999, 39.27}, {2.004, 34.57}, {2.009, 34.4}, {2.02, 33.85}, {2.035, 32.2}, {2.036, 32.7}, {2.039, 30.3}, {2.055, 31.5}, {2.064, 31.07}, {2.066, 31.4}, {2.071, 32.6}};
    std::vector<Data> sample_points(data_points.begin() + 59, data_points.begin()+260);

    std::size_t const N_amplitudes = scattering_amplitudes.size();

    double const epsilon = 0.001;

    Timer epoch_timer;
    for( std::size_t epoch = 0; epoch < N_EPOCHS; ++epoch ){
        epoch_timer.start();
        double loss = 0.;

        std::vector<Amplitude<1>> scattering_left(scattering_amplitudes.begin(), scattering_amplitudes.end());
        std::vector<Amplitude<1>> scattering_mid(scattering_amplitudes.begin(), scattering_amplitudes.end());
        std::vector<Amplitude<1>> scattering_right(scattering_amplitudes.begin(), scattering_amplitudes.end());

        for( std::size_t i = 1; i < N_amplitudes - 1; ++i ){
            auto scale_left = sigmoid(sigm_branching_ratios[i] - epsilon);
            auto scale_mid = sigmoid(sigm_branching_ratios[i]);
            auto scale_right = sigmoid(sigm_branching_ratios[i] + epsilon);
            scattering_left[i].scale([scale_left](double){ return scale_left; });
            scattering_mid[i].scale([scale_mid ](double){ return scale_mid;});
            scattering_right[i].scale([scale_right](double){ return scale_right; });
        }

        func_t<1> y_hat;
        std::vector<func_t<1>> y_hat0(N_amplitudes);
        std::vector<func_t<1>> y_hat1(N_amplitudes);


        Amplitude<1> interference = scattering_mid[0];
        for( std::size_t i = 1; i < N_amplitudes - 1; ++i ){
            interference = interference + scattering_mid[i];
            Amplitude<1> interference_left = scattering_left[0];
            Amplitude<1> interference_right = scattering_right[0];
            for( std::size_t j = 1; j < N_amplitudes - 1; ++j ){
                if( i == j ){
                    interference_left = interference_left + scattering_left[j];
                    interference_right = interference_right + scattering_right[j];
                } else{
                    interference_left = interference_left + scattering_mid[j];
                    interference_right = interference_right + scattering_mid[j];
                }
            }
            interference_left = interference_left + scattering_mid[N_amplitudes - 1];
            interference_right = interference_right + scattering_mid[N_amplitudes - 1];
            y_hat0[i] = interference_left.scattering();
            y_hat1[i] = interference_right.scattering();
        }
        interference = interference + scattering_mid[N_amplitudes - 1];
        y_hat = interference.scattering();


        std::vector<double> derivatives(N_amplitudes, 0.);


        for( auto const& point : sample_points ){
            double diff = point.cross - y_hat(point.srt).real() * static_cast<double>(1._2mbarn);
            loss += diff*diff / (point.dcros * point.dcros);
            for( std::size_t k = 1; k < N_amplitudes-1; ++k ){
                double const loss0 = std::pow(point.cross - y_hat0[k](point.srt).real() * static_cast<double>(1._2mbarn), 2) / (point.dcros * point.dcros);
                double const loss1 = std::pow(point.cross - y_hat1[k](point.srt).real() * static_cast<double>(1._2mbarn), 2) / (point.dcros * point.dcros);
                double const d = (loss1 - loss0) / (2. * epsilon);
                derivatives[k] += d;
            }
        }

        for( std::size_t k = 1; k < N_amplitudes-1; ++k ){
            double temp = rate * derivatives[k] / sample_points.size();
            sigm_branching_ratios[k] -= std::abs(temp) > 0.2? sgn(temp) * 0.2 : temp;
        }

        epoch_timer.stop();
        std::cout << "Epoch: " << epoch << " [" << epoch_timer.time<std::chrono::milliseconds>()/1000. << "] loss: [" << loss / sample_points.size() << "]\n";

        for( std::size_t i = 0; i < resonances.size(); ++i ){
            std::cout << FORMAT("{}: {} {}\n", resonances[i].name, sigm_branching_ratios[i], sigmoid(sigm_branching_ratios[i]));
        }

        if( loss / sample_points.size() < 1.e-3 ){
            break;
        }
    }

//    std::ofstream out("plot.txt");

    Amplitude<1> interference = scattering_amplitudes[0];
    for( std::size_t i = 1; i < N_amplitudes - 1; ++i ){
        std::cout << i << ": " << resonances[i].sigm(sigm_branching_ratios[i]).real() << "\n";
        scattering_amplitudes[i].scale([=](double){return resonances[i].sigm(sigm_branching_ratios[i]).real();});
        interference = interference + scattering_amplitudes[i];
    }

    interference = interference + scattering_amplitudes[N_amplitudes-1];
    auto result = interference.scattering();

    std::cout << "interference={";
    for( auto& point : data_points ){
        std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
    }
    std::cout << "\b};";

    for( std::size_t i = 0; i < resonances.size(); ++i ){
        std::cout << "\nresonance" << resonances[i].name << "={";
        result = scattering_amplitudes[i].scattering();
        for( auto& point : data_points ){
            std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
        }
        std::cout << "\b};";
    }

    result = scattering_amplitudes.back().scattering();
    std::cout << "\ntchannel={";
    for( auto& point : data_points ){
        std::cout << "{" << point.srt << "," << FORMAT("{:.9f}", result(point.srt).real() * static_cast<double>(1._2mbarn)) << "},";
    }
    std::cout << "\b};\n";

    stopwatch.stop();
    std::cout << "Minimization done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
}