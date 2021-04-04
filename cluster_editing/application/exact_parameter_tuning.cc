
#include <iostream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/exact/star_bound.h>

#include <cluster_editing/utils/timer.h>
#include <boost/program_options.hpp>

using namespace std;


int main(int argc, char *argv[]) {
    namespace po = boost::program_options;

    int num;
    int max_num_unchanged = global_star_bound_config.max_num_unchanged;
    float probability = global_star_bound_config.min_degree_candidate_vs_random_candidate_probability;
    int time_limit = 10;

    po::options_description options("Options");
    options.add_options()
            ("help", "show help message")
            ("instance", po::value<int>(&num)->required())
            ("max_num_unchanged", po::value<int>(&max_num_unchanged)->default_value(max_num_unchanged))
            ("probability", po::value<float>(&probability)->default_value(probability))
            ("time_limit", po::value<int>(&time_limit)->default_value(time_limit));

    po::variables_map cmd_vm;
    po::store(po::parse_command_line(argc, argv, options), cmd_vm);
    if (cmd_vm.count("help") != 0 || argc == 1) {
        std::cout << options << std::endl;
        exit(0);
    }

    po::notify(cmd_vm);

    global_star_bound_config.max_num_unchanged = max_num_unchanged;
    global_star_bound_config.min_degree_candidate_vs_random_candidate_probability = probability;

    auto inst = load_exact_instance(num);

    auto &timer = cluster_editing::utils::Timer::instance();
    timer.start_timer("solving", "solving instance");

    ExactSolver solver;
    solver.verbose = false;
    solver.time_limit = chrono::steady_clock::now() + chrono::seconds(time_limit);
    auto solution = solver.solve(inst);

    timer.stop_timer("solving");
    cout << std::boolalpha << solution.worked << "\t" << timer.get("solving") << endl;

    return 0;
}
