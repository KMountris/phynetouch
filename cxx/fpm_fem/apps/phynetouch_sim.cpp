/*
 * PHYNETOUCH. RF ablation simulation software.
 * Copyright (C) 2021  Konstantinos A. Mountris <konstantinos.mountris@gmail.com>
 * ALL RIGHTS RESERVED
 */


#include "PHYNETOUCH/PHYNETOUCH"

#include "tools/parser.hpp"
#include "tools/config_sim.hpp"

#include <termcolor/termcolor.hpp>
#include <nlohmann/json.hpp>

#include <petscsys.h>
#include <mpi.h>

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>
#include <algorithm>

using json = nlohmann::json;
using namespace PNT;
using namespace PNTSIM;

int main(int argc, char *argv[]) {

    // Initialize MPI and PETSc.
    MPI_Init(&argc, &argv);
    PetscCall(PetscInitialize(&argc, &argv, 0, ""));

    try {

        MpiHandler mpi_handler;
        PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_handler.rank_id));
        PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_handler.procs_num));

        if (mpi_handler.rank_id == 0) {
            std::cout << termcolor::bold << "\nWelcome to: \n";
            std::cout << termcolor::green << "  _____  _           _   _   _______               _    " << termcolor::magenta << "  _____ _ \n";
            std::cout << termcolor::green << " |  __ \\| |         | \\ | | |__   __|             | |   " << termcolor::magenta << " / ____(_) \n";
            std::cout << termcolor::green << " | |__) | |__  _   _|  \\| | ___| | ___  _   _  ___| |__ " << termcolor::magenta << "| (___  _ _ __ ___  \n";
            std::cout << termcolor::green << " |  ___/| '_ \\| | | | . ` |/ _ \\ |/ _ \\| | | |/ __| '_ \\ " << termcolor::magenta << "\\___ \\| | '_ ` _ \\ \n";
            std::cout << termcolor::green << " | |    | | | | |_| | |\\  |  __/ | (_) | |_| | (__| | | |" << termcolor::magenta << "____) | | | | | | | \n";
            std::cout << termcolor::green << " |_|    |_| |_|\\__, |_| \\_|\\___|_|\\___/ \\__,_|\\___|_| |_|" << termcolor::magenta << "_____/|_|_| |_| |_| \n";
            std::cout << termcolor::green << "                __/ | \n";
            std::cout << termcolor::green << "               |___/  \n";
            std::cout << "\n" << termcolor::reset;

            std::cout << termcolor::bold << "version: " << termcolor::reset << PHYNETOUCH_VERSION << "\n";
            std::cout << termcolor::bold << "licence: " << termcolor::reset << "To be defined\n";
            std::cout << termcolor::bold << "author: " << termcolor::reset << "Konstantinos A. Mountris\n";
            std::cout << termcolor::bold << "email: " << termcolor::reset << "konstantinos.mountris@gmail.com\n";
            std::cout << termcolor::bold << "web: " << termcolor::reset << "https://www.mountris.org\n";
            std::cout << "\n";
        }


        // Get configuration file.
        std::string pnt_filename = "";
        if (argc == 1) {
            std::cout << termcolor::yellow << Logger::Warning("Configuration file was not provided during launching the PhyNeTouchSim app.\n"
                "                  Give path to configuration filename.") << termcolor::reset;
            // Read configuration file given by the user.
            std::cout << "Path: ";
            std::cin >> pnt_filename;
        }
        else { pnt_filename = argv[1]; }

        // Parse input json file.
        Parser sim_parser(pnt_filename);

        // Set up simulation for given scale.
        short dims = sim_parser.GetValue<short>("tissue.geometry.dimensions");
        short cell_nodes = sim_parser.GetValue<short>("tissue.geometry.cell vertices");

        // if (dims == 1 && cell_nodes == 2) {
        //     ConfigSim<1, 2> config_sim; config_sim.Launch(sim_parser, std::cout);
        // } else if (dims == 2 && cell_nodes == 3) {
        //     ConfigSim<2, 3> config_sim; config_sim.Launch(sim_parser, std::cout);
        // } else if (dims == 2 && cell_nodes == 4) {
        //     ConfigSim<2, 4> config_sim; config_sim.Launch(sim_parser, std::cout);
        // } else if (dims == 3 && cell_nodes == 4) {
        //     ConfigSim<3, 4> config_sim; config_sim.Launch(sim_parser, std::cout);
        // } else if (dims == 3 && cell_nodes == 8) {
        //     ConfigSim<3, 8> config_sim; config_sim.Launch(sim_parser, std::cout);
        // }

        if (dims == 2 && cell_nodes == 3) {
            ConfigSim<2, 3> config_sim; config_sim.Launch(sim_parser, mpi_handler, std::cout);
        } else if (dims == 3 && cell_nodes == 4) {
            ConfigSim<3, 4> config_sim; config_sim.Launch(sim_parser, mpi_handler, std::cout);
        }

        if (mpi_handler.rank_id == 0) {
            std::cout << termcolor::magenta << termcolor::bold;
            std::cout << Logger::Message("The simulation finished successfully. Thank you for using the PhyNeTouchSim app.\n") << termcolor::reset;
        }

    } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (const std::runtime_error &e) {
        std::cerr << "Runtime error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (const std::out_of_range &e) {
        std::cerr << "Out of Range error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (const std::bad_alloc &e) {
        std::cerr << "Bad allocation error:" << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (...) {
        std::cerr << termcolor::red << "PhyNeTouch unknown exception. Check for typing errors in the provided *.json file..." << termcolor::reset << std::endl;
    }

    // Finalize MPI and PETSc.
    PetscFinalize();
    MPI_Finalize();

    return EXIT_SUCCESS;
}