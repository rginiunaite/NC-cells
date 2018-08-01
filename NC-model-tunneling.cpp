/*
IBM model for NC cells coupled to reaction-diffusion model for chemoattractant, model described in SI McLennan et al. (2018). Followers are guided by the leaders by contact and tunnelling mechanism.
 
*********************  MATLAB TDR SYSTEM  ************************************
* Date created  : 2018, Jul 31
* Author(s)     : Rasa Giniunaite (giniunaite@maths.ox.ac.uk)
* Version       : 1.0
* Revisions     : 1.0 initial version (Rasa Giniunaite)
*
*********************  COPYRIGHT NOTICE  *************************************
* Copyright (C) 2018 Rasa Giniunaite
*                         University of Oxford
*                         United Kingdom
******************************************************************************
*/


#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> 
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Aboria;
using namespace Eigen; 

int main() {


    // model parameters

    /*
     * For computational efficiency we scale all the distance parameters by ten and time parameters are adapt to the timestep, which is one minute.
     */

    int n_seed = 1; // seeds are required when multiple simulations of the code are executed


    int length_x = 30; // initial domain length in x direction
    double domain_length = 30; // this variable is for the actual domain length, since it will be increasing
    double old_length = 30; // this variable will be used when updating the change in cell position due to domain growth
    const int length_y = 12; // initial domain length in y direction

    double cell_radius = 0.75;// // radius of a cell
    const double diameter = 2 * cell_radius; // diameter of a cell
    const int N_steps = 1440; // number of timesteps
    const size_t N = 5; // initial number of cells
    double l_filo_y = 2.75; // // sensing radius, filopodia + cell radius in y direction 
    double l_filo_x = 2.75; // sensing radius in x direction, it will have to be rescaled when domain grows
    double l_filo_x_in = l_filo_x; // this value is used for rescaling when domain grows based on initial value
    double l_filo_max = 4.5; // this is the length when two cells which were previously in a chain become dettached

    double diff_conc = 0.1; // sensing accuracy
    int insertion_freq = 1; // determines how frequently new cells are inserted, regulates the density of population

    double speed_l = 0.09; // speed of a leader cell, control case
    double increase_fol_speed = 1.3; // ratio leader over follower
    double speed_f = increase_fol_speed * speed_l; // speed of a follower cell

    double eps = 1; // distance a follower has to be ahead of a leader to swap phenotypes
    const int filo_number = 3; // the number of filopodia sent
    int same_dir = 0; // number of additional steps in the same direction, filopodia stabilisation
    bool random_pers = true; // persistent movement also when the cell moves randomly, polarity measure

    int count_dir = 0; // this is to count the number of times the cell moved the same direction, up to same_dir for each cell

    // tunnelling parameters
    double dist_thres = 1; // threshold distance to tunnel to enter it
    double track_spacing = 2; // spacing between positions on the track
    int track_length = 60; // number of leader positions tracked


    // domain growth parameters

    double L_0 = 30;
    double a = 0.288;
    double L_inf = 86.76;
    double t_s = 12.77; 
    double constant = 29.12 ;

    double domain_len_der = 0; // initialise derivative of the domain growth function


    // parameters for the dynamics of chemoattractant concentration

    double D = 0.0001; //diffusion coefficient
    double t = 0; // initialise time, redundant
    double dt = 1; // time step
    double dt_init = dt; // this is siginificant if dt < 1, for correct number of iterations in chemoattractant concentration
    int number_time = int(1/dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1; // space step in x direction
    double dy = 1; // space step in y direction
    double kai = 0.00001; // production rate of chemoattractant


    // parameters for internalisation

    double R = cell_radius; // cell radius
    double lam = 0.00035; // chemoattractant internalisation

/*
******************************************************************************
* initialise a matrix that stores values of concentration of chemoattractant
******************************************************************************
*/


    MatrixXf chemo = MatrixXf::Zero(length_x, length_y);
    MatrixXf chemo_new = MatrixXf::Zero(length_x, length_y);

    // initialise internalisation matrix
    MatrixXf intern = MatrixXf::Zero(length_x, length_y);

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1; // uniform concentration initially
            chemo_new(i, j) = 1; // this is for later updates
        }
    }


    // four columns for x, y, z, u (z is necessary for paraview, visualisation tool we use)

    // form a matrix which would store x,y,z,u

    MatrixXf chemo_3col(length_x * length_y, 4), chemo_3col_ind(length_x * length_y,
                                                                2); // need four because that is how paraview accepts data, third dimension is just zeros

    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < length_x * length_y) {
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                chemo_3col_ind(k, 0) = i;
                chemo_3col_ind(k, 1) = j;
                chemo_3col(k, 2) = 0;
                k += 1;
            }
        }
    }


    // y and x (initially) column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);
    }

    // create an array to keep track of all the positions of leader cells

    vdouble2 track_position [track_length][N] = {vdouble2(0,0)};
    for (int i =0; i < track_length; i++){
        for(int j = 0; j < N; j++){
            track_position[i][j] = vdouble2(0,0);
        }
    }
    int track_time[N] = {0}; // vector that stores the time values for each leader when there was a sufficiently big change in the position




/*
******************************************************************************
* initialise ABORIA variables
******************************************************************************
*/

	
    ABORIA_VARIABLE(radius, double, "radius") // cell radius
    ABORIA_VARIABLE(direction, vdouble2, "direction") // the direction a cell moved in the last step
    ABORIA_VARIABLE(persistence_extent, int, "persistence extent")// stores whether cell moves only one step in current direction or in a process of moving persistently
    ABORIA_VARIABLE(same_dir_step, int, "same dir step") // the number which stores how many steps in the same direction are made
    ABORIA_VARIABLE(attached_to_id, int, "attached_to_id") // ID of the cell that a follower is attached to
    ABORIA_VARIABLE(type, int, "type") // 0 if a cell is a leader, 1 if follower
    ABORIA_VARIABLE(chain_type, int, "chain_type") // leaders form different chain types, since we have fixed number of leaders, this is a simplification
    ABORIA_VARIABLE(chain, int, "chain") // stores whether a follower is part of the chain or not, 0 if it is not part of
    // the chain and then increasing integer if it is. If it is attached to a leader, it is 1, and then increasing order.

    // tunneling model
    ABORIA_VARIABLE(attached_at_time_step,int,"attached_at_time_step") // stores the timestep when leader was at particular tunnel position
    ABORIA_VARIABLE(attached_leader_nr,int,"attached_leader_nr") // stores leader ID
    ABORIA_VARIABLE(in_track, int, "in_track") // stores whether a follower is in a tunnel


    // create particle type with properties described above 
    typedef Particles<std::tuple<radius, type, attached_to_id, attached_at_time_step,in_track, attached_leader_nr, direction, chain, chain_type, persistence_extent, same_dir_step>, 2> particle_type; 

    // will use stored value of the position of a particle
    typedef particle_type::position position;

    // initialise the number of particles
    particle_type particles(N);

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(2, length_y - 1);


    /*
     * compact initialisation of particles
     */

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i]) = cell_radius;
        get<type>(particles[i]) = 0; // initially all cells are leaders


        get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // uniformly in y
        get<persistence_extent>(particles[i]) = 0;
        get<same_dir_step>(particles[i]) = 0;

    }

    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    vtkWriteGrid("particles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(t*n_seed); // choose different seeds to obtain different random numbers
    std::uniform_real_distribution<double> uniformpi(0, 2 * M_PI);


/*
******************************************************************************
* for each timestep
******************************************************************************
*/


    for (int t = 0; t < N_steps; t++) {


        // insert new cells at the start of the domain at insertion time 

        if (t % insertion_freq == 0) {
            bool free_position = false;
            particle_type::value_type f;

            get<position>(f) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y
            free_position = true;

            /*
             * loop over all neighbouring leaders within "diameter" distance
             */

            for (auto tpl = euclidean_search(particles.get_query(), get<position>(f), diameter); tpl != false; ++tpl) {

                vdouble2 diffx = tpl.dx();

                if (diffx.norm() < diameter) {
                    free_position = false;
                    break;
                }
            }

            // our assumption that all new cells are followers
            get<type>(f) = 1;


            if (free_position) {
                get<chain>(f) = 0;
                get<chain_type>(f) = -1;
                get<attached_to_id>(f) = -1;
                particles.push_back(f);
            }

        }

	// update particle positions
        particles.update_positions();


	/*
	******************************************************************************
	* Domain growth
	******************************************************************************
	*/


	// this loop is to make sure that I do correct number of timestep to solve chemoattractant equation
        for (int nt =0; nt < number_time; nt++){

		
            // update the length of the domain and its derivative
	    // we need to divide dt by 60 because the parameters for domain growth were inferred for hours

            domain_length = ((L_inf * exp(a * (dt/60.0 - t_s))) / (L_inf / L_0 + exp(a * (dt/60.0 - t_s)) - 1)) + constant;

            domain_len_der = ((a * L_inf/60 * exp(a * (dt/60.0 - t_s))) / (L_inf / L_0 + exp(a * (dt/60.0 - t_s)) - 1) -
                              (a * L_inf/60 * exp(2 * a * (dt/60.0 - t_s))) / ((L_inf / L_0 + exp(a * (dt/60.0 - t_s)) - 1)*(L_inf / L_0 + exp(a * (dt/60.0 - t_s)) - 1)));


            // chemoattractant internalisation

            for (int i = 0; i < length_x; i++) {
                for (int j = 0; j < length_y; j++) {
                    //go through all the cells
                    for (int k = 0; k < particles.size(); k++) {
                        vdouble2 x;
                        x = get<position>(particles[k]);
                        intern(i, j) = intern(i, j) + exp(-(((domain_length / length_x) * i - x[0]) *
                                                            ((domain_length / length_x) * i - x[0]) +
                                                            (j - x[1]) * (j - x[1])) /
                                                          (2 * R * R)); // mapping to fixed domain
                    }
                }
            }


            // reaction-diffusion equation for chemoattractant
            for (int i = 1; i < length_x - 1; i++) {
                for (int j = 1; j < length_y - 1; j++) {


                    // logistic production rate

                    chemo_new(i, j) = dt_init * (D * ((1 / ((domain_length / length_x) * (domain_length / length_x))) *
                                                      (chemo(i + 1, j) - 2 * chemo(i, j) + chemo(i - 1, j)) / (dx * dx) +
                                                      (chemo(i, j + 1) - 2 * chemo(i, j) + chemo(i, j - 1)) / (dy * dy)) -
                                                 (chemo(i, j) * lam / (2 * M_PI * R * R)) * intern(i, j) +
                                                 kai * chemo(i, j) * (1 - chemo(i, j)) -
                                                 double(domain_len_der) / double(domain_length) * chemo(i, j)) + chemo(i, j);

                }
            }



            // zero flux boundary conditions

            for (int i = 0; i < length_y; i++) {
                chemo_new(0, i) = chemo_new(1, i);
                chemo_new(length_x - 1, i) = chemo_new(length_x - 2, i);

            }

            for (int i = 0; i < length_x; i++) {
                chemo_new(i, 0) = chemo_new(i, 1);
                chemo_new(i, length_y - 1) = chemo_new(i, length_y - 2);
            }


            chemo = chemo_new; // update chemo concentration

	    // update timestep
            dt += dt_init;

        }



        // save the chemoattractant concentration with properly rescaled coordinates
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
        }


        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }


        // save data to plot chemoattractant concentration
        ofstream output("matrix_growing_domain" + to_string(t) + ".csv");

        output << "x, y, z, u" << "\n" << endl;


        for (int i = 0; i < length_x * length_y; i++) {
            for (int j = 0; j < 4; j++) {
                output << chemo_3col(i, j) << ", ";
            }
            output << "\n" << endl;
        }


        /// update positions uniformly based on the domain growth


            for (int i = 0; i < particles.size(); i++) {
                get<position>(particles)[i] *= vdouble2((domain_length / old_length), 1);
            }
            old_length = domain_length;
 

	/*
	******************************************************************************
         * Update the position of particles
	******************************************************************************
	*/


        //  create a random list of cell ids
        int check_rep = 0; // check for repetitions, 0 no rep, 1 rep


        std::default_random_engine gen2;
        gen2.seed(t*n_seed); // different seeds
        std::uniform_real_distribution<double> uniform_particles(0, particles.size()); // can only move forward

        VectorXi particle_id = VectorXi::Zero(particles.size());

        for (int i = 0; i < particles.size(); i++) {

            check_rep = 1; // set to 1 to enter the while loop
            while (check_rep == 1) {
                check_rep = 0; // it is initially zero and then will be changed to 1 if it is equivalent to others
                particle_id(i) = uniform_particles(gen2);


                for (int j = 0; j < i; j++) {
                    if (particle_id(i) == particle_id(j)) { check_rep = 1; }
                }
            }
        }

        // go through all the particles in a random order created above

        for (int j = 0; j < particles.size(); j++) {


		/*
		******************************************************************************
                * if a particle is a leader
		******************************************************************************
		*/

            if (get<type>(particles[particle_id(j)]) == 0) {

                vdouble2 x; // use variable x for the position of cells
                x = get<position>(particles[particle_id(j)]);

                double x_in; // x coordinate in initial domain length scale
                x_in = (length_x / domain_length) * x[0];
                l_filo_x = (length_x / domain_length) * l_filo_x_in; // rescale the length of filopodia as well



                // if it is still in the process of moving in the same direction
                if (get<persistence_extent>(particles[particle_id(j)]) == 1) {


                    x += get<direction>(particles[particle_id(j)]);

                    x_in = (length_x / domain_length) * x[0];// scale to initial coordinates

                    bool free_position = true; // check if the neighbouring position is free

                    // check if there are other particles in the position where the particle wants to move
                    for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {
                        if (get<id>(*k) != get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            free_position = false;
                        }
                    }

                    // check that the position they want to move to is free and not out of bounds
                    if (free_position && (x_in) > 0 && (x_in) < length_x - 1 && (x[1]) > 0 &&
                        (x[1]) < length_y - 1) {
                        // if that is the case, move into that position
                        get<position>(particles)[particle_id(j)] +=
                                get<direction>(particles)[particle_id(j)];
                    }
                    get<same_dir_step>(particles)[particle_id(
                            j)] += 1; // add regardless whether the step happened or no to that count of the number of movement in the same direction

                }


                // if a particle is not in a sequence of persistent steps
                if (get<persistence_extent>(particles[particle_id(j)]) == 0) {

                    // create an array to store random directions
                    array<double, filo_number + 1> random_angle;

                    // choose the number of angles where the filopodia is sent
                    for (int k = 0; k < filo_number + 1; k++) {

                        double random_angle_tem = uniformpi(gen1);
                        int sign_x_tem, sign_y_tem;

                        random_angle_tem = uniformpi(gen1);
                        random_angle[k] = random_angle_tem;

                    }


                    // choose which direction to move


                    // store variables for concentration at new locations


                    double old_chemo = chemo((round(x_in)), round(x)[1]);
                    array<double, filo_number> new_chemo;


                    for (int i = 0; i < filo_number; i++) {

                        if (round((x_in + sin(random_angle[i]) * l_filo_x)) < 0 ||
                            round((x_in + sin(random_angle[i]) * l_filo_x)) >
                            length_x - 1 || round(x[1] + cos(random_angle[i]) * l_filo_y) < 0 ||
                            round(x[1] + cos(random_angle[i]) * l_filo_y) > length_y - 1) {
                            new_chemo[i] = 0;
                        } else {
                            new_chemo[i] = chemo(round((x_in + sin(random_angle[i]) * l_filo_x)),
                                                 round(x[1] + cos(random_angle[i]) * l_filo_y));
                        }

                    }


                    // find maximum concentration of chemoattractant

                    int chemo_max_number = 0;

                    for (int i = 1; i < filo_number; i++) {
                        if (new_chemo[chemo_max_number] < new_chemo[i]) {
                            chemo_max_number = i;
                        }
                    }

                    // if the concentration in a new place is relatively higher than the old one, move that way

                    if ((new_chemo[chemo_max_number] - old_chemo) / sqrt(old_chemo) > diff_conc) {


                        count_dir += 1;

                        x += speed_l *
                             vdouble2(sin(random_angle[chemo_max_number]), cos(random_angle[chemo_max_number]));

                        x_in = (length_x / domain_length) * x[0];//scale to initial coordinates


                        bool free_position = true; // check if the neighbouring position is free

                        // check if the position the particle wants to move is free
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // if the position they want to move to is free and not out of bounds, move that direction
                        if (free_position && (x_in) > 0 && (x_in) < length_x - 1 && (x[1]) > 0 &&
                            (x[1]) < length_y - 1) {
                            get<position>(particles)[particle_id(j)] +=
                                    speed_l * vdouble2(sin(random_angle[chemo_max_number]),
                                                       cos(random_angle[chemo_max_number])); // update if nothing is in the next position
                            get<direction>(particles)[particle_id(j)] =
                                    speed_l * vdouble2(sin(random_angle[chemo_max_number]),
                                                       cos(random_angle[chemo_max_number]));

                            // if there is some kind of tendency to move persistently
                            if (same_dir > 0) {
                                get<persistence_extent>(particles[particle_id(
                                        j)]) = 1; // assume for now that it also becomes peristent in random direction

                            }

                        }


                    }

                        // if the concentration is not higher, move in random direction

                    else {

                        x += speed_l * vdouble2(sin(random_angle[filo_number]), cos(random_angle[filo_number]));

                        x_in = (length_x / domain_length) * x[0];// scale to initial length

                        bool free_position = true; // check if the neighbouring position is free

                        // if this loop is entered, it means that there is another cell where I want to move
                        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // update the position if the place they want to move to is free and not out of bounds
                        if (free_position && (x_in) > 0 && (x_in) < length_x - 1 && (x[1]) > 0 &&
                            (x[1]) < length_y - 1) {
                            get<position>(particles)[particle_id(j)] +=
                                    speed_l * vdouble2(sin(random_angle[filo_number]),
                                                       cos(random_angle[filo_number])); // update if nothing is in the next position
                            get<direction>(particles)[particle_id(j)] =
                                    speed_l * vdouble2(sin(random_angle[filo_number]),
                                                       cos(random_angle[filo_number]));
                            // if particles start moving persistently in all directions
                            if (random_pers) {
                                if (same_dir > 0) {
                                    get<persistence_extent>(particles[particle_id(
                                            j)]) = 1; // assume for now that it also becomes peristent in random direction

                                }
                            }

                        }

                    }

                }

                // check if it is not the end of moving in the same direction
                if (get<same_dir_step>(particles)[particle_id(j)] > same_dir) {
                    get<persistence_extent>(particles)[particle_id(j)] = 0;
                    get<same_dir_step>(particles[particle_id(j)]) = 0;
                }

            }

	    /*
	    ******************************************************************************
            * if a particle is a follower
            ******************************************************************************
	    */


            if (get<type>(particles[particle_id(j)]) == 1) {

                vdouble2 x;
                x = get<position>(particles[particle_id(j)]);

                double x_in; // x coordinate in initial domain length scale
                x_in = (length_x / domain_length) * x[0];

	    	/*
	    	******************************************************************************
            	* if a particle is part of the chain
            	******************************************************************************
	    	*/

                if (get<chain>(particles[particle_id(j)]) > 0) {


                    // check if it is not too far from the cell it was following

                    vdouble2 dist;

                    dist = get<position>(particles[particle_id(j)]) -
                           get<position>(particles[get<attached_to_id>(particles[particle_id(j)])]);

                    // if it is sufficiently far dettach the cell
                    if (dist.norm() > l_filo_max) {
                        get<chain>(particles[particle_id(j)]) = 0;
                        //dettach also all the cells that are behind it, so that other cells would not be attached to this chain
                        for (int i = 0; i < particles.size(); ++i) {
                            if (get<chain_type>(particles[i]) == get<chain_type>(particles)[particle_id(j)]) {

                                get<chain>(particles[i]) = 0;
                            }

                        }
                    }


                    // direction the same as of the cell it is attached to
                    get<direction>(particles)[particle_id(j)] = get<direction>(particles)[get<attached_to_id>(
                            particles[particle_id(j)])];

                    // try to move in the same direction as the cell it is attached to
                    vdouble2 x_chain = x + increase_fol_speed * get<direction>(particles)[particle_id(j)];

                    double x_in_chain; // scaled coordinate

                    x_in_chain =
                            (length_x / domain_length) * x_chain[0];//uniform growth in the first part of the domain


                    bool free_position = true;

                    // check if the position it wants to move to is free
                    for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter);
                         pos != false; ++pos) {
                        if (get<id>(*pos) !=
                            get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                            free_position = false;
                        }
                    }


                    // update the position if the place they want to move to is free and not out of bounds
                    if (free_position && (x_in_chain) > 0 && (x_in_chain) < length_x - 1 && (x_chain[1]) > 0 &&
                        (x_chain[1]) < length_y - 1) {
                        get<position>(particles)[particle_id(j)] +=
                                increase_fol_speed * get<direction>(particles[particle_id(j)]);

                    }

                }

                /*
                 * Model with tunnels
                 * */

	    	/*
	    	******************************************************************************
            	* if a particle is not part of the chain
            	******************************************************************************
	    	*/


                if (get<chain>(particles)[particle_id(j)] == 0) {

                    // try to find a close leader
                    for (auto k = euclidean_search(particles.get_query(), x, l_filo_x_in); k != false; ++k) {


                        if (get<type>(*k) == 0) { // if it is close to a leader
                            get<direction>(particles)[particle_id(j)] = get<direction>(*k); // set the same direction
                            get<chain>(particles)[particle_id(j)] = 1; // note that it is directly attached to a leader
                            get<attached_to_id>(particles)[particle_id(j)] = get<id>(
                                    *k); // note the id of the particle it is attached to
                            get<chain_type>(particles)[particle_id(j)] = get<id>(
                                    *k); // chain type is the id of the leader
                        }

                    }


                    // if it hasn't found a leader nearby,
                    // try to find a close follower which is in a chain contact with a leader

                    if (get<chain>(particles)[particle_id(j)] != 1) {
                        for (auto k = euclidean_search(particles.get_query(), x, l_filo_y); k != false; ++k) {

                            // if it is close to a follower that is part of the chain
                            if (get<type>(*k) == 1 && get<chain>(*k) > 0) {

                                if (get<id>(*k) != get<id>(particles[particle_id(j)])) {
                                  

                                    get<direction>(particles)[particle_id(j)] = get<direction>(*k);
                                    get<chain>(particles)[particle_id(j)] =
                                            get<chain>(*k) + 1; // it is subsequent member of the chain
                                    get<attached_to_id>(particles)[particle_id(j)] = get<id>(
                                            *k); // id of the particle it is attached to
                                    get<chain_type>(particles)[particle_id(j)] = get<chain_type>(*k); // chain type is
                                    // the same as the one of the particle it is attached to
                                }
                            }
                        }
                    }

                    // if it is in the chain
                    if (get<chain>(particles[particle_id(j)]) > 0){
                        //try to move in the same direction as the cell it is attached to
                        vdouble2 x_chain = x + increase_fol_speed * get<direction>(particles)[particle_id(j)];

                        // Non-uniform domain growth
                        double x_in_chain;

                        x_in_chain =
                                (length_x / domain_length) * x_chain[0];//uniform growth in the first part of the domain



                        bool free_position = true;


                        // check if the position it wants to move is free
                        for (auto pos = euclidean_search(particles.get_query(), x_chain, diameter); pos != false; ++pos) {

                            if (get<id>(*pos) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }


                        // if the position is free and not out of bounds, move that direction
                        if (free_position &&
                            (x_in_chain) > 0 &&
                            (x_in_chain) < length_x - 1 && (x_chain[1]) > 0 &&
                            (x_chain[1]) < length_y - 1) {
                            get<position>(particles)[particle_id(j)] +=
                                    increase_fol_speed * get<direction>(particles[particle_id(j)]);

                        }

                    }


		        /*	
		    	******************************************************************************
		    	* if a cell has not find close cells in chains, try to find a tunnel or move randomly
		    	******************************************************************************
		    	*/

                    if (get<chain>(particles[particle_id(j)]) == 0) {


                        x = get<position>(particles[particle_id(j)]);


                        /*
                         * Check if there is a tunnel nearby
                         * */

                        // check what the closest neighbour is

                        vdouble2 diff; // difference in position
                        double previous_diff = dist_thres; // initialise distance threshold


                        for (int k = 0; k < N; k++) {
                            for (int i = 0; i < track_length; i++) {
                                diff = x - track_position[i][k];
                                if (diff.norm() < dist_thres) {
                                    if (diff.norm() < previous_diff) {
                                        for (int f = 0;
                                             f <
                                             particles.size(); f++) { // make sure other particles are not in that position


                                            if (get<attached_at_time_step>(particles[f]) != j &&
                                                get<attached_leader_nr>(particles[f]) != k) {
                                                get<attached_at_time_step>(particles[particle_id(j)]) = i;
                                                get<attached_leader_nr>(particles[particle_id(j)]) = k;

                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // move in the direction where track goes

                        vdouble2 direc_follow = track_position[get<attached_at_time_step>(
                                particles[particle_id(j)])+1][get<attached_leader_nr>(
                                particles[particle_id(j)])] - track_position[get<attached_at_time_step>(
                                particles[particle_id(j)])][get<attached_leader_nr>(
                                particles[particle_id(j)])];

                        double normalise = direc_follow.norm();
                        //normalise
                        direc_follow = direc_follow/normalise;
                        vdouble2 x_can = x + direc_follow * speed_f;



                        bool free_position = true; // check if the neighbouring position is free

                        // if this loop is entered, it means that there is another cell where I want to move
                        for (auto k = euclidean_search(particles.get_query(), x_can, diameter); k != false; ++k) {

                            if (get<id>(*k) !=
                                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                free_position = false;
                            }
                        }



                        // if the cell has found a track, update its position

                        if (free_position &&
                            ((x_can[0] * (length_x / domain_length))) > 0 &&
                            ((x_can[0] * (length_x / domain_length))) < length_x - 1 && (x_can[1]) > 0 &&
                            (x_can[1]) < length_y - 1) {

                            get<position>(particles[particle_id(j)]) = x_can;
                            get<attached_at_time_step>(particles[particle_id(j)]) += 1;
                            get<in_track>(particles[particle_id(j)]) = 1;
                        }
                        else{
                            get<in_track>(particles[particle_id(j)]) = 0;
                        }



                        // if a cell hasn't found a track, move randomly

                        if (get<in_track>(particles[particle_id(j)]) == 0) {

                            x = get<position>(particles)[particle_id(j)]; // make sure that x is set correctly


                            double random_angle = uniformpi(gen1);


                            while (round((x_in + sin(random_angle) * l_filo_x)) < 0 ||
                                   round((x_in + sin(random_angle) * l_filo_x)) >
                                   length_x - 1 || round(x[1] + cos(random_angle) * l_filo_y) < 0 ||
                                   round(x[1] + cos(random_angle) * l_filo_y) > length_y - 1) {

                                random_angle = uniformpi(gen1);


                            }

                            x += speed_f * vdouble2(sin(random_angle), cos(random_angle));
                          

                            x_in = (length_x / domain_length) * x[0];


                            bool free_position = true; // check if the neighbouring position is free

                            // if this loop is entered, it means that there is another cell where I want to mov

                            for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {



                                //for (int i=0; i < particles.size(); i++) {
                                if (get<id>(*k) !=
                                    get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                                    free_position = false;
                                    break;
                                }
                            }



                            // check that the position they want to move to is free and not out of bounds
                            if (free_position && (x_in) > 0 &&
                                (x_in) < length_x - 1 && (x[1]) > 0 &&
                                (x[1]) < length_y - 1) {
                                get<position>(particles)[particle_id(j)] += speed_f * vdouble2(sin(random_angle),
                                                                                               cos(random_angle)); // update if nothing is in the next position
                            }
                        }
                    }


                }


		/*
		******************************************************************************
                * phenotype switching
		******************************************************************************
		*/


                // find the closest leader


                // so that I would not go through all the cells I will choose the ones that are closer to the front
                // minimum position in x of the leaders

                int min_index = 0;

                for (int i = 1; i < N; ++i) {
                    if (get<position>(particles[i])[0] < get<position>(particles[min_index])[0]) {
                        min_index = i;
                    }

                }

                // if a follower is eps further in front than the leader, swap their types
                if (get<position>(particles[particle_id(j)])[0] > get<position>(particles[min_index])[0] + eps) {
                    // find distance to all the leaders
                    double distances[N];
                    vdouble2 dist_vector;
                    //check which one is the closest
                    for (int i = 0; i < N; ++i) {
                        dist_vector = get<position>(particles[particle_id(j)]) - get<position>(particles[i]);
                        distances[i] = dist_vector.norm();

                        int winning_index = 0;
                        for (int i = 1; i < N; ++i) {
                            if (distances[i] < distances[winning_index]) {
                                winning_index = i;
                            }
                        }

                        // if this closest leader is behind that follower, swap them
                        if (get<position>(particles[particle_id(j)])[0] >
                            get<position>(particles[winning_index])[0] + eps) {
                            particle_type::value_type tmp = particles[winning_index];


                            // their position swap

                            vdouble2 temp = get<position>(particles[winning_index]);
                            get<position>(particles[winning_index]) = get<position>(particles[particle_id(j)]);
                            get<position>(particles[particle_id(j)]) = temp;


                        }

                    }


                }

            }

        }

        // update positions
        particles.update_positions();


        // save at every time step
#ifdef HAVE_VTK
        vtkWriteGrid("particles", t, particles.get_grid(true));
#endif


    }

}

