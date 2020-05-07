//#include "test_cases.h"
//#include <iostream>
//#include "Mesh.h"
//#include "vector_var.h"
//#include <stdio.h>      /* printf */
//#include <iostream>
//#include <math.h>
//#include "Boundary_Conditions.h"
//#include "Solution.h"
//#include "Solver.h"
//#include "quad_bcs.h"
//#include "external_forces.h"
//#include "global_variables.h"
//#include <algorithm>
//#include <string>
//#include <sstream>
//#include <boost/algorithm/string/replace.hpp>
//#include <boost/filesystem.hpp>
//#include <cstdio>
//#include <ctime>
//
//using namespace std;
//using namespace boost::filesystem;
//
//test_cases::test_cases()
//{
//    //ctor
//}
//
//test_cases::~test_cases()
//{
//    //dtor
//}
//
//void test_cases::west_to_east_poiseuille_flow(){
//
//    double X,Y,dx,dy,dt; // dt is streaming time step
//    double kine_viscosity,tau;
//    double reynolds, umax;
//    double simulation_length;
//    double delta_t; // time stepping step
//    quad_bcs_plus bcs;
//    double cs;
//    double pressure_grad;
//    global_variables globals;
//    std::string output_file;
//    double average_pressure;
//    double pre_condition_gamma;
//    std::clock_t start;
//    double duration;
//
//    pre_condition_gamma = 1;
//    //average_pressure =1*3*pre_condition_gamma;
//    average_pressure =1*3;
//
//    //vector_var pressure_gradient(-0.01,0,0), origin(1.1,0,0), origin_loc(0,0,0) ;
//    vector_var pressure_gradient(0,0,0), origin(average_pressure,0,0), origin_loc(0,0,0) ;
//
//    /// Parameters unique to test case
//
//    X= 16;
//    Y= 4.4;
//    dx= 0.4; // grid spacing
//    dy = 0.4;  // grid spacing
//    dt = 0.2;  // streaming time step -> dictates mach number -> grid spacing /2
//
//
//    /// Error :: let dt = l_dx = l_dy i.e. lattice spacing
//
//    simulation_length = 5000;
//    //kine_viscosity = U * X/ reynolds;
//    kine_viscosity = 0.0833333;
//
//    delta_t = 0.3;  // time marching step
//    cs = 1/sqrt(3);
//
//    reynolds = 33.7094;
//    umax = 0.17557;
//    // tau = 0.5 + Umax* Length / reynolds/dt  --- assumes c = 1.
//
//    tau = 0.5 + umax*X /reynolds /dt;
//    tau = 0.5 + umax*X /reynolds /dt *pre_condition_gamma;
//    //tau = 1/pre_condition_gammma*(tau-0.5) + 0.5;
//
//
//
//    pressure_grad =0;
//
//    //tau =0.75
//    output_file = create_output_directory(globals.tolerance,dx,delta_t);
//
//
//    // set boundary conditions for this test case
//    //bcs.w_rho = 1.1*3*pre_condition_gammma;
//    bcs.w_rho = 1.1*3*pre_condition_gamma;
//    bcs.w_u = 0;
//    bcs.w_v = 0;
//    bcs.w_w = 0;
//    bcs.w_type_vel = globals.neumann;
//    bcs.w_type_rho = globals.dirichlet;
//
//    bcs.s_rho = 0;
//    bcs.s_u = 0;
//    bcs.s_v = 0;
//    bcs.s_w = 0;
//    bcs.s_type_vel= globals.dirichlet;
//    bcs.s_type_rho = globals.neumann;
//
//    //bcs.e_rho = 1*3*pre_condition_gammma;
//    bcs.e_rho = 1*3*pre_condition_gamma;
//    bcs.e_u = 0;
//    bcs.e_v = 0;
//    bcs.e_w = 0;
//    bcs.e_type_vel= globals.neumann;
//    bcs.e_type_rho = globals.dirichlet;
//
//    bcs.n_rho = 0;
//    bcs.n_u = 0;
//    bcs.n_v = 0;
//    bcs.n_w = 0;
//    bcs.n_type_vel= globals.dirichlet;
//    bcs.n_type_rho = globals.neumann;
//
//
//
//
//
//    /// Methods to run the test case
//    //vector_var_tests();
//
//    // create Mesh
//    Mesh mesh(X,Y,dx,dy);
//
//    // create boundary conditions
//    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
//    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);
//
//    // assign external force terms
//    external_forces source_term(mesh.get_total_nodes());
//    source_term.set_uniform_force(pressure_grad);
//
//    //create solution
//    Solution soln(mesh.get_total_nodes());
//    soln.assign_pressure_gradient(pressure_gradient,origin_loc,origin,mesh);
//    soln.set_average_rho(average_pressure);
//    // Solvec
//
//    Solver solve;
//
//    solve.Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t,dx ,output_file,source_term,
//                              pre_condition_gamma);
//
//    soln.post_process(pre_condition_gamma);
//    soln.output(output_file);
//    tau = 1;
//
//}
//
//
//
//void test_cases::lid_driven_cavity_N(){
//
//    double X,Y,dx,dy,dt; // dt is streaming time step
//    double reynolds,kine_viscosity,tau;
//    double U;
//
//
//    quad_bcs bcs;
//    double cs,delta_t;
//
//    /// Parameters unique to test case
//    reynolds = 1000;
//    X= 100;
//    Y=100;
//    dx=1; // grid spacing
//    dy = 1;  // grid spacing
//    dt = 0.05;  // streaming time step let dt =dx = dy i.e. lattice spacing
//    U = 1;
//    simulation_length = 200;
//    kine_viscosity = U * X/ reynolds;
//    //kine_viscosity = 20;
//    delta_t = 0.1;
//    cs = 1/sqrt(3);
//    tau = kine_viscosity + 0.5* pow(cs,2) *dt;
//    //tau =0.75;
//
//
//    // set boundary conditions for this test case
//    bcs.w_rho = 1;
//    bcs.w_u = 0;
//    bcs.w_v = 0;
//    bcs.w_w = 0;
//
//    bcs.s_rho = 1;
//    bcs.s_u = 0;
//    bcs.s_v = 0;
//    bcs.s_w = 0;
//
//    bcs.e_rho = 1;
//    bcs.e_u = 0;
//    bcs.e_v = 0;
//    bcs.e_w = 0;
//
//    bcs.n_rho = 1;
//    bcs.n_u = 1;
//    bcs.n_v = 0;
//    bcs.n_w = 0;
//
//    /// Methods to run the test case
//    //vector_var_tests();
//
//    // create Mesh
//    Mesh mesh(X,Y,dx,dy);
//
//    // create boundary conditions
//    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
//    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);
//
//
//    //create solution
//    Solution soln(mesh.get_total_nodes());
//
//    // Solve
//
//    Solver solve;
//    //solve.Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t, dx,tolerance);
//
//    tau = 1;
//
//}
//void test_cases::west_to_east_couette_flow(){
//
//    double X,Y,dx,dy,dt; // dt is streaming time step
//    double kine_viscosity,tau;
//    double U;
//    double simulation_length;
//    double delta_t; // time stepping step
//    quad_bcs bcs;
//    double cs;
//    double pre_conditioned_gamma;
//    pre_conditioned_gamma = 1.0;
//    std::string output_file;
//
//    /// Parameters unique to test case
//
//    X= 0.3;
//    Y=1;
//    dx=0.1; // grid spacing
//    dy = 0.1;  // grid spacing
//    dt = 0.05;  // streaming time step -> dictates mach number -> grid spacing /2
//    /// Error :: let dt =dx = dy i.e. lattice spacing
//    U = 2;
//    simulation_length = 2500;
//    //kine_viscosity = U * X/ reynolds;
//    kine_viscosity = 0.6;
//
//    delta_t = 0.1;  // time marching step
//    cs = 1/sqrt(3);
//    tau = kine_viscosity + 0.5* pow(cs,2) *dt;
//
//
//    //output_file = "/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/";
//
//    output_file = "C:/Users/brendan/Dropbox/PhD/Test Cases/Couette Flow/";
//    //tau =0.75;
//
//    // set boundary conditions for this test case
//    bcs.w_rho = 1;
//    bcs.w_u = 1;
//    bcs.w_v = 0;
//    bcs.w_w = 0;
//    bcs.w_type = 3;
//
//    bcs.s_rho = 1;
//    bcs.s_u = 0;
//    bcs.s_v = 0;
//    bcs.s_w = 0;
//    bcs.s_type = 1;
//
//    bcs.e_rho = 0;
//    bcs.e_u = 0;
//    bcs.e_v = 0;
//    bcs.e_w = 0;
//    bcs.e_type =2;
//
//    bcs.n_rho = 1;
//    bcs.n_u = U;
//    bcs.n_v = 0;
//    bcs.n_w = 0;
//    bcs.n_type = 1;
//
//
//
//    /// Methods to run the test case
//    //vector_var_tests();
//
//    // create Mesh
//    Mesh mesh(X,Y,dx,dy);
//
//    // create boundary conditions
//    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
//    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);
//
//    // assign external force terms
//    external_forces source_term(mesh.get_total_nodes());
//    source_term.set_uniform_force(0.0);
//
//
//    //create solution
//    Solution soln(mesh.get_total_nodes());
//    soln.set_average_rho(1.0);
//
//    // Solve
//
//    Solver solve;
//    solve.Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t,dx ,output_file, source_term,pre_conditioned_gamma);
//    soln.output(output_file);
//    tau = 1;
//
//}
//
//void test_cases::east_to_west_couette_flow(){
//
//    double X,Y,dx,dy,dt; // dt is streaming time step
//    double kine_viscosity,tau;
//    double U;
//    double simulation_length;
//    double delta_t; // time stepping step
//    quad_bcs bcs;
//    double cs;
//
//    double pre_conditioned_gamma;
//    pre_conditioned_gamma = 1.0;
//
//    std::string output_file;
//     output_file = "/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/";
//    /// Parameters unique to test case
//
//    X= 0.3;
//    Y=1;
//    dx=0.1; // grid spacing
//    dy = 0.1;  // grid spacing
//    dt = 0.05;  // streaming time step -> dictates mach number -> grid spacing /2
//    /// Error :: let dt =dx = dy i.e. lattice spacing
//    U = -2;
//    simulation_length = 2500;
//    //kine_viscosity = U * X/ reynolds;
//    kine_viscosity = 0.6;
//
//    delta_t = 0.1;  // time marching step
//    cs = 1/sqrt(3);
//    tau = kine_viscosity + 0.5* pow(cs,2) *dt;
//    //tau =0.75;
//
//
//    output_file = "/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/";
//
//    // set boundary conditions for this test case
//    bcs.w_rho = 0;
//    bcs.w_u = 0;
//    bcs.w_v = 0;
//    bcs.w_w = 0;
//    bcs.w_type = 2;
//
//    bcs.s_rho = 1;
//    bcs.s_u = 0;
//    bcs.s_v = 0;
//    bcs.s_w = 0;
//    bcs.s_type = 1;
//
//    bcs.e_rho = 1;
//    bcs.e_u = 0;
//    bcs.e_v = 0;
//    bcs.e_w = 0;
//    bcs.e_type =3;
//
//    bcs.n_rho = 1;
//    bcs.n_u = U;
//    bcs.n_v = 0;
//    bcs.n_w = 0;
//    bcs.n_type = 1;
//
//
//
//    /// Methods to run the test case
//    //vector_var_tests();
//
//    // create Mesh
//    Mesh mesh(X,Y,dx,dy);
//
//    // create boundary conditions
//    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
//    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);
//
//    // assign external force terms
//    external_forces source_term(mesh.get_total_nodes());
//    source_term.set_uniform_force(0.0);
//
//    //create solution
//    Solution soln(mesh.get_total_nodes());
//
//    // Solve
//
//
//    Solver solve;
//    solve.Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t,dx, output_file,source_term,pre_conditioned_gamma);
//    soln.output(output_file);
//
//    tau = 1;
//
//}
//void test_cases::north_to_south_couette_flow(){
//
//    double X,Y,dx,dy,dt; // dt is streaming time step
//    double kine_viscosity,tau;
//    double U;
//    double simulation_length;
//    double delta_t; // time stepping step
//    quad_bcs bcs;
//    double cs;
//    std::string output_file;
//
//    double pre_conditioned_gamma;
//    pre_conditioned_gamma = 1.0;
//     output_file = "/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/";
//
//    /// Parameters unique to test case
//
//    X= 1;
//    Y=0.3;
//    dx=0.1; // grid spacing
//    dy = 0.1;  // grid spacing
//    dt = 0.05;  // streaming time step -> dictates mach number -> grid spacing /2
//    /// Error :: let dt =dx = dy i.e. lattice spacing
//    U = -2;
//    simulation_length = 2500;
//    //kine_viscosity = U * X/ reynolds;
//    kine_viscosity = 0.6;
//
//    delta_t = 0.1;  // time marching step
//    cs = 1/sqrt(3);
//    tau = kine_viscosity + 0.5* pow(cs,2) *dt;
//    //tau =0.75;
//
//
//
//    // set boundary conditions for this test case
//    bcs.w_rho = 1;
//    bcs.w_u = 0;
//    bcs.w_v = 0;
//    bcs.w_w = 0;
//    bcs.w_type = 1;
//
//    bcs.s_rho = 0;
//    bcs.s_u = 0;
//    bcs.s_v = 0;
//    bcs.s_w = 0;
//    bcs.s_type = 2;
//
//    bcs.e_rho = 1;
//    bcs.e_u = 0;
//    bcs.e_v = U;
//    bcs.e_w = 0;
//    bcs.e_type =1;
//
//    bcs.n_rho = 1;
//    bcs.n_u = 0;
//    bcs.n_v = 0;
//    bcs.n_w = 0;
//    bcs.n_type = 3;
//
//
//
//    /// Methods to run the test case
//    //vector_var_tests();
//
//    // create Mesh
//    Mesh mesh(X,Y,dx,dy);
//
//    // create boundary conditions
//    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
//    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);
//
//    // assign external force terms
//    external_forces source_term(mesh.get_total_nodes());
//    source_term.set_uniform_force(0.0);
//
//    //create solution
//    Solution soln(mesh.get_total_nodes());
//
//    // Solve
//
//    Solver solve;
//    solve.Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t,dx,output_file,source_term,pre_conditioned_gamma);
//
//    tau = 1;
//
//}
//
//void test_cases::south_to_north_couette_flow(){
//
//    double X,Y,dx,dy,dt; // dt is streaming time step
//    double kine_viscosity,tau;
//    double U;
//    double simulation_length;
//    double delta_t; // time stepping step
//    quad_bcs bcs;
//    double cs;
//    std::string output_file;
//
//    double pre_conditioned_gamma;
//    pre_conditioned_gamma = 1.0;
//     output_file = "/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/";
//
//    /// Parameters unique to test case
//
//    X= 1;
//    Y=0.3;
//    dx=0.1; // grid spacing
//    dy = 0.1;  // grid spacing
//    dt = 0.05;  // streaming time step -> dictates mach number -> grid spacing /2
//    /// Error :: let dt =dx = dy i.e. lattice spacing
//    U = 2;
//    simulation_length = 2500;
//    //kine_viscosity = U * X/ reynolds;
//    kine_viscosity = 0.6;
//
//    delta_t = 0.1;  // time marching step
//    cs = 1/sqrt(3);
//    tau = kine_viscosity + 0.5* pow(cs,2) *dt;
//    //tau =0.75;
//
//    // set boundary conditions for this test case
//    bcs.w_rho = 1;
//    bcs.w_u = 0;
//    bcs.w_v = 0;
//    bcs.w_w = 0;
//    bcs.w_type = 1;
//
//    bcs.s_rho = 1;
//    bcs.s_u = 0;
//    bcs.s_v = 0;
//    bcs.s_w = 0;
//    bcs.s_type = 3;
//
//    bcs.e_rho = 1;
//    bcs.e_u = 0;
//    bcs.e_v = U;
//    bcs.e_w = 0;
//    bcs.e_type =1;
//
//    bcs.n_rho = 0;
//    bcs.n_u = 0;
//    bcs.n_v = 0;
//    bcs.n_w = 0;
//    bcs.n_type = 2;
//
//
//
//    /// Methods to run the test case
//    //vector_var_tests();
//
//    // create Mesh
//    Mesh mesh(X,Y,dx,dy);
//
//    // create boundary conditions
//    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
//    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);
//
//    // assign external force terms
//    external_forces source_term(mesh.get_total_nodes());
//    source_term.set_uniform_force(0.0);
//
//    //create solution
//    Solution soln(mesh.get_total_nodes());
//
//    // Solve
//
//    Solver solve;
//    solve.Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t,dx, output_file, source_term,pre_conditioned_gamma);
//
//    tau = 1;
//
//}
//
//std::string test_cases::create_output_directory(double tol, double dx, double dt){
//
//    std::string output_file;
//    std::string folder;
//    std::ostringstream s;
//    //output_file = "C:/Users/brendan/Dropbox/PhD/Test Cases/Poiseuille Flow/";
//
//    output_file = "/home/brendan/Dropbox/PhD/Test Cases/Poiseuille Flow/";
//    s << "tol " << tol << " x " << dx
//     << " t " << dt;
//     folder = s.str();
//    //folder.replace(folder.begin(),folder.end(), ".",  "_");
//    boost::replace_all(folder, "." , "_");
//    output_file = output_file + folder;
//
//    boost::filesystem::path dir(output_file);
//    boost::filesystem::create_directories(dir);
//
//    return output_file;
//
//}
//
//void test_cases::vector_var_tests(){
//     vector_var a,b,c,d,e,f;
//
//    a.x = 3;
//    a.y = 4;
//    a.z = 0;
//
//    b.x = 2;
//    b.y = 2;
//    b.z = 2;
//
//    c.x = 4;
//    c.y = -3;
//    c.z = 0;
//
//    d.x= 0;
//    d.y = 0;
//    d.z = 0;
//
//    f.x = 5;
//    f.y = 5;
//    f.z = 0;
//
//
//
//
//    double mag , dp, ang;
//    mag = a.Magnitude();
//
//    cout <<  "Magnitude of a:" << mag << endl ;
//    dp = a.Dot_Product(b);
//    cout << "Dot product of a and b: " << dp << endl ;
//    ang = a.Angle_Between_Vectors(c);
//    ang = ang *360/2/M_PI;
//    cout << "angle between a and c: " << ang << endl ;
//    e.Get_Gradient(10,20,d,b);
//    cout << "gradient vector between d and b =" << e.x << "," << e.y << "," << e.z << endl ;
//
//
//
//}
